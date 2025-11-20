/*! @file WriteArray.c
    @author Debojyoti Ghosh
    @brief Write a vector field, stored as an array, to file

    Contains functions to write out a vector field, stored
    as an array to a file.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mpivars.h>
#include <hypar.h>

/* Function declarations */
static int WriteArraySerial   (int,int,int*,int*,int,double*,double*,void*,void*,char*);
#ifndef serial
static int WriteArrayParallel (int,int,int*,int*,int,double*,double*,void*,void*,char*);
#endif

/*! Write out a vector field, stored as an array, to file: wrapper function that calls
    the appropriate function depending on output mode (#HyPar::m_output_mode). The
    output file format is determined by #HyPar::m_op_file_format
*/
int WriteArray(
                int     a_ndims,        /*!< Number of spatial dimensions */
                int     a_nvars,        /*!< Number of variables per grid point */
                int     *a_dim_global,  /*!< Integer array of size ndims with global grid size in each dimension */
                int     *a_dim_local,   /*!< Integer array of size ndims with local  grid size in each dimension */
                int     a_ghosts,       /*!< Number of ghost points */
                double  *a_x,           /*!< Array of spatial coordinates (i.e. the grid) */
                double  *a_u,           /*!< Vector field to write */
                void    *a_s,           /*!< Solver object of type #HyPar */
                void    *a_m,           /*!< MPI object of type #MPIVariables */
                char*   a_fname_root    /*!< Filename root (extension is added automatically). For unsteady output,
                                           a numerical index is added that is the same as for the solution output files. */
              )
{
  HyPar         *solver = (HyPar*)       a_s;
  MPIVariables  *mpi    = (MPIVariables*)a_m;
  _DECLARE_IERR_;

  /* if WriteOutput() is NULL, then return */
  if (!solver->WriteOutput) return(0);

#ifndef serial
  if (!strcmp(solver->m_output_mode,"serial")) {
#endif
    IERR WriteArraySerial(a_ndims,a_nvars,a_dim_global,a_dim_local,a_ghosts,a_x,a_u,
                          solver,mpi,a_fname_root); CHECKERR(ierr);
#ifndef serial
  } else {
    IERR WriteArrayParallel(a_ndims,a_nvars,a_dim_global,a_dim_local,a_ghosts,a_x,a_u,
                            solver,mpi,a_fname_root); CHECKERR(ierr);
  }
#endif

  return(0);
}

/*!
  Function to write out a vector field, stored as an array, and
  its associated Cartesian grid to a file in serial mode. It will
  allocate the global domain on rank 0, so do not use for big
  problems for which the entire global domain will not fit on one
  node. This approach is also not very scalable.
  + Supports binary and ASCII formats (specified through
    #HyPar::m_op_file_format).

  \sa WriteBinary(), WriteText(), WriteTecplot2D(), WriteTecplot3D()
*/
int WriteArraySerial(
                      int     a_ndims,        /*!< Number of spatial dimensions */
                      int     a_nvars,        /*!< Number of variables per grid point */
                      int     *a_dim_global,  /*!< Integer array of size ndims with global grid size in each dimension */
                      int     *a_dim_local,   /*!< Integer array of size ndims with local  grid size in each dimension */
                      int     a_ghosts,       /*!< Number of ghost points */
                      double  *a_x,           /*!< Array of spatial coordinates (i.e. the grid) */
                      double  *a_u,           /*!< Vector field to write */
                      void    *a_s,           /*!< Solver object of type #HyPar */
                      void    *a_m,           /*!< MPI object of type #MPIVariables */
                      char*   a_fname_root    /*!< Filename root (extension is added automatically). For unsteady output,
                                                 a numerical index is added that is the same as for the solution output files. */
                    )
{
  HyPar         *solver = (HyPar*)       a_s;
  MPIVariables  *mpi    = (MPIVariables*)a_m;
  int           d;
  _DECLARE_IERR_;

  /* root process: allocate global output arrays */
  double *ug, *xg;
  if (!mpi->m_rank) {
    int size_global;

    size_global = 1;
    for (d=0; d<a_ndims; d++) size_global *= a_dim_global[d];
    ug = (double*) calloc (size_global*a_nvars,sizeof(double));
    _ArraySetValue_(ug,size_global*a_nvars,0.0);

    size_global = 0;
    for (d=0; d<a_ndims; d++) size_global += a_dim_global[d];
    xg = (double*) calloc (size_global,sizeof(double));
    _ArraySetValue_(xg,size_global,0.0); CHECKERR(ierr);

  } else {

    /* null pointers on non-root processes */
    ug = xg = NULL;

  }

  /* Assemble the local output arrays into the global output arrays */
  IERR MPIGatherArraynD(a_ndims,mpi,ug,a_u,a_dim_global,a_dim_local,
                        a_ghosts,a_nvars);  CHECKERR(ierr);
  int offset_global, offset_local;
  offset_global = offset_local = 0;
  for (d=0; d<a_ndims; d++) {
    IERR MPIGatherArray1D(mpi,(mpi->m_rank?NULL:&xg[offset_global]),
                            &a_x[offset_local+a_ghosts],
                            mpi->m_is[d],mpi->m_ie[d],a_dim_local[d],0); CHECKERR(ierr);
    offset_global += a_dim_global[d];
    offset_local  += a_dim_local [d] + 2*a_ghosts;
  }

  if (!mpi->m_rank) {
    /* write output file to disk */
    char filename[_MAX_STRING_SIZE_] = "";
    strcat(filename,a_fname_root);
    if (!strcmp(solver->m_op_overwrite,"no")) {
      strcat(filename,"_");
      strcat(filename,solver->m_filename_index);
    }
    strcat(filename,solver->m_solnfilename_extn);
    printf("Writing solution file %a_s.\n",filename);
    IERR solver->WriteOutput(a_ndims,a_nvars,a_dim_global,xg,ug,filename,
                             solver->m_index); CHECKERR(ierr);

    /* Clean up output arrays */
    free(xg);
    free(ug);
  }

  return(0);
}

/*! Write a vector field, stored as an array, and its associated Cartesian grid
    to a file in parallel. All the MPI ranks are divided into IO groups,
    each with a group leader that writes out the data for the ranks in
    its group. The number of groups (and thus, the number of ranks
    participating in file I/O) is given by #MPIVariables::N_IORanks.
    The group leader receives the local data from each rank in its
    group, and writes it out to the corresponding file.
    + The data is written in binary format only.
    + The number of files written is equal to the number of IO groups
      (#MPIVariables::N_IORanks), and are named as <fname_root>.bin.nnnn
      where "nnnn" is a string of format "%04d" representing n,
      0 <= n < MPIVariables::N_IORanks.
    + Each file contains the following blocks of data:\n
      {\n
        x0_i (0 <= i < dim_local[0])\n
        x1_i (0 <= i < dim_local[1])\n
        ...\n
        x{ndims-1}_i (0 <= i < dim_local[ndims-1])\n
        [u0,u1,...,u{nvars-1}]_p (0 <= p < N) (with no commas)\n
        \n
        where \n
        x0, x1, ..., x{ndims-1} represent the spatial dimensions (for a 3D problem, x0 = x, x1 = y, x2 = z),\n
        u0, u1, ..., u{nvars-1} are each component of the vector u at a grid point,\n
        N = dim_local[0]*dim_local[1]*...*dim_local[ndims-1] is the total number of points,\n
        and p = i0 + dim_local[0]*( i1 + dim_local[1]*( i2 + dim_local[2]*( ... _ i{ndims-1} )))\n
        with i0, i1, i2, etc representing grid indices along each spatial dimension, i.e.,\n
        0 <= i0 < dim_local[0]-1\n
        0 <= i1 < dim_local[1]-1\n
        ...\n
        0 <= i{ndims-1} < dim_local[ndims=1]-1\n
      }\n
      for each rank in the corresponding group (i.e., there are #MPIVariables::nproc divided by
      #MPIVariables::N_IORanks such blocks in each file).
    + To stitch all the local data in these files into the global solution, and write that out to
      a binary file can be done by Extras/ParallelOutput.c.
    + If #HyPar::m_op_overwrite is set to 0, the vector field at the various simulation times are appended
      to each of the <fname_root>.bin.nnnn. The code Extras/ParallelOutput.c will take care of writing
      the global solution at each simulation time to a different file (in binary format) (the same files
      that WriteArraySerial() would have written out if serial file output was chosen).

    This approach has been observed to be very scalable (with up to ~500,000 MPI ranks). The number of MPI
    ranks participating in file I/O (#MPIVariables::N_IORanks) should be set to the number of I/O nodes
    available on a HPC platform, given the number of compute nodes the simulation is running on.
*/
#ifndef serial
int WriteArrayParallel(
                        int     a_ndims,        /*!< Number of spatial dimensions */
                        int     a_nvars,        /*!< Number of variables per grid point */
                        int     *a_dim_global,  /*!< Integer array of size ndims with global grid size in each dimension */
                        int     *a_dim_local,   /*!< Integer array of size ndims with local  grid size in each dimension */
                        int     a_ghosts,       /*!< Number of ghost points */
                        double  *a_x,           /*!< Array of spatial coordinates (i.e. the grid) */
                        double  *a_u,           /*!< Vector field to write */
                        void    *a_s,           /*!< Solver object of type #HyPar */
                        void    *a_m,           /*!< MPI object of type #MPIVariables */
                        char*   a_fname_root    /*!< Filename root (extension is added automatically). For unsteady output,
                                                   a numerical index is added that is the same as for the solution output files. */
                      )
{
  HyPar         *solver = (HyPar*)        a_s;
  MPIVariables  *mpi    = (MPIVariables*) a_m;
  int           proc,d;
  _DECLARE_IERR_;

  static int count = 0;

  char filename_root[_MAX_STRING_SIZE_];
  strcpy(filename_root,a_fname_root);
  strcat(filename_root,solver->m_solnfilename_extn);
  if (!mpi->m_rank) printf("Writing solution file %a_s.xxxx (parallel mode).\n",filename_root);

  /* calculate size of the local grid on this rank */
  int sizex = 0;     for (d=0; d<a_ndims; d++) sizex += a_dim_local[d];
  int sizeu = a_nvars; for (d=0; d<a_ndims; d++) sizeu *= a_dim_local[d];

  /* allocate buffer arrays to write grid and solution */
  double *buffer = (double*) calloc (sizex+sizeu, sizeof(double));

  /* copy the grid to buffer */
  int offset1 = 0, offset2 = 0;
  for (d = 0; d < a_ndims; d++) {
    _ArrayCopy1D_((a_x+offset1+a_ghosts),(buffer+offset2),a_dim_local[d]);
    offset1 += (a_dim_local[d]+2*a_ghosts);
    offset2 +=  a_dim_local[d];
  }

  /* copy the solution */
  int index[a_ndims];
  IERR ArrayCopynD(a_ndims,a_u,(buffer+sizex),a_dim_local,a_ghosts,0,index,a_nvars); CHECKERR(ierr);

  if (mpi->m_IOParticipant) {

    /* if this rank is responsible for file I/O */
    double *write_buffer = NULL;
    int     write_size_x, write_size_u, write_total_size;
    int     is[a_ndims], ie[a_ndims];

    /* open the file */
    FILE *out;
    int  bytes;
    char filename[_MAX_STRING_SIZE_];
    MPIGetFilename(filename_root,&mpi->m_IOWorld,filename);

    if (!strcmp(solver->m_op_overwrite,"no")) {
      if ((!count) && (!solver->m_restart_iter)) {
        /* open a new file, since this function is being called the first time
           and this is not a restart run*/
        out = fopen(filename,"wb");
        if (!out) {
          fprintf(stderr,"Error in WriteArrayParallel(): File %a_s could not be opened for writing.\n",filename);
          return(1);
        }
      } else {
        /* append to existing file */
        out = fopen(filename,"ab");
        if (!out) {
          fprintf(stderr,"Error in WriteArrayParallel(): File %a_s could not be opened for appending.\n",filename);
          return(1);
        }
      }
    } else {
      /* write a new file / overwrite existing file */
      out = fopen(filename,"wb");
      if (!out) {
        fprintf(stderr,"Error in WriteArrayParallel(): File %a_s could not be opened for writing.\n",filename);
        return(1);
      }
    }
    count++;

    /* Write own data and free buffer */
    bytes = fwrite(buffer,sizeof(double),(sizex+sizeu),out);
    if (bytes != (sizex+sizeu)) {
      fprintf(stderr,"Error in WriteArrayParallel(): Failed to write data to file %a_s.\n",filename);
      return(1);
    }
    free(buffer);

    /* receive and write the data for the other processors in this IO rank'a_s group */
    for (proc=mpi->m_GroupStartRank+1; proc<mpi->m_GroupEndRank; proc++) {
      /* get the local domain limits for process proc */
      IERR MPILocalDomainLimits(a_ndims,proc,mpi,a_dim_global,is,ie);
      /* calculate the size of its local data and allocate write buffer */
      write_size_x = 0;      for (d=0; d<a_ndims; d++) write_size_x += (ie[d]-is[d]);
      write_size_u = a_nvars;  for (d=0; d<a_ndims; d++) write_size_u *= (ie[d]-is[d]);
      write_total_size = write_size_x + write_size_u;
      write_buffer = (double*) calloc (write_total_size, sizeof(double));
      /* receive the data */
      MPI_Request req = MPI_REQUEST_NULL;
      MPI_Irecv(write_buffer,write_total_size,MPI_DOUBLE,proc,1449,mpi->m_world,&req);
      MPI_Wait(&req,MPI_STATUS_IGNORE);
      /* write the data */
      bytes = fwrite(write_buffer,sizeof(double),write_total_size,out);
      if (bytes != write_total_size) {
        fprintf(stderr,"Error in WriteArrayParallel(): Failed to write data to file %a_s.\n",filename);
        return(1);
      }
      free(write_buffer);
    }

    /* close the file */
    fclose(out);

  } else {

    /* all other processes, just send the data to the rank responsible for file I/O */
    MPI_Request req = MPI_REQUEST_NULL;
    MPI_Isend(buffer,(sizex+sizeu),MPI_DOUBLE,mpi->m_IORank,1449,mpi->m_world,&req);
    MPI_Wait(&req,MPI_STATUS_IGNORE);
    free(buffer);

  }

  return(0);
}
#endif

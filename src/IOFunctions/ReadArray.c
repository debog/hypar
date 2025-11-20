/*! @file ReadArray.c
    @author Debojyoti Ghosh
    @brief Read in a vector field from file
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mpivars.h>
#include <hypar.h>

static int ReadArraySerial    (int,int,int*,int*,int,void*,void*,double*,double*,char*,int*);
#ifndef serial
static int ReadArrayParallel  (int,int,int*,int*,int,void*,void*,double*,double*,char*,int*);
static int ReadArrayMPI_IO    (int,int,int*,int*,int,void*,void*,double*,double*,char*,int*);
#endif

/*! Read in a vector field from file: wrapper function that calls
    the appropriate function depending on input mode (#HyPar::m_input_mode).\n\n
    The mode and type of input are specified through #HyPar::m_input_mode and
    #HyPar::m_ip_file_type. A vector field is read from file and stored in an array.
*/
int ReadArray(
              int     a_ndims,        /*!< Number of spatial dimensions */
              int     a_nvars,        /*!< Number of variables per grid point */
              int     *a_dim_global,  /*!< Integer array of size ndims with global grid size in each dimension */
              int     *a_dim_local,   /*!< Integer array of size ndims with local  grid size in each dimension */
              int     a_ghosts,       /*!< Number of ghost points */
              void    *a_s,           /*!< Solver object of type #HyPar */
              void    *a_m,           /*!< MPI object of type #MPIVariables */
              double  *a_x,           /*!< Grid associated with the array (can be NULL) */
              double  *a_u,           /*!< Array to hold the vector field */
              char    *a_fname_root,  /*!< Filename root */
              int     *a_read_flag    /*!< Flag to indicate if the file was read */
             )
{
  HyPar         *solver = (HyPar*) a_s;
  MPIVariables  *mpi    = (MPIVariables*) a_m;
  _DECLARE_IERR_;

  if      (!strcmp(solver->m_input_mode,"serial")) {
    IERR ReadArraySerial(a_ndims,a_nvars,a_dim_global,a_dim_local,a_ghosts,a_s,a_m,a_x,a_u,a_fname_root,a_read_flag);
    CHECKERR(ierr);
#ifndef serial
  } else if (!strcmp(solver->m_input_mode,"parallel")) {
    ReadArrayParallel(a_ndims,a_nvars,a_dim_global,a_dim_local,a_ghosts,a_s,a_m,a_x,a_u,a_fname_root,a_read_flag);
    CHECKERR(ierr);
  } else if (!strcmp(solver->m_input_mode,"mpi-io"  )) {
    ReadArrayMPI_IO(a_ndims,a_nvars,a_dim_global,a_dim_local,a_ghosts,a_s,a_m,a_x,a_u,a_fname_root,a_read_flag);
    CHECKERR(ierr);
#endif
  } else {
    fprintf(stderr,"Error: Illegal value (%a_s) for input_mode.\n",solver->m_input_mode);
    return(1);
  }

  if (a_x) {
    int offset, d;
    /* exchange MPI-boundary values of a_x between processors */
    offset = 0;
    for (d = 0; d < a_ndims; d++) {
      IERR MPIExchangeBoundaries1D(mpi,&a_x[offset],a_dim_local[d],
                                   a_ghosts,d,a_ndims); CHECKERR(ierr);
      offset  += (a_dim_local [d] + 2*a_ghosts);
    }
    /* fill in ghost values of a_x at physical boundaries by extrapolation */
    offset = 0;
    for (d = 0; d < a_ndims; d++) {
      double *X     = &a_x[offset];
      int    *dim   = a_dim_local, i;
      if (mpi->m_ip[d] == 0) {
        /* fill left boundary along this dimension */
        for (i = 0; i < a_ghosts; i++) {
          int delta = a_ghosts - i;
          X[i] = X[a_ghosts] + ((double) delta) * (X[a_ghosts]-X[a_ghosts+1]);
        }
      }
      if (mpi->m_ip[d] == mpi->m_iproc[d]-1) {
        /* fill right boundary along this dimension */
        for (i = dim[d]+a_ghosts; i < dim[d]+2*a_ghosts; i++) {
          int delta = i - (dim[d]+a_ghosts-1);
          X[i] =  X[dim[d]+a_ghosts-1]
                  + ((double) delta) * (X[dim[d]+a_ghosts-1]-X[dim[d]+a_ghosts-2]);
        }
      }
      offset  += (dim[d] + 2*a_ghosts);
    }
  }

  return(0);
}

/*! Read an array in a serial fashion: For multi-processor simulation, only rank 0
    reads in the entire solution from the file, and then distributes the relevant portions
    to each of the processors. This involves memory allocation for the global domain on rank
    0. Thus, do not use for large domains. This approach is not very scalable either, if running
    with a very large number of processors (> 1000). Supports both binary and ASCII formats.
    \n\n
    The name of the file being read is <fname_root>.inp
    \n\n
    \b ASCII format:-\n
    The input file should contain the ASCII data as follows:\n
    \n
    x0_i (0 <= i < dim_global[0])\n
    x1_i (0 <= i < dim_global[1])\n
    ...\n
    x{ndims-1}_i (0 <= i < dim_global[ndims-1])\n
    u0_p (0 <= p < N)\n
    u1_p (0 <= p < N)\n
    ...\n
    u{nvars-1}_p (0 <= p < N)\n
    \n
    where \n
    x0, x1, ..., x{ndims-1} represent the spatial dimensions (for a 3D problem, x0 = x, x1 = y, x2 = z),\n
    u0, u1, ..., u{nvars-1} are each component of the vector u,\n
    N = dim_global[0]*dim_global[1]*...*dim_global[ndims-1] is the total number of points,\n
    and p = i0 + dim_global[0]*( i1 + dim_global[1]*( i2 + dim_global[2]*( ... + dim_global[ndims-2]*i{ndims-1} )))
    (see #_ArrayIndexnD_)\n
    with i0, i1, i2, etc representing grid indices along each spatial dimension, i.e.,\n
    0 <= i0 < dim_global[0]-1\n
    0 <= i1 < dim_global[1]-1\n
    ...\n
    0 <= i{ndims-1} < dim_global[ndims=1]-1\n
    \n\n
    \b Binary format:-\n
    The input file should contain the binary data in as follows:\n
    \n
    x0_i (0 <= i < dim_global[0])\n
    x1_i (0 <= i < dim_global[1])\n
    ...\n
    x{ndims-1}_i (0 <= i < dim_global[ndims-1])\n
    [u0,u1,...,u{nvars-1}]_p (0 <= p < N) (with no commas)\n
    \n
    where \n
    x0, x1, ..., x{ndims-1} represent the spatial dimensions (for a 3D problem, x0 = x, x1 = y, x2 = z),\n
    u0, u1, ..., u{nvars-1} are each component of the vector u at a grid point,\n
    N = dim_global[0]*dim_global[1]*...*dim_global[ndims-1] is the total number of points,\n
    and p = i0 + dim_global[0]*( i1 + dim_global[1]*( i2 + dim_global[2]*( ... + dim_global[ndims-2]*i{ndims-1} )))
    (see #_ArrayIndexnD_)\n
    with i0, i1, i2, etc representing grid indices along each spatial dimension, i.e.,\n
    0 <= i0 < dim_global[0]-1\n
    0 <= i1 < dim_global[1]-1\n
    ...\n
    0 <= i{ndims-1} < dim_global[ndims=1]-1\n
    \n
    For serial runs, this is the only input mode (of course!).
*/
int ReadArraySerial(
                      int     a_ndims,        /*!< Number of spatial dimensions */
                      int     a_nvars,        /*!< Number of variables per grid point */
                      int     *a_dim_global,  /*!< Integer array of size ndims with global grid size in each dimension */
                      int     *a_dim_local,   /*!< Integer array of size ndims with local  grid size in each dimension */
                      int     a_ghosts,       /*!< Number of ghost points */
                      void    *a_s,           /*!< Solver object of type #HyPar */
                      void    *a_m,           /*!< MPI object of type #MPIVariables */
                      double  *a_x,           /*!< Grid associated with the array (can be NULL) */
                      double  *a_u,           /*!< Array to hold the vector field being read */
                      char    *a_fname_root,  /*!< Filename root */
                      int     *a_read_flag    /*!< Flag to indicate if the file was read */
                    )
{
  HyPar         *solver = (HyPar*)        a_s;
  MPIVariables  *mpi    = (MPIVariables*) a_m;
  int           i, d, ferr, index[a_ndims];
  double        *ug = NULL, *xg = NULL;
  _DECLARE_IERR_;

  *a_read_flag = 0;
  /* Only root process reads from the file */
  if (!mpi->m_rank) {

    if (!strcmp(solver->m_ip_file_type,"ascii")) {
      char filename[_MAX_STRING_SIZE_];
      strcpy(filename,a_fname_root);
      strcat(filename,".inp");
      FILE *in; in = fopen(filename,"r");
      if (!in) *a_read_flag = 0;
      else {
        *a_read_flag = 1;
        /* Reading from file */
        printf("Reading array from ASCII file %a_s (Serial mode).\n",filename);
        int size,offset;
        /* allocate global solution array */
        size  = 1; for (d=0; d<a_ndims; d++) size *= a_dim_global[d]; size *= a_nvars;
        ug    = (double*) calloc(size,sizeof(double));
        size  = 0; for (d=0; d<a_ndims; d++) size += a_dim_global[d];
        xg    = (double*) calloc(size,sizeof(double));

        /* read grid */
        offset = 0;
        for (d = 0; d < a_ndims; d++) {
          for (i = 0; i < a_dim_global[d]; i++) {
            ferr = fscanf(in,"%lf",&xg[i+offset]);
            if (ferr != 1) {
              printf("Error in ReadArraySerial(): unable to read data. ferr=%d\n", ferr);
              exit(1);
            }
          }
          offset += a_dim_global[d];
        }

        /* read solution */
        for (i = 0; i < a_nvars; i++) {
          int done = 0; _ArraySetValue_(index,a_ndims,0);
          while (!done) {
            int p; _ArrayIndex1D_(a_ndims,a_dim_global,index,0,p);
            ferr = fscanf(in,"%lf",&ug[p*a_nvars+i]);
            if (ferr != 1) {
              printf("Error in ReadArraySerial(): unable to read data. ferr=%d\n", ferr);
              exit(1);
            }
            _ArrayIncrementIndex_(a_ndims,a_dim_global,index,done);
          }
        }

        fclose(in);
      }
    } else if ((!strcmp(solver->m_ip_file_type,"bin")) || (!strcmp(solver->m_ip_file_type,"binary"))) {

      char filename[_MAX_STRING_SIZE_];
      strcpy(filename,a_fname_root);
      strcat(filename,".inp");
      FILE *in; in = fopen(filename,"rb");
      if (!in) *a_read_flag = 0;
      else {
        *a_read_flag = 1;
        printf("Reading array from binary file %a_s (Serial mode).\n",filename);
        size_t bytes;
        int size;
        /* allocate global solution array */
        size  = 1; for (d=0; d<a_ndims; d++) size *= a_dim_global[d]; size *= a_nvars;
        ug = (double*) calloc(size,sizeof(double));
        size = 0; for (d=0; d<a_ndims; d++) size += a_dim_global[d];
        xg      = (double*) calloc(size,sizeof(double));

        /* read grid */
        size = 0;
        for (d = 0; d < a_ndims; d++) size += a_dim_global[d];
        bytes = fread(xg, sizeof(double), size, in);
        if ((int)bytes != size) {
          fprintf(stderr,"Error in ReadArray(): Unable to read grid. Expected %d, Read %d.\n",
                  size, (int)bytes);
        }

        /* read solution */
        size = 1;
        for (d = 0; d < a_ndims; d++) size *= a_dim_global[d]; size *= a_nvars;
        bytes = fread(ug, sizeof(double), size, in);
        if ((int)bytes != size) {
          fprintf(stderr,"Error in ReadArray(): Unable to read solution. Expected %d, Read %d.\n",
                  size, (int)bytes);
        }

        fclose(in);
      }

    }
  }

  /* Broadcast a_read_flag to all processes */
  IERR MPIBroadcast_integer(a_read_flag,1,0,&mpi->m_world); CHECKERR(ierr);

  if (*a_read_flag) {

    /* partition global array to all processes */
    IERR MPIPartitionArraynD(a_ndims,mpi,(mpi->m_rank?NULL:ug),a_u,a_dim_global,
                             a_dim_local,a_ghosts,a_nvars); CHECKERR(ierr);

    if (a_x) {
      /* partition a_x vector across the processes */
      int offset_global = 0, offset_local = 0;
      for (d=0; d<a_ndims; d++) {
        IERR MPIPartitionArray1D(mpi,(mpi->m_rank?NULL:&xg[offset_global]),
                                 &a_x[offset_local+a_ghosts],
                                 mpi->m_is[d],mpi->m_ie[d],a_dim_local[d],0); CHECKERR(ierr);
        offset_global += a_dim_global[d];
        offset_local  += a_dim_local [d] + 2*a_ghosts;
      }
    }

    /* free global arrays */
    if (!mpi->m_rank) {
      free(ug);
      free(xg);
    }
  }

  return(0);
}

#ifndef serial

/*! Read in a vector field in a parallel fashion: The number of MPI ranks participating in file I/O
    is specified as an input (#MPIVariables::N_IORanks). All the MPI ranks are divided into that many
    I/O groups, with one rank in each group as the "leader" that does the file reading and writing.
    For reading in the solution, the leader of an I/O group reads its own file and distributes the
    solution to the processors in its group. The number of I/O group is typically specified as the
    number of I/O nodes available on the HPC platform, given the number of compute nodes the code is
    running on. This is a good balance between all the processors serially reading from the same file,
    and having as many files (with the local solution) as the number of processors. This approach has
    been observed to be very scalable (up to ~ 100,000 - 1,000,000 processors).
    \n
    + Supports only binary format.

   There should be as many files as the number of IO ranks (#MPIVariables::N_IORanks).
   The files should be named as: <fname_root>_par.inp.<nnnn>, where <nnnn> is the string
   of formast "%04d" corresponding to integer n, 0 <= n < #MPIVariables::N_IORanks.\n
   Each file should contain the following data:
   \n
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
     and p = i0 + dim_local[0]*( i1 + dim_local[1]*( i2 + dim_local[2]*( ... + dim_global[ndims-2]*i{ndims-1} )))
     (see #_ArrayIndexnD_)\n
     with i0, i1, i2, etc representing grid indices along each spatial dimension, i.e.,\n
     0 <= i0 < dim_local[0]-1\n
     0 <= i1 < dim_local[1]-1\n
     ...\n
     0 <= i{ndims-1} < dim_local[ndims=1]-1\n
   }\n
   for each rank in the IO group corresponding to the file being read.
   + The above block represents the local grid and vector field for each rank
   + Each file should contain as many such blocks of data as there are members in the
     corresponding IO group.
   + The ranks that belong to a particular IO group are given as p, where
     #MPIVariables::GroupStartRank <= p < #MPIVariables::GroupEndRank
   + The code Extras/ParallelInput.c can generate the files <fname_root>_par.inp.<nnnn>
     from the file <fname_root>.inp that is read by ReadArraySerial() if input_mode in
     the input file "solver.inp" is set to "parallel n" where n is the number of IO ranks.
*/
int ReadArrayParallel(
                      int     a_ndims,        /*!< Number of spatial dimensions */
                      int     a_nvars,        /*!< Number of variables per grid point */
                      int     *a_dim_global,  /*!< Integer array of size ndims with global grid size in each dimension */
                      int     *a_dim_local,   /*!< Integer array of size ndims with local  grid size in each dimension */
                      int     a_ghosts,       /*!< Number of ghost points */
                      void    *a_s,           /*!< Solver object of type #HyPar */
                      void    *a_m,           /*!< MPI object of type #MPIVariables */
                      double  *a_x,           /*!< Grid associated with the array (can be NULL) */
                      double  *a_u,           /*!< Array to hold the vector field being read */
                      char    *a_fname_root,  /*!< Filename root */
                      int     *a_read_flag    /*!< Flag to indicate if file was read */
                     )
{
  HyPar         *solver = (HyPar*)        a_s;
  MPIVariables  *mpi    = (MPIVariables*) a_m;
  int           proc,d;
  _DECLARE_IERR_;

  *a_read_flag = 1;
  char filename_root[_MAX_STRING_SIZE_];
  strcpy(filename_root,a_fname_root);
  strcat(filename_root,"_par.inp");

  /* check for existence of the file */
  if (mpi->m_IOParticipant) {
    FILE *in;
    char filename[_MAX_STRING_SIZE_];
    MPIGetFilename(filename_root,&mpi->m_IOWorld,filename);
    in = fopen(filename,"rb");
    if (!in)  *a_read_flag = 0;
    else {
      *a_read_flag = 1;
      fclose(in);
    }
  }
  IERR MPIMin_integer(a_read_flag,a_read_flag,1,&mpi->m_world);

  if (*a_read_flag) {

    if (!mpi->m_rank) printf("Reading from binary file %a_s.xxx (parallel mode).\n",filename_root);

    /* calculate size of the local grid on this rank */
    int sizex = 0;     for (d=0; d<a_ndims; d++) sizex += a_dim_local[d];
    int sizeu = a_nvars; for (d=0; d<a_ndims; d++) sizeu *= a_dim_local[d];

    /* allocate buffer arrays to read in grid and solution */
    double *buffer = (double*) calloc (sizex+sizeu, sizeof(double));

    if (mpi->m_IOParticipant) {

      /* if this rank is responsible for file I/O */
      double *read_buffer = NULL;
      int     read_size_x, read_size_u, read_total_size;
      int     is[a_ndims], ie[a_ndims];

      /* open the file */
      FILE *in;
      int  bytes;
      char filename[_MAX_STRING_SIZE_];
      MPIGetFilename(filename_root,&mpi->m_IOWorld,filename);

      in = fopen(filename,"rb");
      if (!in) {
        fprintf(stderr,"Error in ReadArrayParallel(): File %a_s could not be opened.\n",filename);
        return(1);
      }

      /* Read own data */
      bytes = fread(buffer,sizeof(double),(sizex+sizeu),in);
      if (bytes != (sizex+sizeu)) {
        fprintf(stderr,"Error in ReadArrayParallel(): File %a_s contains insufficient data.\n",filename);
        return(1);
      }

      /* read and send the data for the other processors in this IO rank'a_s group */
      for (proc=mpi->m_GroupStartRank+1; proc<mpi->m_GroupEndRank; proc++) {
        /* get the local domain limits for process proc */
        IERR MPILocalDomainLimits(a_ndims,proc,mpi,a_dim_global,is,ie);
        /* calculate the size of its local data and allocate read buffer */
        read_size_x = 0;      for (d=0; d<a_ndims; d++) read_size_x += (ie[d]-is[d]);
        read_size_u = a_nvars;  for (d=0; d<a_ndims; d++) read_size_u *= (ie[d]-is[d]);
        read_total_size = read_size_x + read_size_u;
        read_buffer = (double*) calloc (read_total_size, sizeof(double));
        /* read the data */
        bytes = fread(read_buffer,sizeof(double),read_total_size,in);
        if (bytes != read_total_size) {
          fprintf(stderr,"Error in ReadArrayParallel(): File %a_s contains insufficient data.\n",filename);
          return(1);
        }
        /* send the data */
        MPI_Request req = MPI_REQUEST_NULL;
        MPI_Isend(read_buffer,read_total_size,MPI_DOUBLE,proc,1100,mpi->m_world,&req);
        MPI_Wait(&req,MPI_STATUS_IGNORE);
        free(read_buffer);
      }

      /* close the file */
      fclose(in);

    } else {

      /* all other processes, just receive the data from
       * the rank responsible for file I/O */
      MPI_Request req = MPI_REQUEST_NULL;
      MPI_Irecv(buffer,(sizex+sizeu),MPI_DOUBLE,mpi->m_IORank,1100,mpi->m_world,&req);
      MPI_Wait(&req,MPI_STATUS_IGNORE);

    }

    /* copy the grid */
    if (a_x) {
      int offset1 = 0, offset2 = 0;
      for (d = 0; d < a_ndims; d++) {
        _ArrayCopy1D_((buffer+offset2),(a_x+offset1+a_ghosts),a_dim_local[d]);
        offset1 += (a_dim_local[d]+2*a_ghosts);
        offset2 +=  a_dim_local[d];
      }
    }

    /* copy the solution */
    int index[a_ndims];
    IERR ArrayCopynD(a_ndims,(buffer+sizex),a_u,a_dim_local,0,a_ghosts,index,a_nvars);
    CHECKERR(ierr);

    /* free buffers */
    free(buffer);
  }

  return(0);

}

/*! Read in an array in a parallel fashion using MPI-IO: Similar to ReadArrayParallel(),
    except that the I/O leaders read from the same file using the MPI I/O routines, by
    calculating their respective offsets and reading the correct chunk of data from that
    offset. The MPI-IO functions (part of MPICH) are constantly being developed to be
    scalable on the latest and greatest HPC platforms.
    \n
    + Supports only binary format.

   There should be as one file named as <fname_root>_mpi.inp. It should contain the
   following data:
   \n
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
     and p = i0 + dim_local[0]*( i1 + dim_local[1]*( i2 + dim_local[2]*( ... + dim_global[ndims-2]*i{ndims-1} )))
     (see #_ArrayIndexnD_)\n
     with i0, i1, i2, etc representing grid indices along each spatial dimension, i.e.,\n
     0 <= i0 < dim_local[0]-1\n
     0 <= i1 < dim_local[1]-1\n
     ...\n
     0 <= i{ndims-1} < dim_local[ndims=1]-1\n
   }\n
   for each rank, in the order of rank number (0 to nproc-1).
   + The above block represents the local grid and vector field for each rank
   + The file should contain as many such blocks of data as there are MPI ranks.
   + Each IO rank computes its offset to figure out where the data it's supposed to
     read exists in the file, and then reads it.
   + The code Extras/MPIInput.c can generate the file <fname_root>_mpi.inp
     from the file <fname_root>.inp that is read by ReadArraySerial() if input_mode in
     the input file "solver.inp" is set to "mpi-io n" where n is the number of IO ranks.
*/
int ReadArrayMPI_IO(
                      int     a_ndims,        /*!< Number of spatial dimensions */
                      int     a_nvars,        /*!< Number of variables per grid point */
                      int     *a_dim_global,  /*!< Integer array of size ndims with global grid size in each dimension */
                      int     *a_dim_local,   /*!< Integer array of size ndims with local  grid size in each dimension */
                      int     a_ghosts,       /*!< Number of ghost points */
                      void    *a_s,           /*!< Solver object of type #HyPar */
                      void    *a_m,           /*!< MPI object of type #MPIVariables */
                      double  *a_x,           /*!< Grid associated with the array (can be NULL) */
                      double  *a_u,           /*!< Array to hold the vector field being read */
                      char    *a_fname_root,  /*!< Filename root */
                      int     *a_read_flag    /*!< Flag to indicate if file was read */
                   )
{
  HyPar         *solver = (HyPar*)        a_s;
  MPIVariables  *mpi    = (MPIVariables*) a_m;
  int           proc,d;
  _DECLARE_IERR_;

  *a_read_flag = 0;
  char filename[_MAX_STRING_SIZE_];
  strcpy(filename,a_fname_root);
  strcat(filename,"_mpi.inp");

  /* check for existence of file */
  if (!mpi->m_rank) {
    FILE *in;
    in = fopen(filename,"rb");
    if (!in)  *a_read_flag = 0;
    else {
      *a_read_flag = 1;
      fclose(in);
    }
  }
  IERR MPIBroadcast_integer(a_read_flag,1,0,&mpi->m_world);

  if (*a_read_flag) {

    if (!mpi->m_rank) printf("Reading from binary file %a_s (MPI-IO mode).\n",filename);

    /* calculate size of the local grid on this rank */
    int sizex = 0;     for (d=0; d<a_ndims; d++) sizex += a_dim_local[d];
    int sizeu = a_nvars; for (d=0; d<a_ndims; d++) sizeu *= a_dim_local[d];

    /* allocate buffer arrays to read in grid and solution */
    double *buffer = (double*) calloc (sizex+sizeu, sizeof(double));

    if (mpi->m_IOParticipant) {

      /* if this rank is responsible for file I/O */
      double *read_buffer = NULL;
      int     read_size_x, read_size_u, read_total_size;
      int     is[a_ndims], ie[a_ndims], size;

      /* calculate offset */
      long long offset = 0;
      for (proc=0; proc < mpi->m_rank; proc++) {
        /* get the local domain limits for process proc */
        IERR MPILocalDomainLimits(a_ndims,proc,mpi,a_dim_global,is,ie);
        /* calculate the size of its local grid */
        size = 0; for (d=0; d<a_ndims; d++) size += (ie[d]-is[d]);
        offset += size;
        /* calculate the size of the local solution */
        size = a_nvars; for (d=0; d<a_ndims; d++) size *= (ie[d]-is[d]);
        offset += size;
      }

      /* open the file */
      MPI_Status  status;
      MPI_File    in;
      int         error;
      error = MPI_File_open(mpi->m_IOWorld,filename,MPI_MODE_RDONLY,MPI_INFO_NULL,&in);
      if (error != MPI_SUCCESS) {
        fprintf(stderr,"Error in ReadArrayMPI_IO(): Unable to open %a_s.\n",filename);
        return(1);
      }

      /* set offset */
      MPI_Offset FileOffset = (MPI_Offset) (offset * sizeof(double));
      MPI_File_seek(in,FileOffset,MPI_SEEK_SET);

      /* Read own data */
      MPI_File_read(in,buffer,(sizex+sizeu)*sizeof(double),MPI_BYTE,&status);

      /* read and send the data for the other processors in this IO rank'a_s group */
      for (proc=mpi->m_GroupStartRank+1; proc<mpi->m_GroupEndRank; proc++) {
        /* get the local domain limits for process proc */
        IERR MPILocalDomainLimits(a_ndims,proc,mpi,a_dim_global,is,ie);
        /* calculate the size of its local data and allocate read buffer */
        read_size_x = 0;      for (d=0; d<a_ndims; d++) read_size_x += (ie[d]-is[d]);
        read_size_u = a_nvars;  for (d=0; d<a_ndims; d++) read_size_u *= (ie[d]-is[d]);
        read_total_size = read_size_x + read_size_u;
        read_buffer = (double*) calloc (read_total_size, sizeof(double));
        /* read the data */
        MPI_File_read(in,read_buffer,read_total_size*sizeof(double),MPI_BYTE,&status);
        /* send the data */
        MPI_Request req = MPI_REQUEST_NULL;
        MPI_Isend(read_buffer,read_total_size,MPI_DOUBLE,proc,1100,mpi->m_world,&req);
        MPI_Wait(&req,MPI_STATUS_IGNORE);
        free(read_buffer);
      }

      /* close the file */
      MPI_File_close(&in);

    } else {

      /* all other processes, just receive the data from
       * the rank responsible for file I/O */
      MPI_Request req = MPI_REQUEST_NULL;
      MPI_Irecv(buffer,(sizex+sizeu),MPI_DOUBLE,mpi->m_IORank,1100,mpi->m_world,&req);
      MPI_Wait(&req,MPI_STATUS_IGNORE);

    }

    /* copy the grid */
    if (a_x) {
      int offset1 = 0, offset2 = 0;
      for (d = 0; d < a_ndims; d++) {
        _ArrayCopy1D_((buffer+offset2),(a_x+offset1+a_ghosts),a_dim_local[d]);
        offset1 += (a_dim_local[d]+2*a_ghosts);
        offset2 +=  a_dim_local[d];
      }
    }

    /* copy the solution */
    int index[a_ndims];
    IERR ArrayCopynD(a_ndims,(buffer+sizex),a_u,a_dim_local,0,a_ghosts,index,a_nvars);
    CHECKERR(ierr);

    /* free buffers */
    free(buffer);
  }

  return(0);

}

#endif

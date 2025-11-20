/*! @file ReadArraywInterp.c
    @author Debojyoti Ghosh
    @brief Read in a vector field from file
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mathfunctions.h>
#include <mpivars.h>
#include <hypar.h>

static int ReadArraywInterpSerial(int,int,int*,int*,int*,int,void*,void*,double*,double*,char*,int*);

/*! Read in a vector field from file: This version allows reading data with different
    dimensions than the array being read in. The data is read in and stored in a new global
    array with the appropriate size, and the array to be read is filled by interpolation.
    Currently, the dimensions of the array to be read and the those of the actual data
    can only differ by factors that are integer powers of 2.

    This is a wrapper function that calls
    the appropriate function depending on input mode (#HyPar::m_input_mode).\n\n
    The mode and type of input are specified through #HyPar::m_input_mode and
    #HyPar::m_ip_file_type. A vector field is read from file and stored in an array.
*/
int ReadArraywInterp( int     a_ndims,          /*!< Number of spatial dimensions */
                      int     a_nvars,          /*!< Number of variables per grid point */
                      int     *a_dim_global,    /*!< Integer array of size ndims with global size in each dimension */
                      int     *a_dim_local,     /*!< Integer array of size ndims with local size in each dimension */
                      int     *a_dim_global_src,/*!< Integer array of size ndims with global size of the data in each dimension */
                      int     a_ghosts,         /*!< Number of ghost points */
                      void    *a_s,             /*!< Solver object of type #HyPar */
                      void    *a_m,             /*!< MPI object of type #MPIVariables */
                      double  *a_x,             /*!< Grid associated with the array (can be NULL) */
                      double  *a_u,             /*!< Array to hold the vector field */
                      char    *a_fname_root,    /*!< Filename root */
                      int     *a_read_flag      /*!< Flag to indicate if the file was read */
                     )
{
  HyPar         *solver = (HyPar*) a_s;
  MPIVariables  *mpi    = (MPIVariables*) a_m;
  _DECLARE_IERR_;

  int retval = ReadArraywInterpSerial(  a_ndims,
                                        a_nvars,
                                        a_dim_global,
                                        a_dim_local,
                                        a_dim_global_src,
                                        a_ghosts,
                                        a_s,
                                        a_m,
                                        a_x,
                                        a_u,
                                        a_fname_root,
                                        a_read_flag );
  if (retval) return retval;

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

  return 0;
}

/*! Read an array in a serial fashion: This version allows reading data with different
    dimensions than the array being read in. The data is read in and stored in a new global
    array with the appropriate size, and the array to be read is filled by interpolation.
    Currently, the dimensions of the array to be read and the those of the actual data
    can only differ by factors that are integer powers of 2.

    For multi-processor simulation, only rank 0
    reads in the entire solution from the file, interpolates it to desired resolution,
    and then distributes the relevant portions
    to each of the processors. This involves memory allocation for the global domain on rank
    0 as well as interpolation between two global arrays. Thus, do not use for large domains. T
    his approach is not very scalable either, if running
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

*/
int ReadArraywInterpSerial( int     a_ndims,          /*!< Number of spatial dimensions */
                            int     a_nvars,          /*!< Number of variables per grid point */
                            int     *a_dim_global,    /*!< Integer array of size ndims with global size in each dimension */
                            int     *a_dim_local,     /*!< Integer array of size ndims with local size in each dimension */
                            int     *a_dim_global_src,/*!< Integer array of size ndims with global size of the data in each dimension */
                            int     a_ghosts,         /*!< Number of ghost points */
                            void    *a_s,             /*!< Solver object of type #HyPar */
                            void    *a_m,             /*!< MPI object of type #MPIVariables */
                            double  *a_x,             /*!< Grid associated with the array (can be NULL) */
                            double  *a_u,             /*!< Array to hold the vector field being read */
                            char    *a_fname_root,    /*!< Filename root */
                            int     *a_read_flag      /*!< Flag to indicate if the file was read */
                          )
{
  HyPar         *solver = (HyPar*)        a_s;
  MPIVariables  *mpi    = (MPIVariables*) a_m;
  int           i, d, ferr, index[a_ndims];
  double        *ug_src = NULL, *xg_src = NULL;
  double        *ug = NULL, *xg = NULL;
  _DECLARE_IERR_;

  *a_read_flag = 0;
  /* Only root process reads from the file */
  if (!mpi->m_rank) {

    /* read data from file - this data is of dimensions given by dim_global_sec */
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
        size   = 1; for (d=0; d<a_ndims; d++) size *= a_dim_global_src[d]; size *= a_nvars;
        ug_src = (double*) calloc(size,sizeof(double));
        size   = 0; for (d=0; d<a_ndims; d++) size += a_dim_global_src[d];
        xg_src = (double*) calloc(size,sizeof(double));

        /* read grid */
        offset = 0;
        for (d = 0; d < a_ndims; d++) {
          for (i = 0; i < a_dim_global_src[d]; i++) {
            ferr = fscanf(in,"%lf",&xg_src[i+offset]);
            if (ferr != 1) {
              printf("Error in ReadArraySerial(): unable to read data. ferr=%d\n", ferr);
              exit(1);
            }
          }
          offset += a_dim_global_src[d];
        }

        /* read solution */
        for (i = 0; i < a_nvars; i++) {
          int done = 0; _ArraySetValue_(index,a_ndims,0);
          while (!done) {
            int p; _ArrayIndex1D_(a_ndims,a_dim_global_src,index,0,p);
            ferr = fscanf(in,"%lf",&ug_src[p*a_nvars+i]);
            if (ferr != 1) {
              printf("Error in ReadArraySerial(): unable to read data. ferr=%d\n", ferr);
              exit(1);
            }
            _ArrayIncrementIndex_(a_ndims,a_dim_global_src,index,done);
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
        size   = 1; for (d=0; d<a_ndims; d++) size *= a_dim_global_src[d]; size *= a_nvars;
        ug_src = (double*) calloc(size,sizeof(double));
        size   = 0; for (d=0; d<a_ndims; d++) size += a_dim_global_src[d];
        xg_src = (double*) calloc(size,sizeof(double));

        /* read grid */
        size = 0;
        for (d = 0; d < a_ndims; d++) size += a_dim_global_src[d];
        bytes = fread(xg_src, sizeof(double), size, in);
        if ((int)bytes != size) {
          fprintf(stderr,"Error in ReadArray(): Unable to read grid. Expected %d, Read %d.\n",
                  size, (int)bytes);
        }

        /* read solution */
        size = 1;
        for (d = 0; d < a_ndims; d++) size *= a_dim_global_src[d]; size *= a_nvars;
        bytes = fread(ug_src, sizeof(double), size, in);
        if ((int)bytes != size) {
          fprintf(stderr,"Error in ReadArray(): Unable to read solution. Expected %d, Read %d.\n",
                  size, (int)bytes);
        }

        fclose(in);
      }

    }

    double *ug_src_wg;
    long size = a_nvars;
    for (d=0; d<a_ndims; d++) size *= (a_dim_global_src[d]+2*a_ghosts);
    ug_src_wg = (double*) calloc (size, sizeof(double));
    ArrayCopynD(  a_ndims,
                  ug_src,
                  ug_src_wg,
                  a_dim_global_src,
                  0,
                  a_ghosts,
                  index,
                  a_nvars );
    FillGhostCells( a_dim_global_src,
                    a_ghosts,
                    ug_src_wg,
                    a_nvars,
                    a_ndims,
                    solver->m_is_periodic);
    free(ug_src);

    /* interpolate from the data read in to a global array
     * with specified dimensions */
    int ierr = InterpolateGlobalnDVar(  a_dim_global,
                                        &ug,
                                        a_dim_global_src,
                                        ug_src_wg,
                                        a_nvars,
                                        a_ghosts,
                                        a_ndims,
                                        solver->m_is_periodic );
    if (ierr) {
      fprintf(stderr, "Error in ReadArraywInterpSerial()\n");
      fprintf(stderr, "  InterpolateGlobalnDVar() returned with error!\n");
      return ierr;
    }

    if (a_x) {
      fprintf(stderr,"Error in ReadArraywInterpSerial()\n");
      fprintf(stderr,"  Not yet implemented for a_x\n");
      //InterpolateGlobal1DVar(xg_src, a_dim_global_src, xg, a_dim_global);
    }
    if (xg_src) free(xg_src);

    if (ug == NULL) {
      fprintf(stderr,"Error in ReadArraywInterp() - ug is NULL!\n");
      return 1;
    }
    if ((xg == NULL) && (a_x != NULL)) {
      fprintf(stderr,"Error in ReadArraywInterp() - xg is NULL!\n");
      return 1;
    }
  }

  /* Broadcast a_read_flag to all processes */
  IERR MPIBroadcast_integer(a_read_flag,1,0,&mpi->m_world); CHECKERR(ierr);

  if (*a_read_flag) {

    /* partition global array to all processes */
    IERR MPIPartitionArraynDwGhosts(  a_ndims,
                                      mpi,
                                      (mpi->m_rank?NULL:ug),
                                      a_u,a_dim_global,
                                      a_dim_local,
                                      a_ghosts,
                                      a_nvars ); CHECKERR(ierr);

//    if (a_x) {
//      /* partition a_x vector across the processes */
//      int offset_global = 0, offset_local = 0;
//      for (d=0; d<a_ndims; d++) {
//        IERR MPIPartitionArray1D(mpi,(mpi->m_rank?NULL:&xg[offset_global]),
//                                 &a_x[offset_local+a_ghosts],
//                                 mpi->m_is[d],mpi->m_ie[d],a_dim_local[d],0); CHECKERR(ierr);
//        offset_global += a_dim_global[d];
//        offset_local  += a_dim_local [d] + 2*a_ghosts;
//      }
//    }

    /* free global arrays */
    if (!mpi->m_rank) {
      free(ug);
      free(xg);
    }
  }

  return(0);
}


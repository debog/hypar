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
    the appropriate function depending on input mode (#HyPar::input_mode).\n\n
    The mode and type of input are specified through #HyPar::input_mode and
    #HyPar::ip_file_type. A vector field is read from file and stored in an array.
*/
int ReadArraywInterp( int     ndims,          /*!< Number of spatial dimensions */
                      int     nvars,          /*!< Number of variables per grid point */
                      int     *dim_global,    /*!< Integer array of size ndims with global size in each dimension */
                      int     *dim_local,     /*!< Integer array of size ndims with local size in each dimension */
                      int     *dim_global_src,/*!< Integer array of size ndims with global size of the data in each dimension */
                      int     ghosts,         /*!< Number of ghost points */
                      void    *s,             /*!< Solver object of type #HyPar */
                      void    *m,             /*!< MPI object of type #MPIVariables */
                      double  *x,             /*!< Grid associated with the array (can be NULL) */
                      double  *u,             /*!< Array to hold the vector field */
                      char    *fname_root,    /*!< Filename root */
                      int     *read_flag      /*!< Flag to indicate if the file was read */
                     )
{
  HyPar         *solver = (HyPar*) s;
  MPIVariables  *mpi    = (MPIVariables*) m;
  _DECLARE_IERR_;

  int retval = ReadArraywInterpSerial(  ndims,
                                        nvars,
                                        dim_global,
                                        dim_local,
                                        dim_global_src,
                                        ghosts,
                                        s,
                                        m,
                                        x,
                                        u,
                                        fname_root,
                                        read_flag );
  if (retval) return retval;

  if (x) {
    int offset, d;
    /* exchange MPI-boundary values of x between processors */
    offset = 0;
    for (d = 0; d < ndims; d++) {
      IERR MPIExchangeBoundaries1D(mpi,&x[offset],dim_local[d],
                                   ghosts,d,ndims); CHECKERR(ierr);
      offset  += (dim_local [d] + 2*ghosts);
    }
    /* fill in ghost values of x at physical boundaries by extrapolation */
    offset = 0;
    for (d = 0; d < ndims; d++) {
      double *X     = &x[offset];
      int    *dim   = dim_local, i;
      if (mpi->ip[d] == 0) {
        /* fill left boundary along this dimension */
        for (i = 0; i < ghosts; i++) {
          int delta = ghosts - i;
          X[i] = X[ghosts] + ((double) delta) * (X[ghosts]-X[ghosts+1]);
        }
      }
      if (mpi->ip[d] == mpi->iproc[d]-1) {
        /* fill right boundary along this dimension */
        for (i = dim[d]+ghosts; i < dim[d]+2*ghosts; i++) {
          int delta = i - (dim[d]+ghosts-1);
          X[i] =  X[dim[d]+ghosts-1]
                  + ((double) delta) * (X[dim[d]+ghosts-1]-X[dim[d]+ghosts-2]);
        }
      }
      offset  += (dim[d] + 2*ghosts);
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
int ReadArraywInterpSerial( int     ndims,          /*!< Number of spatial dimensions */
                            int     nvars,          /*!< Number of variables per grid point */
                            int     *dim_global,    /*!< Integer array of size ndims with global size in each dimension */
                            int     *dim_local,     /*!< Integer array of size ndims with local size in each dimension */
                            int     *dim_global_src,/*!< Integer array of size ndims with global size of the data in each dimension */
                            int     ghosts,         /*!< Number of ghost points */
                            void    *s,             /*!< Solver object of type #HyPar */
                            void    *m,             /*!< MPI object of type #MPIVariables */
                            double  *x,             /*!< Grid associated with the array (can be NULL) */
                            double  *u,             /*!< Array to hold the vector field being read */
                            char    *fname_root,    /*!< Filename root */
                            int     *read_flag      /*!< Flag to indicate if the file was read */
                          )
{
  HyPar         *solver = (HyPar*)        s;
  MPIVariables  *mpi    = (MPIVariables*) m;
  int           i, d, ferr, index[ndims];
  double        *ug_src = NULL, *xg_src = NULL;
  double        *ug = NULL, *xg = NULL;
  _DECLARE_IERR_;

  *read_flag = 0;
  /* Only root process reads from the file */
  if (!mpi->rank) {

    /* read data from file - this data is of dimensions given by dim_global_sec */
    if (!strcmp(solver->ip_file_type,"ascii")) {
      char filename[_MAX_STRING_SIZE_];
      strcpy(filename,fname_root);
      strcat(filename,".inp");
      FILE *in; in = fopen(filename,"r");
      if (!in) *read_flag = 0;
      else {
        *read_flag = 1;
        /* Reading from file */
        printf("Reading array from ASCII file %s (Serial mode).\n",filename);
        int size,offset;
        /* allocate global solution array */
        size   = 1; for (d=0; d<ndims; d++) size *= dim_global_src[d]; size *= nvars;
        ug_src = (double*) calloc(size,sizeof(double));
        size   = 0; for (d=0; d<ndims; d++) size += dim_global_src[d];
        xg_src = (double*) calloc(size,sizeof(double));

        /* read grid */
        offset = 0;
        for (d = 0; d < ndims; d++) {
          for (i = 0; i < dim_global_src[d]; i++) {
            ferr = fscanf(in,"%lf",&xg_src[i+offset]);
            if (ferr != 1) {
              printf("Error in ReadArraySerial(): unable to read data. ferr=%d\n", ferr);
              exit(1);
            }
          }
          offset += dim_global_src[d];
        }

        /* read solution */
        for (i = 0; i < nvars; i++) {
          int done = 0; _ArraySetValue_(index,ndims,0);
          while (!done) {
            int p; _ArrayIndex1D_(ndims,dim_global_src,index,0,p);
            ferr = fscanf(in,"%lf",&ug_src[p*nvars+i]);
            if (ferr != 1) {
              printf("Error in ReadArraySerial(): unable to read data. ferr=%d\n", ferr);
              exit(1);
            }
            _ArrayIncrementIndex_(ndims,dim_global_src,index,done);
          }
        }

        fclose(in);
      }
    } else if ((!strcmp(solver->ip_file_type,"bin")) || (!strcmp(solver->ip_file_type,"binary"))) {

      char filename[_MAX_STRING_SIZE_];
      strcpy(filename,fname_root);
      strcat(filename,".inp");
      FILE *in; in = fopen(filename,"rb");
      if (!in) *read_flag = 0;
      else {
        *read_flag = 1;
        printf("Reading array from binary file %s (Serial mode).\n",filename);
        size_t bytes;
        int size;
        /* allocate global solution array */
        size   = 1; for (d=0; d<ndims; d++) size *= dim_global_src[d]; size *= nvars;
        ug_src = (double*) calloc(size,sizeof(double));
        size   = 0; for (d=0; d<ndims; d++) size += dim_global_src[d];
        xg_src = (double*) calloc(size,sizeof(double));

        /* read grid */
        size = 0;
        for (d = 0; d < ndims; d++) size += dim_global_src[d];
        bytes = fread(xg_src, sizeof(double), size, in);
        if ((int)bytes != size) {
          fprintf(stderr,"Error in ReadArray(): Unable to read grid. Expected %d, Read %d.\n",
                  size, (int)bytes);
        }

        /* read solution */
        size = 1;
        for (d = 0; d < ndims; d++) size *= dim_global_src[d]; size *= nvars;
        bytes = fread(ug_src, sizeof(double), size, in);
        if ((int)bytes != size) {
          fprintf(stderr,"Error in ReadArray(): Unable to read solution. Expected %d, Read %d.\n",
                  size, (int)bytes);
        }

        fclose(in);
      }

    }

    double *ug_src_wg;
    long size = nvars;
    for (d=0; d<ndims; d++) size *= (dim_global_src[d]+2*ghosts);
    ug_src_wg = (double*) calloc (size, sizeof(double));
    ArrayCopynD(  ndims,
                  ug_src,
                  ug_src_wg,
                  dim_global_src,
                  0,
                  ghosts,
                  index,
                  nvars );
    fillGhostCells( dim_global_src,
                    ghosts,
                    ug_src_wg,
                    nvars,
                    ndims,
                    solver->isPeriodic);
    free(ug_src);

    /* interpolate from the data read in to a global array
     * with specified dimensions */
    int ierr = InterpolateGlobalnDVar(  dim_global,
                                        &ug,
                                        dim_global_src,
                                        ug_src_wg,
                                        nvars,
                                        ghosts,
                                        ndims,
                                        solver->isPeriodic );
    if (ierr) {
      fprintf(stderr, "Error in ReadArraywInterpSerial()\n");
      fprintf(stderr, "  InterpolateGlobalnDVar() returned with error!\n");
      return ierr;
    }

    if (x) {
      fprintf(stderr,"Error in ReadArraywInterpSerial()\n");
      fprintf(stderr,"  Not yet implemented for x\n");
      //InterpolateGlobal1DVar(xg_src, dim_global_src, xg, dim_global);
    }
    if (xg_src) free(xg_src);

    if (ug == NULL) {
      fprintf(stderr,"Error in ReadArraywInterp() - ug is NULL!\n");
      return 1;
    }
    if ((xg == NULL) && (x != NULL)) {
      fprintf(stderr,"Error in ReadArraywInterp() - xg is NULL!\n");
      return 1;
    }
  }

  /* Broadcast read_flag to all processes */
  IERR MPIBroadcast_integer(read_flag,1,0,&mpi->world); CHECKERR(ierr);

  if (*read_flag) {

    /* partition global array to all processes */
    IERR MPIPartitionArraynDwGhosts(  ndims,
                                      mpi,
                                      (mpi->rank?NULL:ug),
                                      u,dim_global,
                                      dim_local,
                                      ghosts,
                                      nvars ); CHECKERR(ierr);

//    if (x) {
//      /* partition x vector across the processes */
//      int offset_global = 0, offset_local = 0;
//      for (d=0; d<ndims; d++) {
//        IERR MPIPartitionArray1D(mpi,(mpi->rank?NULL:&xg[offset_global]),
//                                 &x[offset_local+ghosts],
//                                 mpi->is[d],mpi->ie[d],dim_local[d],0); CHECKERR(ierr);
//        offset_global += dim_global[d];
//        offset_local  += dim_local [d] + 2*ghosts;
//      }
//    }

    /* free global arrays */
    if (!mpi->rank) {
      free(ug);
      free(xg);
    }
  }

  return(0);
}


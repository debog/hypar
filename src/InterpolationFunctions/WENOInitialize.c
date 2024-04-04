/*! @file WENOInitialize.c
    @brief Initializes the WENO-type schemes
    @author Debojyoti Ghosh, Youngdae Kim
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#if defined(HAVE_CUDA)
#include <arrayfunctions_gpu.h>
#else
#include <arrayfunctions.h>
#endif
#include <interpolation.h>
#include <mpivars.h>
#include <hypar.h>

int WENOFifthOrderInitializeWeights   (double* const,double* const,double* const,
                                       const int* const, int,void*,void*);
int WENOFifthOrderCalculateWeights    (double*,double*,double*,int,void*,void*);
int WENOFifthOrderCalculateWeightsChar(double*,double*,double*,int,void*,void*);

#if defined(HAVE_CUDA)
int gpuWENOFifthOrderCalculateWeights   (double*,double*,double*,int,void*,void*);
#endif

/*!
  This function initializes the WENO-type methods.
  + Sets the parameters to default values.
  + Reads in the parameters from optional input file "weno.inp", if available.
  + Allocates memory for and initializes the nonlinear weights used by WENO-type
    schemes.
*/
int WENOInitialize(
                    void *s,      /*!< Solver object of type #HyPar */
                    void *m,      /*!< MPI object of type #MPIVariables */
                    char *scheme, /*!< Name of scheme */
                    char *type    /*!< Type of interpolation */
                  )
{
  HyPar           *solver = (HyPar*) s;
  MPIVariables    *mpi    = (MPIVariables*) m;
  WENOParameters  *weno   = (WENOParameters*) solver->interp;

  static int count = 0;

  int nvars = solver->nvars;
  int ndims = solver->ndims;

  /* default parameters */
  weno->mapped      = 0;
  weno->borges      = 0;
  weno->yc          = 0;
  weno->no_limiting = 0;
  weno->eps         = 1e-6;
  weno->p           = 2.0;

  weno->rc          = 0.3;
  weno->xi          = 0.001;
  weno->tol         = 1e-16;

  if (!mpi->rank) {
    FILE *in;
    int ferr;
    in = fopen("weno.inp","r");
    if (!in) printf("Warning: File weno.inp not found. Using default parameters for WENO5/CRWENO5/HCWENO5 scheme.\n");
    else {
      if (!count) printf("Reading WENO parameters from weno.inp.\n");
      char word[_MAX_STRING_SIZE_];
      ferr = fscanf(in,"%s",word); if (ferr != 1) return(1);
      if (!strcmp(word, "begin")){
        while (strcmp(word, "end")){
          ferr = fscanf(in,"%s",word); if (ferr != 1) return(1);
          if      (!strcmp(word,"mapped"     )) { ferr = fscanf(in,"%d" ,&weno->mapped     ); if (ferr != 1) return(1); }
          else if (!strcmp(word,"borges"     )) { ferr = fscanf(in,"%d" ,&weno->borges     ); if (ferr != 1) return(1); }
          else if (!strcmp(word,"yc"         )) { ferr = fscanf(in,"%d" ,&weno->yc         ); if (ferr != 1) return(1); }
          else if (!strcmp(word,"no_limiting")) { ferr = fscanf(in,"%d" ,&weno->no_limiting); if (ferr != 1) return(1); }
          else if (!strcmp(word,"epsilon"    )) { ferr = fscanf(in,"%lf",&weno->eps        ); if (ferr != 1) return(1); }
          else if (!strcmp(word,"p"          )) { ferr = fscanf(in,"%lf",&weno->p          ); if (ferr != 1) return(1); }
          else if (!strcmp(word,"rc"         )) { ferr = fscanf(in,"%lf",&weno->rc         ); if (ferr != 1) return(1); }
          else if (!strcmp(word,"xi"         )) { ferr = fscanf(in,"%lf",&weno->xi         ); if (ferr != 1) return(1); }
          else if (!strcmp(word,"tol"        )) { ferr = fscanf(in,"%lf",&weno->tol        ); if (ferr != 1) return(1); }
          else if (strcmp(word,"end")) {
            char useless[_MAX_STRING_SIZE_];
            ferr = fscanf(in,"%s",useless); if (ferr != 1) return(ferr);
            printf("Warning: keyword %s in file \"weno.inp\" with value %s not ",word,useless);
            printf("recognized or extraneous. Ignoring.\n");
          }
        }
      } else {
        fprintf(stderr,"Error: Illegal format in file \"weno.inp\".\n");
        return(1);
      }
      fclose(in);
    }
  }

  int     integer_data[4];
  double  real_data[5];
  if (!mpi->rank) {
    integer_data[0] = weno->mapped;
    integer_data[1] = weno->borges;
    integer_data[2] = weno->yc;
    integer_data[3] = weno->no_limiting;
    real_data[0]    = weno->eps;
    real_data[1]    = weno->p;
    real_data[2]    = weno->rc;
    real_data[3]    = weno->xi;
    real_data[4]    = weno->tol;
  }
  MPIBroadcast_integer(integer_data,4,0,&mpi->world);
  MPIBroadcast_double (real_data   ,5,0,&mpi->world);

  weno->mapped      = integer_data[0];
  weno->borges      = integer_data[1];
  weno->yc          = integer_data[2];
  weno->no_limiting = integer_data[3];
  weno->eps         = real_data   [0];
  weno->p           = real_data   [1];
  weno->rc          = real_data   [2];
  weno->xi          = real_data   [3];
  weno->tol         = real_data   [4];

  /* WENO weight calculation is hard-coded for p=2, so return error if p != 2 in
   * user input file, so that there's no confusion */
  if (weno->p != 2.0) {
    if (!mpi->rank && !count) printf("Warning from WENOInitialize(): \"p\" parameter is 2.0. Any other value will be ignored!\n");
  }

  weno->offset = NULL;
  weno->w1 = NULL;
  weno->w2 = NULL;
  weno->w3 = NULL;

  weno->offset = (int*) calloc (ndims,sizeof(int));
  int dir,d;
  for (dir=0; dir<ndims; dir++) {
    weno->offset[dir] = 0;
    for (d=0; d<dir; d++) {
      int size = nvars, i;
      for (i=0; i<ndims; i++)
        size *= ( i==d ? solver->dim_local[i]+1 : solver->dim_local[i] );
      weno->offset[dir] += size;
    }
  }

  int total_size = 0;
  for (d=0; d<ndims; d++) {
    int size = nvars, i;
    for (i=0; i<ndims; i++)
      size *= ( i==d ? solver->dim_local[i]+1 : solver->dim_local[i] );
    total_size += size;
  }
  weno->size = total_size;

  if ((!strcmp(type,_CHARACTERISTIC_)) && (nvars > 1))
    solver->SetInterpLimiterVar = WENOFifthOrderCalculateWeightsChar;
  else {
#if defined(HAVE_CUDA)
    if (solver->use_gpu) {
      solver->SetInterpLimiterVar = gpuWENOFifthOrderCalculateWeights;
    } else {
#endif
      solver->SetInterpLimiterVar = WENOFifthOrderCalculateWeights;
#if defined(HAVE_CUDA)
    }
#endif
  }

  /* initialize WENO weights to their optimal values */
  double* tmp_w1 = (double*) calloc (4*total_size,sizeof(double));
  double* tmp_w2 = (double*) calloc (4*total_size,sizeof(double));
  double* tmp_w3 = (double*) calloc (4*total_size,sizeof(double));
  for (d=0; d<ndims; d++) WENOFifthOrderInitializeWeights(  tmp_w1,
                                                            tmp_w2,
                                                            tmp_w3,
                                                            weno->offset,
                                                            d,
                                                            solver,
                                                            mpi);
  count++;

#if defined(HAVE_CUDA)
  if (solver->use_gpu) {
    //gpuMalloc((void**)&weno->gpu_offset, solver->ndims*sizeof(int));
    gpuMalloc((void**)&weno->w1, 4*total_size*sizeof(double));
    gpuMalloc((void**)&weno->w2, 4*total_size*sizeof(double));
    gpuMalloc((void**)&weno->w3, 4*total_size*sizeof(double));

    //gpuMemcpy(weno->gpu_offset, weno->offset, solver->ndims*sizeof(int), gpuMemcpyHostToDevice);
    gpuMemcpy(weno->w1, tmp_w1, 4*total_size*sizeof(double), gpuMemcpyHostToDevice);
    gpuMemcpy(weno->w2, tmp_w2, 4*total_size*sizeof(double), gpuMemcpyHostToDevice);
    gpuMemcpy(weno->w3, tmp_w3, 4*total_size*sizeof(double), gpuMemcpyHostToDevice);
  } else {
#endif
    weno->w1 = (double*) calloc (4*total_size,sizeof(double));
    weno->w2 = (double*) calloc (4*total_size,sizeof(double));
    weno->w3 = (double*) calloc (4*total_size,sizeof(double));

    _ArrayCopy1D_(tmp_w1, weno->w1, 4*total_size);
    _ArrayCopy1D_(tmp_w2, weno->w2, 4*total_size);
    _ArrayCopy1D_(tmp_w3, weno->w3, 4*total_size);
#if defined(HAVE_CUDA)
  }
#endif

  free(tmp_w1);
  free(tmp_w2);
  free(tmp_w3);

  return 0;
}

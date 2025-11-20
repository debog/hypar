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
  WENOParameters  *weno   = (WENOParameters*) solver->m_interp;

  static int count = 0;

  int nvars = solver->m_nvars;
  int ndims = solver->m_ndims;

  /* default parameters */
  weno->m_mapped      = 0;
  weno->m_borges      = 0;
  weno->m_yc          = 0;
  weno->m_no_limiting = 0;
  weno->m_eps         = 1e-6;
  weno->m_p           = 2.0;

  weno->m_rc          = 0.3;
  weno->m_xi          = 0.001;
  weno->m_tol         = 1e-16;

  if (!mpi->m_rank) {
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
          if      (!strcmp(word,"mapped"     )) { ferr = fscanf(in,"%d" ,&weno->m_mapped     ); if (ferr != 1) return(1); }
          else if (!strcmp(word,"borges"     )) { ferr = fscanf(in,"%d" ,&weno->m_borges     ); if (ferr != 1) return(1); }
          else if (!strcmp(word,"yc"         )) { ferr = fscanf(in,"%d" ,&weno->m_yc         ); if (ferr != 1) return(1); }
          else if (!strcmp(word,"no_limiting")) { ferr = fscanf(in,"%d" ,&weno->m_no_limiting); if (ferr != 1) return(1); }
          else if (!strcmp(word,"epsilon"    )) { ferr = fscanf(in,"%lf",&weno->m_eps        ); if (ferr != 1) return(1); }
          else if (!strcmp(word,"p"          )) { ferr = fscanf(in,"%lf",&weno->m_p          ); if (ferr != 1) return(1); }
          else if (!strcmp(word,"rc"         )) { ferr = fscanf(in,"%lf",&weno->m_rc         ); if (ferr != 1) return(1); }
          else if (!strcmp(word,"xi"         )) { ferr = fscanf(in,"%lf",&weno->m_xi         ); if (ferr != 1) return(1); }
          else if (!strcmp(word,"tol"        )) { ferr = fscanf(in,"%lf",&weno->m_tol        ); if (ferr != 1) return(1); }
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
  if (!mpi->m_rank) {
    integer_data[0] = weno->m_mapped;
    integer_data[1] = weno->m_borges;
    integer_data[2] = weno->m_yc;
    integer_data[3] = weno->m_no_limiting;
    real_data[0]    = weno->m_eps;
    real_data[1]    = weno->m_p;
    real_data[2]    = weno->m_rc;
    real_data[3]    = weno->m_xi;
    real_data[4]    = weno->m_tol;
  }
  MPIBroadcast_integer(integer_data,4,0,&mpi->m_world);
  MPIBroadcast_double (real_data   ,5,0,&mpi->m_world);

  weno->m_mapped      = integer_data[0];
  weno->m_borges      = integer_data[1];
  weno->m_yc          = integer_data[2];
  weno->m_no_limiting = integer_data[3];
  weno->m_eps         = real_data   [0];
  weno->m_p           = real_data   [1];
  weno->m_rc          = real_data   [2];
  weno->m_xi          = real_data   [3];
  weno->m_tol         = real_data   [4];

  /* WENO weight calculation is hard-coded for p=2, so return error if p != 2 in
   * user input file, so that there's no confusion */
  if (weno->m_p != 2.0) {
    if (!mpi->m_rank && !count) printf("Warning from WENOInitialize(): \"p\" parameter is 2.0. Any other value will be ignored!\n");
  }

  weno->m_offset = NULL;
  weno->m_w1 = NULL;
  weno->m_w2 = NULL;
  weno->m_w3 = NULL;

  weno->m_offset = (int*) calloc (ndims,sizeof(int));
  int dir,d;
  for (dir=0; dir<ndims; dir++) {
    weno->m_offset[dir] = 0;
    for (d=0; d<dir; d++) {
      int size = nvars, i;
      for (i=0; i<ndims; i++)
        size *= ( i==d ? solver->m_dim_local[i]+1 : solver->m_dim_local[i] );
      weno->m_offset[dir] += size;
    }
  }

  int total_size = 0;
  for (d=0; d<ndims; d++) {
    int size = nvars, i;
    for (i=0; i<ndims; i++)
      size *= ( i==d ? solver->m_dim_local[i]+1 : solver->m_dim_local[i] );
    total_size += size;
  }
  weno->m_size = total_size;

  if ((!strcmp(type,_CHARACTERISTIC_)) && (nvars > 1))
    solver->SetInterpLimiterVar = WENOFifthOrderCalculateWeightsChar;
  else {
#if defined(HAVE_CUDA)
    if (solver->m_use_gpu) {
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
                                                            weno->m_offset,
                                                            d,
                                                            solver,
                                                            mpi);
  count++;

#if defined(HAVE_CUDA)
  if (solver->m_use_gpu) {
    //gpuMalloc((void**)&weno->gpu_offset, solver->m_ndims*sizeof(int));
    gpuMalloc((void**)&weno->m_w1, 4*total_size*sizeof(double));
    gpuMalloc((void**)&weno->m_w2, 4*total_size*sizeof(double));
    gpuMalloc((void**)&weno->m_w3, 4*total_size*sizeof(double));

    //gpuMemcpy(weno->gpu_offset, weno->m_offset, solver->m_ndims*sizeof(int), gpuMemcpyHostToDevice);
    gpuMemcpy(weno->m_w1, tmp_w1, 4*total_size*sizeof(double), gpuMemcpyHostToDevice);
    gpuMemcpy(weno->m_w2, tmp_w2, 4*total_size*sizeof(double), gpuMemcpyHostToDevice);
    gpuMemcpy(weno->m_w3, tmp_w3, 4*total_size*sizeof(double), gpuMemcpyHostToDevice);
  } else {
#endif
    weno->m_w1 = (double*) calloc (4*total_size,sizeof(double));
    weno->m_w2 = (double*) calloc (4*total_size,sizeof(double));
    weno->m_w3 = (double*) calloc (4*total_size,sizeof(double));

    _ArrayCopy1D_(tmp_w1, weno->m_w1, 4*total_size);
    _ArrayCopy1D_(tmp_w2, weno->m_w2, 4*total_size);
    _ArrayCopy1D_(tmp_w3, weno->m_w3, 4*total_size);
#if defined(HAVE_CUDA)
  }
#endif

  free(tmp_w1);
  free(tmp_w2);
  free(tmp_w3);

  return 0;
}

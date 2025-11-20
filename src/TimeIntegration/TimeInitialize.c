/*! @file TimeInitialize.c
    @brief Initialize time integration
    @author Debojyoti Ghosh, Youngdae Kim
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <basic.h>
#if defined(HAVE_CUDA)
#include <arrayfunctions_gpu.h>
#else
#include <arrayfunctions.h>
#endif
#include <simulation_object.h>
#include <timeintegration.h>

int TimeRHSFunctionExplicit(double*,double*,void*,void*,double);

/*!
  Initialize time integration: This function initializes all that is required
  for time integration.
  + It sets the number of iterations, time step size, simulation time, etc.
  + It allocates solution, right-hand-side, and stage solution arrays needed
    by specific time integration methods.
  + It calls the method-specific initialization functions.
*/
int TimeInitialize( void  *a_s,     /*!< Array of simulation objects of type #SimulationObject */
                    int   a_nsims,  /*!< number of simulation objects */
                    int   a_rank,   /*!< MPI rank of this process */
                    int   a_nproc,  /*!< number of MPI processes */
                    void  *a_ts     /*!< Time integration object of type #TimeIntegration */
                  )
{
  TimeIntegration*  TS  = (TimeIntegration*) a_ts;
  SimulationObject* sim = (SimulationObject*) a_s;
  int ns, d, i;

  if (!sim) return(1);


  TS->m_simulation    = sim;
  TS->m_nsims         = a_nsims;

  TS->m_rank          = a_rank;
  TS->m_nproc         = a_nproc;

  TS->m_n_iter        = sim[0].solver.m_n_iter;
  TS->m_restart_iter  = sim[0].solver.m_restart_iter;
  TS->m_dt            = sim[0].solver.m_dt;

  TS->m_waqt          = (double) TS->m_restart_iter * TS->m_dt;
  TS->m_max_cfl       = 0.0;
  TS->m_norm          = 0.0;
  TS->TimeIntegrate = sim[0].solver.TimeIntegrate;
  TS->m_iter_wctime_total = 0.0;

  TS->m_u_offsets = (long*) calloc (a_nsims, sizeof(long));
  TS->m_u_sizes = (long*) calloc (a_nsims, sizeof(long));

  TS->m_u_size_total = 0;
  for (ns = 0; ns < a_nsims; ns++) {
    TS->m_u_offsets[ns] = TS->m_u_size_total;
    TS->m_u_sizes[ns] = sim[ns].solver.m_npoints_local_wghosts * sim[ns].solver.m_nvars ;
    TS->m_u_size_total += TS->m_u_sizes[ns];
  }

  TS->m_u   = (double*) calloc (TS->m_u_size_total,sizeof(double));
  TS->m_rhs = (double*) calloc (TS->m_u_size_total,sizeof(double));
  _ArraySetValue_(TS->m_u  ,TS->m_u_size_total,0.0);
  _ArraySetValue_(TS->m_rhs,TS->m_u_size_total,0.0);

  /* initialize arrays to NULL, then allocate as necessary */
  TS->m_U             = NULL;
  TS->m_Udot          = NULL;
  TS->m_BoundaryFlux  = NULL;

  TS->m_bf_offsets = (int*) calloc (a_nsims, sizeof(int));
  TS->m_bf_sizes = (int*) calloc (a_nsims, sizeof(int));
  TS->m_bf_size_total = 0;
  for (ns = 0; ns < a_nsims; ns++) {
    TS->m_bf_offsets[ns] = TS->m_bf_size_total;
    TS->m_bf_sizes[ns] =  2 * sim[ns].solver.m_ndims * sim[ns].solver.m_nvars;
    TS->m_bf_size_total += TS->m_bf_sizes[ns];
  }

#if defined(HAVE_CUDA)
  if (sim[0].solver.m_use_gpu) {

    if (!strcmp(sim[0].solver.m_time_scheme,_RK_)) {

      /* explicit Runge-Kutta methods */
      ExplicitRKParameters  *params = (ExplicitRKParameters*)  sim[0].solver.m_msti;
      int nstages = params->nstages;

      gpuMalloc((void**)&TS->m_gpu_U, TS->m_u_size_total*nstages*sizeof(double));
      gpuMalloc((void**)&TS->m_gpu_Udot, TS->m_u_size_total*nstages*sizeof(double));
      gpuMalloc((void**)&TS->m_gpu_BoundaryFlux, TS->m_bf_size_total*nstages*sizeof(double));
      gpuMemset(TS->m_gpu_U, 0, TS->m_u_size_total*nstages*sizeof(double));
      gpuMemset(TS->m_gpu_Udot, 0, TS->m_u_size_total*nstages*sizeof(double));
      gpuMemset(TS->m_gpu_BoundaryFlux, 0, TS->m_bf_size_total*nstages*sizeof(double));

    } else {

      fprintf(stderr,"ERROR in TimeInitialize(): %a_s is not yet implemented on GPUs.\n",
              sim[0].solver.m_time_scheme );
      return 1;

    }

  } else {
#endif

    if (!strcmp(sim[0].solver.m_time_scheme,_RK_)) {

      /* explicit Runge-Kutta methods */
      ExplicitRKParameters  *params = (ExplicitRKParameters*)  sim[0].solver.m_msti;
      int nstages = params->nstages;
      TS->m_U     = (double**) calloc (nstages,sizeof(double*));
      TS->m_Udot  = (double**) calloc (nstages,sizeof(double*));
      for (i = 0; i < nstages; i++) {
        TS->m_U[i]    = (double*) calloc (TS->m_u_size_total,sizeof(double));
        TS->m_Udot[i] = (double*) calloc (TS->m_u_size_total,sizeof(double));
      }

      TS->m_BoundaryFlux = (double**) calloc (nstages,sizeof(double*));
      for (i=0; i<nstages; i++) {
        TS->m_BoundaryFlux[i] = (double*) calloc (TS->m_bf_size_total,sizeof(double));
      }

    } else if (!strcmp(sim[0].solver.m_time_scheme,_FORWARD_EULER_)) {

      int nstages = 1;
      TS->m_BoundaryFlux = (double**) calloc (nstages,sizeof(double*));
      for (i=0; i<nstages; i++) {
        TS->m_BoundaryFlux[i] = (double*) calloc (TS->m_bf_size_total,sizeof(double));
      }

    } else if (!strcmp(sim[0].solver.m_time_scheme,_GLM_GEE_)) {

      /* General Linear Methods with Global Error Estimate */
      GLMGEEParameters *params = (GLMGEEParameters*) sim[0].solver.m_msti;
      int nstages = params->nstages;
      int r       = params->r;
      TS->m_U     = (double**) calloc (2*r-1  ,sizeof(double*));
      TS->m_Udot  = (double**) calloc (nstages,sizeof(double*));
      for (i=0; i<2*r-1; i++)   TS->m_U[i]    = (double*) calloc (TS->m_u_size_total,sizeof(double));
      for (i=0; i<nstages; i++) TS->m_Udot[i] = (double*) calloc (TS->m_u_size_total,sizeof(double));

      TS->m_BoundaryFlux = (double**) calloc (nstages,sizeof(double*));
      for (i=0; i<nstages; i++) {
        TS->m_BoundaryFlux[i] = (double*) calloc (TS->m_bf_size_total,sizeof(double));
      }


      if (!strcmp(params->ee_mode,_GLM_GEE_YYT_)) {
        for (ns = 0; ns < a_nsims; ns++) {
          for (i=0; i<r-1; i++) {
            _ArrayCopy1D_(  (sim[ns].solver.m_u),
                            (TS->m_U[r+i] + TS->m_u_offsets[ns]),
                            (TS->m_u_sizes[ns])  );
          }
        }
      } else {
        for (i=0; i<r-1; i++) {
          _ArraySetValue_(TS->m_U[r+i],TS->m_u_size_total,0.0);
        }
      }
    }
#if defined(HAVE_CUDA)
  }
#endif

  /* set right-hand side function pointer */
  TS->RHSFunction = TimeRHSFunctionExplicit;

  /* open files for writing */
  if (!a_rank) {
    if (sim[0].solver.m_write_residual) TS->m_ResidualFile = (void*) fopen("residual.out","w");
    else                              TS->m_ResidualFile = NULL;
  } else                              TS->m_ResidualFile = NULL;

  for (ns = 0; ns < a_nsims; ns++) {
    sim[ns].solver.m_time_integrator = TS;
  }

  return 0;
}


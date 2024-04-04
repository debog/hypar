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
int TimeInitialize( void  *s,     /*!< Array of simulation objects of type #SimulationObject */
                    int   nsims,  /*!< number of simulation objects */
                    int   rank,   /*!< MPI rank of this process */
                    int   nproc,  /*!< number of MPI processes */
                    void  *ts     /*!< Time integration object of type #TimeIntegration */
                  )
{
  TimeIntegration*  TS  = (TimeIntegration*) ts;
  SimulationObject* sim = (SimulationObject*) s;
  int ns, d, i;

  if (!sim) return(1);


  TS->simulation    = sim;
  TS->nsims         = nsims;

  TS->rank          = rank;
  TS->nproc         = nproc;

  TS->n_iter        = sim[0].solver.n_iter;
  TS->restart_iter  = sim[0].solver.restart_iter;
  TS->dt            = sim[0].solver.dt;

  TS->waqt          = (double) TS->restart_iter * TS->dt;
  TS->max_cfl       = 0.0;
  TS->norm          = 0.0;
  TS->TimeIntegrate = sim[0].solver.TimeIntegrate;
  TS->iter_wctime_total = 0.0;

  TS->u_offsets = (long*) calloc (nsims, sizeof(long));
  TS->u_sizes = (long*) calloc (nsims, sizeof(long));

  TS->u_size_total = 0;
  for (ns = 0; ns < nsims; ns++) {
    TS->u_offsets[ns] = TS->u_size_total;
    TS->u_sizes[ns] = sim[ns].solver.npoints_local_wghosts * sim[ns].solver.nvars ;
    TS->u_size_total += TS->u_sizes[ns];
  }

  TS->u   = (double*) calloc (TS->u_size_total,sizeof(double));
  TS->rhs = (double*) calloc (TS->u_size_total,sizeof(double));
  _ArraySetValue_(TS->u  ,TS->u_size_total,0.0);
  _ArraySetValue_(TS->rhs,TS->u_size_total,0.0);

  /* initialize arrays to NULL, then allocate as necessary */
  TS->U             = NULL;
  TS->Udot          = NULL;
  TS->BoundaryFlux  = NULL;

  TS->bf_offsets = (int*) calloc (nsims, sizeof(int));
  TS->bf_sizes = (int*) calloc (nsims, sizeof(int));
  TS->bf_size_total = 0;
  for (ns = 0; ns < nsims; ns++) {
    TS->bf_offsets[ns] = TS->bf_size_total;
    TS->bf_sizes[ns] =  2 * sim[ns].solver.ndims * sim[ns].solver.nvars;
    TS->bf_size_total += TS->bf_sizes[ns];
  }

#if defined(HAVE_CUDA)
  if (sim[0].solver.use_gpu) {

    if (!strcmp(sim[0].solver.time_scheme,_RK_)) {

      /* explicit Runge-Kutta methods */
      ExplicitRKParameters  *params = (ExplicitRKParameters*)  sim[0].solver.msti;
      int nstages = params->nstages;

      gpuMalloc((void**)&TS->gpu_U, TS->u_size_total*nstages*sizeof(double));
      gpuMalloc((void**)&TS->gpu_Udot, TS->u_size_total*nstages*sizeof(double));
      gpuMalloc((void**)&TS->gpu_BoundaryFlux, TS->bf_size_total*nstages*sizeof(double));
      gpuMemset(TS->gpu_U, 0, TS->u_size_total*nstages*sizeof(double));
      gpuMemset(TS->gpu_Udot, 0, TS->u_size_total*nstages*sizeof(double));
      gpuMemset(TS->gpu_BoundaryFlux, 0, TS->bf_size_total*nstages*sizeof(double));

    } else {

      fprintf(stderr,"ERROR in TimeInitialize(): %s is not yet implemented on GPUs.\n",
              sim[0].solver.time_scheme );
      return 1;

    }

  } else {
#endif

    if (!strcmp(sim[0].solver.time_scheme,_RK_)) {

      /* explicit Runge-Kutta methods */
      ExplicitRKParameters  *params = (ExplicitRKParameters*)  sim[0].solver.msti;
      int nstages = params->nstages;
      TS->U     = (double**) calloc (nstages,sizeof(double*));
      TS->Udot  = (double**) calloc (nstages,sizeof(double*));
      for (i = 0; i < nstages; i++) {
        TS->U[i]    = (double*) calloc (TS->u_size_total,sizeof(double));
        TS->Udot[i] = (double*) calloc (TS->u_size_total,sizeof(double));
      }

      TS->BoundaryFlux = (double**) calloc (nstages,sizeof(double*));
      for (i=0; i<nstages; i++) {
        TS->BoundaryFlux[i] = (double*) calloc (TS->bf_size_total,sizeof(double));
      }

    } else if (!strcmp(sim[0].solver.time_scheme,_FORWARD_EULER_)) {

      int nstages = 1;
      TS->BoundaryFlux = (double**) calloc (nstages,sizeof(double*));
      for (i=0; i<nstages; i++) {
        TS->BoundaryFlux[i] = (double*) calloc (TS->bf_size_total,sizeof(double));
      }

    } else if (!strcmp(sim[0].solver.time_scheme,_GLM_GEE_)) {

      /* General Linear Methods with Global Error Estimate */
      GLMGEEParameters *params = (GLMGEEParameters*) sim[0].solver.msti;
      int nstages = params->nstages;
      int r       = params->r;
      TS->U     = (double**) calloc (2*r-1  ,sizeof(double*));
      TS->Udot  = (double**) calloc (nstages,sizeof(double*));
      for (i=0; i<2*r-1; i++)   TS->U[i]    = (double*) calloc (TS->u_size_total,sizeof(double));
      for (i=0; i<nstages; i++) TS->Udot[i] = (double*) calloc (TS->u_size_total,sizeof(double));

      TS->BoundaryFlux = (double**) calloc (nstages,sizeof(double*));
      for (i=0; i<nstages; i++) {
        TS->BoundaryFlux[i] = (double*) calloc (TS->bf_size_total,sizeof(double));
      }


      if (!strcmp(params->ee_mode,_GLM_GEE_YYT_)) {
        for (ns = 0; ns < nsims; ns++) {
          for (i=0; i<r-1; i++) {
            _ArrayCopy1D_(  (sim[ns].solver.u),
                            (TS->U[r+i] + TS->u_offsets[ns]),
                            (TS->u_sizes[ns])  );
          }
        }
      } else {
        for (i=0; i<r-1; i++) {
          _ArraySetValue_(TS->U[r+i],TS->u_size_total,0.0);
        }
      }
    }
#if defined(HAVE_CUDA)
  }
#endif

  /* set right-hand side function pointer */
  TS->RHSFunction = TimeRHSFunctionExplicit;

  /* open files for writing */
  if (!rank) {
    if (sim[0].solver.write_residual) TS->ResidualFile = (void*) fopen("residual.out","w");
    else                              TS->ResidualFile = NULL;
  } else                              TS->ResidualFile = NULL;

  for (ns = 0; ns < nsims; ns++) {
    sim[ns].solver.time_integrator = TS;
  }

  return 0;
}


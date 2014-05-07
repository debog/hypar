#include <basic.h>
#include <arrayfunctions.h>
#include <mpivars.h>
#include <hypar.h>
#include <timeintegration.h>

int TimeRK(void *ts)
{
  TimeIntegration *TS     = (TimeIntegration*) ts;
  HyPar           *solver = (HyPar*)           TS->solver;
  MPIVariables    *mpi    = (MPIVariables*)    TS->mpi;
  MSTIParameters  *params = (MSTIParameters*)  solver->msti;
  int             d, stage, i;
  _DECLARE_IERR_;

  int size = 1;
  for (d=0; d<solver->ndims; d++) size *= (solver->dim_local[d]+2*solver->ghosts);

  if (solver->PreStep)  { IERR solver->PreStep(solver->u,solver,mpi,TS->waqt); CHECKERR(ierr); }

  /* Calculate stage values */
  for (stage = 0; stage < params->nstages; stage++) {
    double stagetime = TS->waqt + params->c[stage]*TS->dt;
    if (solver->PreStage) { 
      if (stage) { IERR solver->PreStage(stage,TS->U ,solver,mpi,stagetime); CHECKERR(ierr); }
      else       { IERR solver->PreStage(1,&solver->u,solver,mpi,stagetime); CHECKERR(ierr); }
    }
    _ArrayCopy1D_(solver->u,TS->U[stage],size*solver->nvars);
    for (i = 0; i < stage; i++) {
      _ArrayAXPY_(TS->Udot[i],solver->dt*params->A[stage*params->nstages+i],
                  TS->U[stage],size*solver->nvars); 
    }
    IERR TS->RHSFunction(TS->Udot[stage],TS->U[stage],solver,mpi,stagetime);
    if (solver->PostStage) 
      { IERR solver->PostStage(stage,TS->U,solver,mpi,stagetime); CHECKERR(ierr); }

    _ArraySetValue_(TS->BoundaryFlux[stage],2*solver->ndims*solver->nvars,0.0);
    _ArrayCopy1D_(solver->StageBoundaryIntegral,TS->BoundaryFlux[stage],2*solver->ndims*solver->nvars);
  }

  /* Step completion */
  for (stage = 0; stage < params->nstages; stage++) {
    _ArrayAXPY_(TS->Udot[stage],solver->dt*params->b[stage],solver->u,size*solver->nvars);
    _ArrayAXPY_(TS->BoundaryFlux[stage],solver->dt*params->b[stage],solver->StepBoundaryIntegral,
                2*solver->ndims*solver->nvars);
  }

  if (solver->PostStep)  { IERR solver->PostStep(solver->u,solver,mpi,TS->waqt); CHECKERR(ierr); }

  return(0);
}


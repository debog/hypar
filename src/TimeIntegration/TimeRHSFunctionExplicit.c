/*! @file TimeRHSFunctionExplicit.c
    @brief Right-hand-side computation for explicit time integration
    @author Debojyoti Ghosh, Youngdae Kim
*/

#include <basic.h>
#if defined(HAVE_CUDA)
#include <arrayfunctions_gpu.h>
#else
#include <arrayfunctions.h>
#endif
#include <mpivars.h>
#include <hypar.h>

#include <time.h>

/*!
  This function computes the right-hand-side of the ODE given by
  \f{equation}{
    \frac {{\bf u}}{dt} = {\bf F}\left({\bf u}\right)
  \f}
  for explicit time integration methods, i.e., where
  \f{equation}{
    {\bf F}\left({\bf u}\right) = - {\bf F}_{\rm hyperbolic}\left({\bf u}\right)
                                  + {\bf F}_{\rm parabolic} \left({\bf u}\right)
                                  + {\bf F}_{\rm source}    \left({\bf u}\right),
  \f}
  given the solution \f${\bf u}\f$ and the current simulation time.
*/
int TimeRHSFunctionExplicit(
                              double  *rhs, /*!< Array to hold the computed right-hand-side */
                              double  *u,   /*!< Array holding the solution */
                              void    *s,   /*!< Solver object of type #HyPar */
                              void    *m,   /*!< MPI object of type #MPIVariables */
                              double  t     /*!< Current simulation time */
                           )
{
  HyPar           *solver = (HyPar*)        s;
  MPIVariables    *mpi    = (MPIVariables*) m;
  int             d;

  int size = 1;
  for (d=0; d<solver->ndims; d++) size *= (solver->dim_local[d]+2*solver->ghosts);

  /* apply boundary conditions and exchange data over MPI interfaces */
  solver->ApplyBoundaryConditions(solver,mpi,u,NULL,t);
  solver->ApplyIBConditions(solver,mpi,u,t);

  /* Evaluate hyperbolic, parabolic and source terms  and the RHS */
#if defined(HAVE_CUDA)
  if (solver->use_gpu) {
    gpuMPIExchangeBoundariesnD( solver->ndims,
                                solver->nvars,
                                solver->gpu_dim_local,
                                solver->ghosts,
                                mpi,
                                u );
  } else {
#endif
    MPIExchangeBoundariesnD(  solver->ndims,
                              solver->nvars,
                              solver->dim_local,
                              solver->ghosts,
                              mpi,
                              u);
#if defined(HAVE_CUDA)
  }
#endif

    solver->HyperbolicFunction( solver->hyp,
                                u,
                                solver,
                                mpi,
                                t,
                                1,
                                solver->FFunction,
                                solver->Upwind );
    solver->ParabolicFunction(solver->par,u,solver,mpi,t);
    solver->SourceFunction(solver->source,u,solver,mpi,t);

#if defined(HAVE_CUDA)
  if (solver->use_gpu) {
    gpuArraySetValue(rhs, size*solver->nvars, 0.0);
    gpuArrayAXPY(solver->hyp,    -1.0, rhs, size*solver->nvars);
    gpuArrayAXPY(solver->par,     1.0, rhs, size*solver->nvars);
    gpuArrayAXPY(solver->source,  1.0, rhs, size*solver->nvars);
  } else {
#endif
    _ArraySetValue_(rhs,size*solver->nvars,0.0);
    _ArrayAXPY_(solver->hyp   ,-1.0,rhs,size*solver->nvars);
    _ArrayAXPY_(solver->par   , 1.0,rhs,size*solver->nvars);
    _ArrayAXPY_(solver->source, 1.0,rhs,size*solver->nvars);
#if defined(HAVE_CUDA)
  }
#endif

  return(0);
}

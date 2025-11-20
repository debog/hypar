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
                              double  *a_rhs, /*!< Array to hold the computed right-hand-side */
                              double  *a_u,   /*!< Array holding the solution */
                              void    *a_s,   /*!< Solver object of type #HyPar */
                              void    *a_m,   /*!< MPI object of type #MPIVariables */
                              double  a_t     /*!< Current simulation time */
                           )
{
  HyPar           *solver = (HyPar*)        a_s;
  MPIVariables    *mpi    = (MPIVariables*) a_m;
  int             d;

  int size = 1;
  for (d=0; d<solver->m_ndims; d++) size *= (solver->m_dim_local[d]+2*solver->m_ghosts);

  /* apply boundary conditions and exchange data over MPI interfaces */
  solver->ApplyBoundaryConditions(solver,mpi,a_u,NULL,a_t);
  solver->ApplyIBConditions(solver,mpi,a_u,a_t);

  /* Evaluate hyperbolic, parabolic and source terms  and the RHS */
#if defined(HAVE_CUDA)
  if (solver->m_use_gpu) {
    gpuMPIExchangeBoundariesnD( solver->m_ndims,
                                solver->m_nvars,
                                solver->m_gpu_dim_local,
                                solver->m_ghosts,
                                mpi,
                                a_u );
  } else {
#endif
    MPIExchangeBoundariesnD(  solver->m_ndims,
                              solver->m_nvars,
                              solver->m_dim_local,
                              solver->m_ghosts,
                              mpi,
                              a_u);
#if defined(HAVE_CUDA)
  }
#endif

    solver->HyperbolicFunction( solver->m_hyp,
                                a_u,
                                solver,
                                mpi,
                                a_t,
                                1,
                                solver->FFunction,
                                solver->Upwind );
    solver->ParabolicFunction(solver->m_par,a_u,solver,mpi,a_t);
    solver->SourceFunction(solver->m_source,a_u,solver,mpi,a_t);

#if defined(HAVE_CUDA)
  if (solver->m_use_gpu) {
    gpuArraySetValue(a_rhs, size*solver->m_nvars, 0.0);
    gpuArrayAXPY(solver->m_hyp,    -1.0, a_rhs, size*solver->m_nvars);
    gpuArrayAXPY(solver->m_par,     1.0, a_rhs, size*solver->m_nvars);
    gpuArrayAXPY(solver->m_source,  1.0, a_rhs, size*solver->m_nvars);
  } else {
#endif
    _ArraySetValue_(a_rhs,size*solver->m_nvars,0.0);
    _ArrayAXPY_(solver->m_hyp   ,-1.0,a_rhs,size*solver->m_nvars);
    _ArrayAXPY_(solver->m_par   , 1.0,a_rhs,size*solver->m_nvars);
    _ArrayAXPY_(solver->m_source, 1.0,a_rhs,size*solver->m_nvars);
#if defined(HAVE_CUDA)
  }
#endif

  return(0);
}

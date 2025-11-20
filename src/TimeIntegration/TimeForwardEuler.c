/*! @file TimeForwardEuler.c
    @brief Forward Euler method
    @author Debojyoti Ghosh
*/

#include <basic.h>
#include <arrayfunctions.h>
#include <simulation_object.h>
#include <timeintegration.h>

/*!
  Advance the ODE given by
  \f{equation}{
    \frac{d{\bf u}}{dt} = {\bf F} \left({\bf u}\right)
  \f}
  by one time step of size #HyPar::m_dt using the forward Euler method
  given by
  \f{equation}{
    {\bf u}^{n+1} = {\bf u}^n + \Delta t {\bf F}\left( {\bf u}^n \right)
  \f}
  where the superscript represents the time level, \f$\Delta t\f$ is the
  time step size #HyPar::m_dt, and \f${\bf F}\left({\bf u}\right)\f$ is
  computed by #TimeIntegration::RHSFunction.
*/
int TimeForwardEuler(
                      void *ts /*!< Time integrator object of type #TimeIntegration */
                    )
{
  TimeIntegration* TS = (TimeIntegration*) ts;
  SimulationObject* sim = (SimulationObject*) TS->m_simulation;
  int ns, nsims = TS->m_nsims;
  _DECLARE_IERR_;

  for (ns = 0; ns < nsims; ns++) {

    HyPar* solver = &(sim[ns].solver);
    MPIVariables* mpi = &(sim[ns].mpi);

    double* rhs = TS->m_rhs + TS->m_u_offsets[ns];

    /* Evaluate right-hand side */
    IERR TS->RHSFunction( rhs,
                          solver->m_u,
                          solver,
                          mpi,
                          TS->m_waqt  );   CHECKERR(ierr);

    /* update solution */
    _ArrayAXPY_(  rhs,
                  TS->m_dt,
                  solver->m_u,
                  TS->m_u_sizes[ns] );

    _ArrayScaleCopy1D_( solver->m_stage_boundary_integral,
                        TS->m_dt,
                        solver->m_step_boundary_integral,
                        TS->m_bf_sizes[ns] );
  }

  return(0);
}

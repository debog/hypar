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
  by one time step of size #HyPar::dt using the forward Euler method
  given by
  \f{equation}{
    {\bf u}^{n+1} = {\bf u}^n + \Delta t {\bf F}\left( {\bf u}^n \right)
  \f}
  where the superscript represents the time level, \f$\Delta t\f$ is the
  time step size #HyPar::dt, and \f${\bf F}\left({\bf u}\right)\f$ is
  computed by #TimeIntegration::RHSFunction.
*/
int TimeForwardEuler(
                      void *ts /*!< Time integrator object of type #TimeIntegration */
                    )
{
  TimeIntegration* TS = (TimeIntegration*) ts;
  SimulationObject* sim = (SimulationObject*) TS->simulation;
  int ns, nsims = TS->nsims;
  _DECLARE_IERR_;

  for (ns = 0; ns < nsims; ns++) {

    HyPar* solver = &(sim[ns].solver);
    MPIVariables* mpi = &(sim[ns].mpi);

    double* rhs = TS->rhs + TS->u_offsets[ns];

    /* Evaluate right-hand side */
    IERR TS->RHSFunction( rhs,
                          solver->u,
                          solver,
                          mpi,
                          TS->waqt  );   CHECKERR(ierr);

    /* update solution */
    _ArrayAXPY_(  rhs,
                  TS->dt,
                  solver->u,
                  TS->u_sizes[ns] );

    _ArrayScaleCopy1D_( solver->StageBoundaryIntegral,
                        TS->dt,
                        solver->StepBoundaryIntegral,
                        TS->bf_sizes[ns] );
  }

  return(0);
}

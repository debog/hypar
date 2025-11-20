/*! @file TimeGLMGEE.c
    @brief General Linear Methods with Global Error Estimators
    @author Debojyoti Ghosh
*/

#include <string.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <simulation_object.h>
#include <timeintegration.h>

/*!
  Advance the ODE given by
  \f{equation}{
    \frac{d{\bf u}}{dt} = {\bf F} \left({\bf u}\right)
  \f}
  by one time step of size #HyPar::m_dt using the \f$s\f$-stage General Linear Method with
  Global Error Estimation (GLM-GEE), given by
  \f{align}{
    {\bf U}^{\left(i\right)} &= c_{i0}{\bf u}_n + \sum_{j=1}^{r-1} c_{ij} \tilde{\bf u}_n^j
                               + \Delta t \sum_{j=1}^{i-1} a_{ij} {\bf F}\left({\bf U}^{\left(j\right)}\right), i=1,\cdots,s\\
    {\bf u}_{n+1} &= d_{00} {\bf u}_n + \sum_{j=1}^{r-1} d_{0j} \tilde{\bf u}_n^j
                      + \Delta t \sum_{j=1}^s b_{0j} {\bf F}\left({\bf U}^{\left(j\right)}\right), \\
    \tilde{\bf u}_{n+1}^i &= d_{i0} {\bf u}_n + \sum_{j=1}^{r-1} d_{ij} \tilde{\bf u}_n^j
                      + \Delta t \sum_{j=1}^s b_{ij} {\bf F}\left({\bf U}^{\left(j\right)}\right), i=1,\cdots,r-1
  \f}
  where the superscripts in parentheses represent stages, the subscripts represent the time level, the
  superscripts without parentheses represent auxiliary solutions, \f$\Delta t\f$ is the
  time step size #HyPar::m_dt, \f$\tilde{\bf u}^i, i=1,\cdots,r-1\f$ are the auxiliary solutions,
  and \f${\bf F}\left({\bf u}\right)\f$ is computed by #TimeIntegration::RHSFunction. The coefficients
  defining this methods are:
  + \f$a_{ij}, i=1,\cdots,s, j=1,\cdots,s\f$ (#GLMGEEParameters::A)
  + \f$b_{ij}, i=0,\cdots,r-1, j=1,\cdots,s\f$ (#GLMGEEParameters::B)
  + \f$c_{ij}, i=1,\cdots,s, j=0,\cdots,r-1\f$ (#GLMGEEParameters::C)
  + \f$d_{ij}, i=0,\cdots,r-1, j=0,\cdots,r-1\f$ (#GLMGEEParameters::D)

  where \f$s\f$ is the number of stages (#GLMGEEParameters::nstages) and \f$r\f$ is the number of auxiliary solutions
  propagated with the solution (#GLMGEEParameters::r).

  Note: In the code #TimeIntegration::Udot is equivalent to \f${\bf F}\left({\bf u}\right)\f$.

  References:
  + Constantinescu, E. M., "Estimating Global Errors in Time Stepping.", Submitted, 2015 (http://arxiv.org/abs/1503.05166).
*/
int TimeGLMGEE(void *ts /*!< Object of type #TimeIntegration */)
{
  TimeIntegration* TS = (TimeIntegration*) ts;
  SimulationObject* sim = (SimulationObject*) TS->m_simulation;
  GLMGEEParameters* params = (GLMGEEParameters*) sim[0].solver.m_msti;
  int ns, stage, i, j, nsims = TS->m_nsims;
  _DECLARE_IERR_;

  int    s  = params->nstages;
  int    r  = params->r;
  double dt = TS->m_dt;
  double  *A =params->A,
          *B =params->B,
          *C=params->C,
          *D=params->D,
          *c=params->c,
          **U = TS->m_U,
          **Udot = TS->m_Udot,
          **Uaux = &TS->m_U[r];

  /* Calculate stage values */
  for (j=0; j<s; j++) {

    double stagetime = TS->m_waqt + c[j]*dt;

    for (ns = 0; ns < nsims; ns++) {
      _ArrayScaleCopy1D_( sim[ns].solver.m_u,
                          C[j*r+0],
                          (U[0] + TS->m_u_offsets[ns]),
                          (TS->m_u_sizes[ns]) );
    }

    for (i=1;i<r;i++) _ArrayAXPY_(Uaux[i-1], C[j*r+i]   , U[0], TS->m_u_size_total);
    for (i=0;i<j;i++) _ArrayAXPY_(Udot[i]  , dt*A[j*s+i], U[0], TS->m_u_size_total);

    for (ns = 0; ns < nsims; ns++) {
      if (sim[ns].solver.PreStage) {
        fprintf(stderr,"Call to solver->PreStage() commented out in TimeGLMGEE()!\n");
        return 1;
//        IERR sim[ns].solver.PreStage( stage,
//                                      (U),
//                                      &(sim[ns].solver),
//                                      &(sim[ns].mpi),
//                                      stagetime ); CHECKERR(ierr);
      }
    }

    for (ns = 0; ns < nsims; ns++) {
      IERR TS->RHSFunction( (Udot[j] + TS->m_u_offsets[ns]),
                            (U[0] + TS->m_u_offsets[ns]),
                            &(sim[ns].solver),
                            &(sim[ns].mpi),
                            stagetime);
    }

    for (ns = 0; ns < nsims; ns++) {
      if (sim[ns].solver.PostStage) {
        IERR sim[ns].solver.PostStage(  (U[j] + TS->m_u_offsets[ns]),
                                        &(sim[ns].solver),
                                        &(sim[ns].mpi),
                                        stagetime); CHECKERR(ierr);
      }
    }

    _ArraySetValue_(TS->m_BoundaryFlux[j], TS->m_bf_size_total, 0.0);
    for (ns = 0; ns < nsims; ns++) {
      _ArrayCopy1D_(  sim[ns].solver.m_stage_boundary_integral,
                      (TS->m_BoundaryFlux[j] + TS->m_bf_offsets[ns]),
                      TS->m_bf_sizes[ns] );
    }
  }

  /* Step completion */
  for (j=0; j<r; j++) {

    for (ns = 0; ns < nsims; ns++) {

      _ArrayScaleCopy1D_( sim[ns].solver.m_u,
                          D[j*r+0],
                          (U[j]+TS->m_u_offsets[ns]),
                          (TS->m_u_sizes[ns]) );

    }

    for (i=1; i<r; i++) _ArrayAXPY_(Uaux[i-1], D[j*r+i]   , U[j], TS->m_u_size_total);
    for (i=0; i<s; i++) _ArrayAXPY_(Udot[i]  , dt*B[j*s+i], U[j], TS->m_u_size_total);

  }

  for (ns = 0; ns < nsims; ns++) {
    for (i=0; i<s; i++) {
      _ArrayAXPY_(  (TS->m_BoundaryFlux[i] + TS->m_bf_offsets[ns]),
                    dt*B[0*s+i],
                    sim[ns].solver.m_step_boundary_integral,
                    TS->m_bf_sizes[ns] );
    }

  }

  for (ns = 0; ns < nsims; ns++) {
    _ArrayCopy1D_(  (U[0] + TS->m_u_offsets[ns]),
                    sim[ns].solver.m_u,
                    TS->m_u_sizes[ns] );
  }
  for (i=1; i<r; i++) {
    _ArrayCopy1D_(U[i],Uaux[i-1],TS->m_u_size_total);
  }

  return(0);
}


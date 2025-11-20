/*! @file VlasovFunctions.c
    @author Debojyoti Ghosh
    @brief Misc. functions for the Vlasov equations
*/

#include <physicalmodels/vlasov.h>
#include <hypar.h>
#include <common.h>

int VlasovEField(double*, void*, double);

/*! Vlasov-specific function called at the beginning of each time-step:
    Calls the function to set the electric field
*/
int VlasovPreStep( double  *u,   /*!< Solution (conserved variables) */
                   void    *s,   /*!< Solver object of type #HyPar */
                   void    *m,   /*!< MPI object of type #MPIVariables */
                   double  waqt  /*!< Current solution time */
                 )
{
  HyPar  *solver = (HyPar*) s;
  Vlasov *param  = (Vlasov*) solver->m_physics;

  if (param->m_use_log_form) takeExp(u,solver->m_npoints_local_wghosts);
  int ierr = VlasovEField(u, solver, waqt);
  if (param->m_use_log_form) takeLog(u,solver->m_npoints_local_wghosts);
  if (ierr) return ierr;

  return 0;
}

/*! Vlasov-specific function called at the end of each stage in
    a multistage time integrator:
    Calls the function to set the electric field
*/
int VlasovPostStage( double  *u,   /*!< Solution (conserved variables) */
                     void    *s,   /*!< Solver object of type #HyPar */
                     void    *m,   /*!< MPI object of type #MPIVariables */
                     double  waqt  /*!< Current solution time */
                   )
{
  HyPar  *solver = (HyPar*) s;
  Vlasov *param  = (Vlasov*) solver->m_physics;

  if (param->m_use_log_form) takeExp(u,solver->m_npoints_local_wghosts);
  int ierr = VlasovEField(u, solver, waqt);
  if (param->m_use_log_form) takeLog(u,solver->m_npoints_local_wghosts);
  if (ierr) return ierr;

  return 0;
}

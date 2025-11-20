/*! @file PetscSetInitialGuessROM.cpp
 *  @brief Compute an initial guess for implicit solves using a reduced-order model
 *  @author Debojyoti Ghosh
*/

#ifdef with_petsc

#include <simulation_object.h>
#include <petscinterface.h>
#ifdef with_librom
#include <librom_interface.h>
#endif

#ifdef with_librom

/*! Compute the initial guess for a nonlinear solve
 *  using a trained libROM reduced-order model. */
PetscErrorCode PetscSetInitialGuessROM( SNES  snes, /*!< Nonlinear solver object (see PETSc SNES) */
                                        Vec   X,    /*!< Initial guess vector */
                                        void* ctxt  /*!< Application context */ )
{
  PETScContext* context = (PETScContext*) ctxt;
  SimulationObject* sim = (SimulationObject*) context->m_simobj;
  int nsims = context->m_nsims;

  PetscFunctionBegin;

  double stage_time = context->m_t_start;
  if (context->m_stage_times.size() > (context->m_stage_index+1)) {
    stage_time += context->m_stage_times[context->m_stage_index+1] * context->m_dt;
  } else {
    stage_time += context->m_dt;
  }

  ((libROMInterface*)context->m_rom_interface)->predict(sim, stage_time);
  if (!context->m_rank) {
    printf(  "libROM: Predicted ROM intial guess (time %1.4e), wallclock time: %f.\n",
             stage_time, ((libROMInterface*)context->m_rom_interface)->predictWallclockTime() );
  }

  for (int ns = 0; ns < nsims; ns++) {
    TransferVecToPETSc(sim[ns].solver.m_u,X,context,ns,context->m_offsets[ns]);
  }


  PetscFunctionReturn(0);
}

#endif

#endif

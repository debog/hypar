/*! @file PetscCleanup.c
    @author Debojyoti Ghosh
    @brief Cleans up allocations in PETSc interface.
*/

#ifdef with_petsc

#include <stdlib.h>
#include <petscinterface_struct.h>

/*! Clean up allocations in the PETSc interface */
int PetscCleanup(void *obj /*!< Object of type #PETScContext */)
{
  PETScContext *ctxt = (PETScContext*) obj;
  for (int i = 0; i < ctxt->m_globalDOF.size(); i++) {
    free(ctxt->m_globalDOF[i]);
  }
  ctxt->m_globalDOF.clear();
  for (int i = 0; i < ctxt->m_points.size(); i++) {
    free(ctxt->m_points[i]);
  }
  ctxt->m_points.clear();
  if (ctxt->m_offsets) free(ctxt->m_offsets);
  return(0);
}

#endif

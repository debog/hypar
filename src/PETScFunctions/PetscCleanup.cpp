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
  for (int i = 0; i < ctxt->globalDOF.size(); i++) {
    free(ctxt->globalDOF[i]);
  }
  ctxt->globalDOF.clear();
  for (int i = 0; i < ctxt->points.size(); i++) {
    free(ctxt->points[i]);
  }
  ctxt->points.clear();
  if (ctxt->offsets) free(ctxt->offsets);
  return(0);
}

#endif

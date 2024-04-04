/*! @file TransferFromPETSc.cpp
    @brief Copy from PETSc vector to HyPar array
    @author Debojyoti Ghosh
*/

#ifdef with_petsc

#include <stdlib.h>
#include <vector>
#include <basic.h>
#include <arrayfunctions.h>
#include <simulation_object.h>
#include <petscinterface_struct.h>

#undef __FUNCT__
#define __FUNCT__ "TransferVecFromPETSc"

/*! Copy data from a PETSc vector (used by PETSc time integrators, and with no
    ghost points) to a HyPar::u array (with ghost points).

    \sa TransferVecToPETSc()
*/
int TransferVecFromPETSc( double* const u, /*!< HyPar::u type array (with ghost points) */
                          const Vec Y, /*!< PETSc vector */
                          void* ctxt, /*!< Object of type #PETScContext */
                          const int sim_idx,/*!< Simulation object index */
                          const int offset  /*!< Offset */ )
{
  PETScContext* context = (PETScContext*) ctxt;
  SimulationObject* sim = (SimulationObject*) context->simobj;
  const double* Yarr;

  PetscFunctionBegin;
  VecGetArrayRead(Y,&Yarr);
  std::vector<int> index(sim[sim_idx].solver.ndims,0);
  ArrayCopynD(  sim[sim_idx].solver.ndims,
                (Yarr+offset),
                u,
                sim[sim_idx].solver.dim_local,
                0,
                sim[sim_idx].solver.ghosts,
                index.data(),
                sim[sim_idx].solver.nvars );
  VecRestoreArrayRead(Y,&Yarr);

  PetscFunctionReturn(0);
}

#endif

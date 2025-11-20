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
int TransferVecFromPETSc( double* const a_u, /*!< HyPar::u type array (with ghost points) */
                          const Vec a_Y, /*!< PETSc vector */
                          void* a_ctxt, /*!< Object of type #PETScContext */
                          const int a_sim_idx,/*!< Simulation object index */
                          const int a_offset  /*!< Offset */ )
{
  PETScContext* context = (PETScContext*) a_ctxt;
  SimulationObject* sim = (SimulationObject*) context->m_simobj;
  const double* Yarr;

  PetscFunctionBegin;
  VecGetArrayRead(a_Y,&Yarr);
  std::vector<int> index(sim[a_sim_idx].solver.m_ndims,0);
  ArrayCopynD(  sim[a_sim_idx].solver.m_ndims,
                (Yarr+a_offset),
                a_u,
                sim[a_sim_idx].solver.m_dim_local,
                0,
                sim[a_sim_idx].solver.m_ghosts,
                index.data(),
                sim[a_sim_idx].solver.m_nvars );
  VecRestoreArrayRead(a_Y,&Yarr);

  PetscFunctionReturn(0);
}

#endif

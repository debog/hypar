/*! @file TransferToPETSc.cpp
    @brief Copy from a HyPar array to a PETSc vector.
    @author Debojyoti Ghosh
*/

#ifdef with_petsc

#include <stdlib.h>
#include <vector>
#include <basic.h>
#include <arrayfunctions.h>
#include <simulation_object.h>
#include <bandedmatrix.h>
#include <petscinterface_struct.h>

#undef __FUNCT__
#define __FUNCT__ "TransferVecToPETSc"

/*! Copy data to a PETSc vector (used by PETSc time integrators, and with no
    ghost points) from a HyPar::u array (with ghost points).

    \sa TransferVecFromPETSc()
*/
int TransferVecToPETSc( const double* const a_u, /*!< HyPar::u type array (with ghost points) */
                        Vec a_Y, /*!< PETSc vector */
                        void* a_ctxt, /*!< Object of type #PETScContext */
                        const int a_sim_idx,/*!< Simulation object index */
                        const int a_offset  /*!< Offset */ )
{
  PETScContext* context = (PETScContext*) a_ctxt;
  SimulationObject* sim = (SimulationObject*) context->m_simobj;
  double* Yarr;

  PetscFunctionBegin;
  VecGetArray(a_Y,&Yarr);
  std::vector<int> index(sim[a_sim_idx].solver.m_ndims,0);
  ArrayCopynD(  sim[a_sim_idx].solver.m_ndims,
                a_u,
                (Yarr+a_offset),
                sim[a_sim_idx].solver.m_dim_local,
                sim[a_sim_idx].solver.m_ghosts,
                0,
                index.data(),
                sim[a_sim_idx].solver.m_nvars );
  VecRestoreArray(a_Y,&Yarr);

  PetscFunctionReturn(0);
}

/*!
  Copy a matrix of type #BandedMatrix to a PETSc matrix.
*/
int TransferMatToPETSc( void *a_J,    /*!< Matrix of type #BandedMatrix */
                        Mat   a_A,    /*!< PETSc matrix */
                        void *a_ctxt  /*!< Object of type #PETScContext */ )
{
  BandedMatrix    *M = (BandedMatrix*) a_J;
  PetscErrorCode  ierr     = 0;
  int             bs = M->m_BlockSize, nbands = M->m_nbands, bs2 = bs*bs;

  for (int i=0; i<M->m_nrows_local; i++) {
    int     colind[nbands];
    double  val[bs][bs*nbands];
    for (int n=0; n<nbands; n++) {
      colind[n] = M->m_ncol[nbands*i+n];
      for (int p=0; p<bs; p++) {
        for (int q = 0; q<bs; q++) {
          val[p][n*bs+q] = M->m_data[i*nbands*bs2+n*bs2+p*bs+q];
        }
      }
    }
    MatSetValuesBlocked(a_A,1,&M->m_nrow[i],M->m_nbands,&colind[0],&val[0][0],INSERT_VALUES);
  }

  MatAssemblyBegin(a_A,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd  (a_A,MAT_FINAL_ASSEMBLY);

  return(0);
}

#endif

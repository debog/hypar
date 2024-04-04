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
int TransferVecToPETSc( const double* const u, /*!< HyPar::u type array (with ghost points) */
                        Vec Y, /*!< PETSc vector */
                        void* ctxt, /*!< Object of type #PETScContext */
                        const int sim_idx,/*!< Simulation object index */
                        const int offset  /*!< Offset */ )
{
  PETScContext* context = (PETScContext*) ctxt;
  SimulationObject* sim = (SimulationObject*) context->simobj;
  double* Yarr;

  PetscFunctionBegin;
  VecGetArray(Y,&Yarr);
  std::vector<int> index(sim[sim_idx].solver.ndims,0);
  ArrayCopynD(  sim[sim_idx].solver.ndims,
                u,
                (Yarr+offset),
                sim[sim_idx].solver.dim_local,
                sim[sim_idx].solver.ghosts,
                0,
                index.data(),
                sim[sim_idx].solver.nvars );
  VecRestoreArray(Y,&Yarr);

  PetscFunctionReturn(0);
}

/*!
  Copy a matrix of type #BandedMatrix to a PETSc matrix.
*/
int TransferMatToPETSc( void *J,    /*!< Matrix of type #BandedMatrix */
                        Mat   A,    /*!< PETSc matrix */
                        void *ctxt  /*!< Object of type #PETScContext */ )
{
  BandedMatrix    *M = (BandedMatrix*) J;
  PetscErrorCode  ierr     = 0;
  int             bs = M->BlockSize, nbands = M->nbands, bs2 = bs*bs;

  for (int i=0; i<M->nrows_local; i++) {
    int     colind[nbands];
    double  val[bs][bs*nbands];
    for (int n=0; n<nbands; n++) {
      colind[n] = M->ncol[nbands*i+n];
      for (int p=0; p<bs; p++) {
        for (int q = 0; q<bs; q++) {
          val[p][n*bs+q] = M->data[i*nbands*bs2+n*bs2+p*bs+q];
        }
      }
    }
    MatSetValuesBlocked(A,1,&M->nrow[i],M->nbands,&colind[0],&val[0][0],INSERT_VALUES);
  }

  MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd  (A,MAT_FINAL_ASSEMBLY);

  return(0);
}

#endif

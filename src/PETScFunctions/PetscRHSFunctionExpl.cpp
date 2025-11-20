/*! @file PetscRHSFunctionExpl.cpp
    @brief Function to compute the right-hand-side for explicit time integration
    @author Debojyoti Ghosh
*/

#ifdef with_petsc

#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mpivars_cpp.h>
#include <simulation_object.h>
#include <petscinterface.h>

#undef __FUNCT__
#define __FUNCT__ "PetscRHSFunctionExpl"

/*!
  Compute the right-hand-side (RHS) for the explicit time integration of the
  governing equations: The spatially discretized ODE can be expressed as
  \f{equation}{
    \frac {d{\bf U}} {dt} = {\bf F}\left({\bf U}\right).
  \f}
  This function computes \f${\bf F}\left({\bf U}\right)\f$, given \f${\bf U}\f$.

  \sa TSSetRHSFunction (https://petsc.org/release/docs/manualpages/TS/TSSetRHSFunction.html)

  Note:
  + \f${\bf U}\f$ is \a Y in the code.
  + See https://petsc.org/release/docs/manualpages/TS/index.html for documentation on
    PETSc's time integrators.
  + All functions and variables whose names start with Vec, Mat, PC, KSP, SNES, and TS are defined by PETSc. Refer to
    the PETSc documentation (https://petsc.org/release/docs/). Usually, googling with the function
    or variable name yields the specific doc page dealing with that function/variable.
*/
PetscErrorCode PetscRHSFunctionExpl(  TS        a_ts,   /*!< Time integration object */
                                      PetscReal a_t,    /*!< Current simulation time */
                                      Vec       a_Y,    /*!< State vector (input) */
                                      Vec       a_F,    /*!< The computed right-hand-side vector */
                                      void *a_ctxt /*!< Object of type #PETScContext */ )
{
  PETScContext* context = (PETScContext*) a_ctxt;
  SimulationObject* sim = (SimulationObject*) context->m_simobj;
  int nsims = context->m_nsims;

  PetscFunctionBegin;

  for (int ns = 0; ns < nsims; ns++) {

    HyPar* solver = &(sim[ns].solver);
    MPIVariables* mpi = &(sim[ns].mpi);

    solver->m_count_rhs_function++;

    int size = solver->m_npoints_local_wghosts;
    double* u = solver->m_u;
    double* rhs = solver->m_rhs;

    /* copy solution from PETSc vector */
    TransferVecFromPETSc(u,a_Y,context,ns,context->m_offsets[ns]);
    /* apply boundary conditions and exchange data over MPI interfaces */
    solver->ApplyBoundaryConditions(solver,mpi,u,NULL,a_t);
    MPIExchangeBoundariesnD(  solver->m_ndims,
                              solver->m_nvars,
                              solver->m_dim_local,
                              solver->m_ghosts,
                              mpi,
                              u );

    /* Evaluate hyperbolic, parabolic and source terms  and the RHS */
    solver->HyperbolicFunction(solver->m_hyp,u,solver,mpi,a_t,1,solver->FFunction,solver->Upwind);

    solver->ParabolicFunction (solver->m_par,u,solver,mpi,a_t);
    solver->SourceFunction    (solver->m_source,u,solver,mpi,a_t);

    _ArraySetValue_(rhs,size*solver->m_nvars,0.0);
    _ArrayAXPY_(solver->m_hyp   ,-1.0,rhs,size*solver->m_nvars);
    _ArrayAXPY_(solver->m_par   , 1.0,rhs,size*solver->m_nvars);
    _ArrayAXPY_(solver->m_source, 1.0,rhs,size*solver->m_nvars);

    /* Transfer RHS to PETSc vector */
    TransferVecToPETSc(rhs,a_F,context,ns,context->m_offsets[ns]);

  }

  PetscFunctionReturn(0);
}

#endif

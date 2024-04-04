/*! @file PetscIFunctionImpl.cpp
    @brief Compute the right-hand-side for implicit time integration.
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
#define __FUNCT__ "PetscIFunctionImpl"

/*!
  Compute the left-hand-side for the implicit time integration of the
  governing equations: The spatially discretized ODE can be expressed as
  \f{equation}{
    \frac {d{\bf U}} {dt} = {\bf F}\left({\bf U}\right).
  \f}
  This function computes \f$\dot{\bf U} - {\bf F}\left({\bf U}\right)\f$,
  given \f${\bf U},\dot{\bf U}\f$.

  \sa TSSetIFunction (https://petsc.org/release/docs/manualpages/TS/TSSetIFunction.html)

  Note:
  + \f${\bf U}\f$ is \a Y in the code.
  + See https://petsc.org/release/docs/manualpages/TS/index.html for documentation on
    PETSc's time integrators.
  + All functions and variables whose names start with Vec, Mat, PC, KSP, SNES, and TS are defined by PETSc. Refer to
    the PETSc documentation (https://petsc.org/release/docs/). Usually, googling with the function
    or variable name yields the specific doc page dealing with that function/variable.
*/
PetscErrorCode PetscIFunctionImpl(  TS        ts,   /*!< Time integration object */
                                    PetscReal t,    /*!< Current simulation time */
                                    Vec       Y,    /*!< State vector (input) */
                                    Vec       Ydot, /*!< Time derivative of the state vector (input) */
                                    Vec       F,    /*!< The computed right-hand-side vector */
                                    void      *ctxt /*!< Object of type #PETScContext */ )
{
  PETScContext* context = (PETScContext*) ctxt;
  SimulationObject* sim = (SimulationObject*) context->simobj;
  int nsims = context->nsims;

  PetscFunctionBegin;
  for (int ns = 0; ns < nsims; ns++) {

    HyPar* solver = &(sim[ns].solver);
    MPIVariables* mpi = &(sim[ns].mpi);

    solver->count_RHSFunction++;

    int size = solver->npoints_local_wghosts;
    double* u = solver->u;
    double* rhs = solver->rhs;

    /* copy solution from PETSc vector */
    TransferVecFromPETSc(u,Y,context,ns,context->offsets[ns]);
    /* apply boundary conditions and exchange data over MPI interfaces */
    solver->ApplyBoundaryConditions(solver,mpi,u,NULL,t);
    MPIExchangeBoundariesnD(  solver->ndims,
                              solver->nvars,
                              solver->dim_local,
                              solver->ghosts,
                              mpi,
                              u );

    /* Evaluate hyperbolic, parabolic and source terms  and the RHS */
    solver->HyperbolicFunction(solver->hyp,u,solver,mpi,t,1,solver->FFunction,solver->Upwind);
    solver->ParabolicFunction (solver->par,u,solver,mpi,t);
    solver->SourceFunction    (solver->source,u,solver,mpi,t);

    _ArraySetValue_(rhs,size*solver->nvars,0.0);
    _ArrayAXPY_(solver->hyp   ,-1.0,rhs,size*solver->nvars);
    _ArrayAXPY_(solver->par   , 1.0,rhs,size*solver->nvars);
    _ArrayAXPY_(solver->source, 1.0,rhs,size*solver->nvars);

    /* save a copy of the solution and RHS for use in IJacobian */
    _ArrayCopy1D_(u  ,solver->uref  ,(size*solver->nvars));
    _ArrayCopy1D_(rhs,solver->rhsref,(size*solver->nvars));

    /* Transfer RHS to PETSc vector */
    TransferVecToPETSc(rhs,F,context,ns,context->offsets[ns]);

  }

  /* LHS = Ydot - F(u) */
  VecAYPX(F,-1.0,Ydot);

  PetscFunctionReturn(0);
}

#endif

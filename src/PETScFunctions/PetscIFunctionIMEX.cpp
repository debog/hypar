/*! @file PetscIFunctionIMEX.cpp
    @brief Compute the implicitly-treated part of the right-hand-side for IMEX time integration.
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
#define __FUNCT__ "PetscIFunctionIMEX"

/*!
  Compute the implicitly-treated part of the right-hand-side for the implicit-explicit (IMEX) time integration
  of the governing equations: The ODE, obtained after discretizing the governing PDE in space,
  is expressed as follows (for the purpose of IMEX time integration):
  \f{eqnarray}{
    \frac {d{\bf U}}{dt} &=& {\bf F}\left({\bf U}\right) + {\bf G}\left({\bf U}\right), \\
    \Rightarrow \dot{\bf U} - {\bf G}\left({\bf U}\right) &=& {\bf F}\left({\bf U}\right),
  \f}
  where \f${\bf F}\f$ is non-stiff and integrated in time explicitly, and \f${\bf G}\f$
  is stiff and integrated in time implicitly, and \f${\bf U}\f$ represents the entire
  solution vector (state vector).

    Note that \f${\bf G}\left({\bf U}\right)\f$ represents all terms that the user has indicated to be
    integrated in time implicitly (#PETScContext::flag_hyperbolic_f, #PETScContext::flag_hyperbolic_df,
    #PETScContext::flag_hyperbolic, #PETScContext::flag_parabolic, and #PETScContext::flag_source).

  This function computes the left-hand-side of the above equation:
  \f{equation}{
    \mathcal{G}\left(\dot{\bf U},{\bf U},t\right) = \dot{\bf U} - {\bf G}\left({\bf U}\right)
  \f}
  given \f$\dot{\bf U}\f$ and \f${\bf U}\f$.

  \sa PetscRHSFunctionIMEX()

  \b Notes:
  + \a Y and \a Ydot in the code are \f${\bf U}\f$ and \f$\dot{\bf U}\f$, respectively. PETsc denotes
    the state vector with \f${\bf Y}\f$ in its time integrators.
  + It is assumed that the reader is familiar with PETSc's implementation of IMEX time integrators, for
    example, TSARKIMEX (https://petsc.org/release/docs/manualpages/TS/TSARKIMEX.html).
  + See https://petsc.org/release/docs/manualpages/TS/index.html for documentation on
    PETSc's time integrators.
  + All functions and variables whose names start with Vec, Mat, PC, KSP, SNES, and TS are defined by PETSc. Refer to
    the PETSc documentation (https://petsc.org/release/docs/). Usually, googling with the function
    or variable name yields the specific doc page dealing with that function/variable.
*/
PetscErrorCode PetscIFunctionIMEX(  TS        ts,     /*!< The time integration object */
                                    PetscReal t,      /*!< Current solution time */
                                    Vec       Y,      /*!< State vector (input) */
                                    Vec       Ydot,   /*!< Time derivative of the state vector (input) */
                                    Vec       F,      /*!< The computed function vector */
                                    void      *ctxt   /*!< Object of type PETScContext */ )
{
  PETScContext* context = (PETScContext*) ctxt;
  SimulationObject* sim = (SimulationObject*) context->m_simobj;
  int nsims = context->m_nsims;

  PetscFunctionBegin;

  context->m_waqt = t;

  for (int ns = 0; ns < nsims; ns++) {

    HyPar* solver = &(sim[ns].solver);
    MPIVariables* mpi = &(sim[ns].mpi);

    solver->m_count_i_function++;

    int size = solver->m_npoints_local_wghosts;
    double *u = solver->m_u;
    double *rhs = solver->m_rhs;

    /* copy solution from PETSc vector */
    TransferVecFromPETSc(u,Y,context,ns,context->m_offsets[ns]);
    /* apply boundary conditions and exchange data over MPI interfaces */
    solver->ApplyBoundaryConditions(solver,mpi,u,NULL,t);
    MPIExchangeBoundariesnD(  solver->m_ndims,
                              solver->m_nvars,
                              solver->m_dim_local,
                              solver->m_ghosts,
                              mpi,
                              u );

    /* initialize right-hand side to zero */
    _ArraySetValue_(rhs,size*solver->m_nvars,0.0);

    /* Evaluate hyperbolic, parabolic and source terms  and the RHS */
    if ((!strcmp(solver->m_split_hyperbolic_flux,"yes")) && solver->m_flag_fdf_specified) {
      if (context->m_flag_hyperbolic_f == _IMPLICIT_) {
        solver->HyperbolicFunction(solver->m_hyp,u,solver,mpi,t,0,solver->FdFFunction,solver->UpwindFdF);
        _ArrayAXPY_(solver->m_hyp,-1.0,rhs,size*solver->m_nvars);
      }
      if (context->m_flag_hyperbolic_df == _IMPLICIT_) {
        solver->HyperbolicFunction(solver->m_hyp,u,solver,mpi,t,0,solver->dFFunction,solver->UpwinddF);
        _ArrayAXPY_(solver->m_hyp,-1.0,rhs,size*solver->m_nvars);
      }
    } else if (!strcmp(solver->m_split_hyperbolic_flux,"yes")) {
      if (context->m_flag_hyperbolic_f == _IMPLICIT_) {
        solver->HyperbolicFunction(solver->m_hyp,u,solver,mpi,t,0,solver->FFunction,solver->Upwind);
        _ArrayAXPY_(solver->m_hyp,-1.0,rhs,size*solver->m_nvars);
        solver->HyperbolicFunction(solver->m_hyp,u,solver,mpi,t,0,solver->dFFunction,solver->UpwinddF);
        _ArrayAXPY_(solver->m_hyp, 1.0,rhs,size*solver->m_nvars);
      }
      if (context->m_flag_hyperbolic_df == _IMPLICIT_) {
        solver->HyperbolicFunction(solver->m_hyp,u,solver,mpi,t,0,solver->dFFunction,solver->UpwinddF);
        _ArrayAXPY_(solver->m_hyp,-1.0,rhs,size*solver->m_nvars);
      }
    } else {
      if (context->m_flag_hyperbolic == _IMPLICIT_) {
        solver->HyperbolicFunction(solver->m_hyp,u,solver,mpi,t,0,solver->FFunction,solver->Upwind);
        _ArrayAXPY_(solver->m_hyp,-1.0,rhs,size*solver->m_nvars);
      }
    }
    if (context->m_flag_parabolic == _IMPLICIT_) {
      solver->ParabolicFunction (solver->m_par,u,solver,mpi,t);
      _ArrayAXPY_(solver->m_par, 1.0,rhs,size*solver->m_nvars);
    }
    if (context->m_flag_source == _IMPLICIT_) {
      solver->SourceFunction    (solver->m_source,u,solver,mpi,t);
      _ArrayAXPY_(solver->m_source, 1.0,rhs,size*solver->m_nvars);
    }

    /* save a copy of the solution and RHS for use in IJacobian */
    _ArrayCopy1D_(u  ,solver->m_uref  ,(size*solver->m_nvars));
    _ArrayCopy1D_(rhs,solver->m_rhsref,(size*solver->m_nvars));

    /* Transfer RHS to PETSc vector */
    TransferVecToPETSc(rhs,F,context,ns,context->m_offsets[ns]);

  }

  /* LHS = Ydot - F(u) */
  VecAYPX(F,-1.0,Ydot);

  PetscFunctionReturn(0);
}

#endif

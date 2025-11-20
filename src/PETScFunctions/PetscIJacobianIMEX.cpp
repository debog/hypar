/*! @file PetscIJacobianIMEX.cpp
    @brief Contains the functions required for Jacobian computations for IMEX time-integration
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
#define __FUNCT__ "PetscIJacobianIMEX"
/*!
    Compute the Jacobian for implicit-explicit (IMEX) time integration of the governing equations:
    The ODE, obtained after discretizing the governing PDE in space, is expressed as follows:
    \f{equation}{
      \frac {d{\bf U}}{dt} = {\bf F}\left({\bf U}\right) + {\bf G}\left({\bf U}\right)
      \Rightarrow \dot{\bf U} - {\bf G}\left({\bf U}\right) = {\bf F}\left({\bf U}\right),
    \f}
    where \f${\bf F}\f$ is the spatially discretized non-stiff right-hand-side (treated explicitly),
    \f${\bf G}\f$ is the spatially discretized stiff right-hand-side (treated implicitly), and \f${\bf U}\f$
    represents the entire solution vector. The Jacobian of the implicit part is thus given by:
    \f{equation}{
      {\bf J} = \left[\alpha{\bf I} - \frac {\partial {\bf G}} {\partial {\bf U}} \right]
    \f}
    where \f$\alpha\f$ is the shift coefficient (#PETScContext::shift) of the time integration method.

    Note that \f${\bf G}\left({\bf U}\right)\f$ represents all terms that the user has indicated to be
    integrated in time implicitly (#PETScContext::flag_hyperbolic_f, #PETScContext::flag_hyperbolic_df,
    #PETScContext::flag_hyperbolic, #PETScContext::flag_parabolic, and #PETScContext::flag_source).

    \b Matrix-free \b representation: The Jacobian is computed using a matrix-free approach, where the
    entire Jacobian is never assembled, computed, or stored. Instead, the action of the Jacobian on
    a vector is defined through functions (PetscJacobianFunctionIMEX_JFNK() for nonlinear systems, and
    PetscJacobianFunctionIMEX_Linear() for linear systems). This approach works well with iterative
    solvers. Thus, this function just does the following:
    + Saves the #PETScContext::shift (\f$\alpha\f$) and #PETScContext::waqt (current simulation time)
      to the application context (so that PetscJacobianFunctionIMEX_JFNK() or PetscJacobianFunctionIMEX_Linear()
      can access these values).
    + If a preconditioner is being used, calls the function to compute the preconditioning matrix.

    \sa PetscIFunctionIMEX()

    \b Notes:
    + The Jacobian is defined as the PETSc type MatShell
      (https://petsc.org/release/docs/manualpages/Mat/MATSHELL.html)
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
PetscErrorCode PetscIJacobianIMEX(TS a_ts,        /*!< Time stepping object (see PETSc TS)*/
                                  PetscReal a_t,  /*!< Current time */
                                  Vec a_Y,        /*!< Solution vector */
                                  Vec a_Ydot,     /*!< Time-derivative of solution vector */
                                  PetscReal a_a,  /*!< Shift */
                                  Mat a_A,        /*!< Jacobian matrix */
                                  Mat a_B,        /*!< Preconditioning matrix */
                                  void *a_ctxt    /*!< Application context */ )
{
  PETScContext* context = (PETScContext*) a_ctxt;

  PetscFunctionBegin;
  for (int ns = 0; ns < context->m_nsims; ns++) {
    ((SimulationObject*)context->m_simobj)[ns].solver.m_count_i_jacobian++;
  }
  context->m_shift = a_a;
  context->m_waqt  = a_t;
  /* Construct preconditioning matrix */
  if (context->m_flag_use_precon) PetscComputePreconMatIMEX(a_B,a_Y,context);

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PetscJacobianFunctionIMEX_JFNK"
/*!
    Computes the action of the Jacobian on a vector: See documentation for PetscIJacobianIMEX()
    for the definition of the Jacobian. This function computes its action on a vector
    for a nonlinear system by taking the directional derivative, as follows:
    \f{equation}{
      {\bf f} = {\bf J}\left({\bf U}_0\right){\bf y} = \frac {\partial \mathcal{F}\left({\bf U}\right)} {\partial {\bf U}}\left({\bf U}_0\right) {\bf y} \approx \frac{1}{\epsilon} \left[ \mathcal{F}\left({\bf U}_0+\epsilon{\bf y}\right)-\mathcal{F}\left({\bf U}_0\right) \right]
    \f}
    In the context of the implicit-explicit (IMEX) time integration of the ODE given by
    \f{equation}{
      \frac {d {\bf U}} {dt} = {\bf F}\left({\bf U}\right) + {\bf G}\left({\bf U}\right)
      \Rightarrow \frac {d {\bf U}} {dt} - {\bf G}\left({\bf U}\right) = {\bf F}\left({\bf U}\right),
    \f}
    where \f${\bf F}\f$ is the spatially discretized non-stiff right-hand-side (treated explicitly),
    \f${\bf G}\f$ is the spatially discretized stiff right-hand-side (treated implicitly), we have that
    \f{equation}{
      \mathcal{F}\left({\bf U}\right) \equiv \dot{\bf U} - {\bf G}\left({\bf U}\right)
      \Rightarrow {\bf J}\left({\bf U}\right) = \left[\alpha{\bf I} - \frac {\partial {\bf G}\left({\bf U}\right)} {\partial {\bf U}} \right].
    \f}
    where \f$\alpha\f$ (#PETScContext::shift) is the shift variable specific to the time integration method.
    So this function computes
    \f{equation}{
      {\bf f} = \alpha {\bf y} - \frac{1}{\epsilon} \left[ {\bf G}\left({\bf U}_0+\epsilon{\bf y}\right)-{\bf G}\left({\bf U}_0\right) \right]
    \f}
    In the code, \f${\bf y}\f$ is \a Y, \f$\bf f\f$ is \a F, and \f${\bf U}_0\f$ is \a #HyPar::m_uref (the reference solution at which the
    nonlinear Jacobian is computed). See papers on Jacobian-free Newton-Krylov (JFNK) methods to understand how \f$\epsilon\f$ is computed.

    Note that \f${\bf G}\left({\bf U}\right)\f$ represents all terms that the user has indicated to be
    integrated in time implicitly (#PETScContext::flag_hyperbolic_f, #PETScContext::flag_hyperbolic_df,
    #PETScContext::flag_hyperbolic, #PETScContext::flag_parabolic, and #PETScContext::flag_source).

    \b Notes:
    + For a nonlinear spatial discretization, the right-hand-side \b must be computed without the nonlinearity
      (i.e. with a previously computed or "frozen" discretization operator). This ensures that the Jacobian being
      computed is consistent.
    + It is assumed that the reader is familiar with PETSc's implementation of IMEX time integrators, for
      example, TSARKIMEX (https://petsc.org/release/docs/manualpages/TS/TSARKIMEX.html).
    + See https://petsc.org/release/docs/manualpages/TS/index.html for documentation on
      PETSc's time integrators.
    + All functions and variables whose names start with Vec, Mat, PC, KSP, SNES, and TS are defined by PETSc. Refer to
      the PETSc documentation (https://petsc.org/release/docs/). Usually, googling with the function
      or variable name yields the specific doc page dealing with that function/variable.
*/
PetscErrorCode PetscJacobianFunctionIMEX_JFNK(Mat a_Jacobian, /*!< Jacobian matrix */
                                              Vec a_Y,/*!< Input vector */
                                              Vec a_F /*!< Output vector (Jacobian times input vector */ )
{
  PETScContext* context(nullptr);

  PetscFunctionBegin;

  MatShellGetContext(a_Jacobian,&context);
  SimulationObject* sim = (SimulationObject*) context->m_simobj;
  int nsims = context->m_nsims;

  double normY;
  VecNorm(a_Y,NORM_2,&normY);

  if (normY < 1e-16) {

    /* a_F = 0 */
    VecZeroEntries(a_F);
    /* [J]a_Y = aY - a_F(a_Y) */
    VecAXPBY(a_F,context->m_shift,0,a_Y);

  } else {

    double epsilon =  context->m_jfnk_eps / normY;
    double t = context->m_waqt; /* current stage/step time */

    for (int ns = 0; ns < nsims; ns++) {

      HyPar* solver = &(sim[ns].solver);
      MPIVariables* mpi = &(sim[ns].mpi);
      solver->m_count_i_jac_function++;

      int size = solver->m_npoints_local_wghosts;

      double *u       = solver->m_u;
      double *uref    = solver->m_uref;
      double *rhsref  = solver->m_rhsref;
      double *rhs     = solver->m_rhs;

      /* copy solution from PETSc vector */
      TransferVecFromPETSc(u,a_Y,context,ns,context->m_offsets[ns]);
      _ArrayAYPX_(uref,epsilon,u,size*solver->m_nvars);
      /* apply boundary conditions and exchange data over MPI interfaces */
      solver->ApplyBoundaryConditions(solver,mpi,u,NULL,t);
      MPIExchangeBoundariesnD(  solver->m_ndims,
                                solver->m_nvars,
                                solver->m_dim_local,
                                solver->m_ghosts,
                                mpi,
                                u );

      /* Evaluate hyperbolic, parabolic and source terms  and the RHS for U+dU */
      _ArraySetValue_(rhs,size*solver->m_nvars,0.0);
      if ((!strcmp(solver->m_split_hyperbolic_flux,"yes")) && solver->m_flag_fdf_specified) {
        if (context->m_flag_hyperbolic_f == _IMPLICIT_) {
          solver->HyperbolicFunction( solver->m_hyp,u,solver,mpi,t,0,
                                      solver->FdFFunction,solver->UpwindFdF);
          _ArrayAXPY_(solver->m_hyp,-1.0,rhs,size*solver->m_nvars);
        }
        if (context->m_flag_hyperbolic_df == _IMPLICIT_) {
          solver->HyperbolicFunction( solver->m_hyp,u,solver,mpi,t,0,
                                      solver->dFFunction,solver->UpwinddF);
          _ArrayAXPY_(solver->m_hyp,-1.0,rhs,size*solver->m_nvars);
        }
      } else if (!strcmp(solver->m_split_hyperbolic_flux,"yes")) {
        if (context->m_flag_hyperbolic_f == _IMPLICIT_) {
          solver->HyperbolicFunction( solver->m_hyp,u,solver,mpi,t,0,
                                      solver->FFunction,solver->Upwind);
          _ArrayAXPY_(solver->m_hyp,-1.0,rhs,size*solver->m_nvars);
          solver->HyperbolicFunction( solver->m_hyp,u,solver,mpi,t,0,
                                      solver->dFFunction,solver->UpwinddF);
          _ArrayAXPY_(solver->m_hyp, 1.0,rhs,size*solver->m_nvars);
        }
        if (context->m_flag_hyperbolic_df == _IMPLICIT_) {
          solver->HyperbolicFunction( solver->m_hyp,u,solver,mpi,t,0,
                                      solver->dFFunction,solver->UpwinddF);
          _ArrayAXPY_(solver->m_hyp,-1.0,rhs,size*solver->m_nvars);
        }
      } else {
        if (context->m_flag_hyperbolic == _IMPLICIT_) {
          solver->HyperbolicFunction( solver->m_hyp,u,solver,mpi,t,0,
                                      solver->FFunction,solver->Upwind);
          _ArrayAXPY_(solver->m_hyp,-1.0,rhs,size*solver->m_nvars);
        }
      }
      if (context->m_flag_parabolic == _IMPLICIT_) {
        solver->ParabolicFunction (solver->m_par,u,solver,mpi,t);
        _ArrayAXPY_(solver->m_par, 1.0,rhs,size*solver->m_nvars);
      }
      if (context->m_flag_source == _IMPLICIT_) {
        solver->SourceFunction (solver->m_source,u,solver,mpi,t);
        _ArrayAXPY_(solver->m_source, 1.0,rhs,size*solver->m_nvars);
      }

      _ArrayAXPY_(rhsref,-1.0,rhs,size*solver->m_nvars);
      /* Transfer RHS to PETSc vector */
      TransferVecToPETSc(rhs,a_F,context,ns,context->m_offsets[ns]);
    }

    /* [J]a_Y = aY - a_F(a_Y) */
    VecAXPBY(a_F,context->m_shift,(-1.0/epsilon),a_Y);

  }

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PetscJacobianFunctionIMEX_Linear"
/*!
    Computes the action of the Jacobian on a vector: See documentation for PetscIJacobianIMEX()
    for the definition of the Jacobian. This function computes its action on a vector
    for a linear system by taking the directional derivative, as follows:
    \f{equation}{
      {\bf f} = {\bf J}{\bf y} = \frac {\partial \mathcal{F}\left({\bf U}\right)} {\partial {\bf U}} {\bf y} = \left[ \mathcal{F}\left({\bf U}_0+{\bf y}\right)-\mathcal{F}\left({\bf U}_0\right) \right],
    \f}
    where \f$\mathcal{F}\f$ is linear, and thus \f${\bf J}\f$ is a constant (\f$\mathcal{F}\left({\bf y}\right) = {\bf J}{\bf y}\f$).
    In the context of the implicit-explicit (IMEX) time integration of the ODE given by
    \f{equation}{
      \frac {d {\bf U}} {dt} = {\bf F}\left({\bf U}\right) + {\bf G}\left({\bf U}\right)
      \Rightarrow \frac {d {\bf U}} {dt} - {\bf G}\left({\bf U}\right) = {\bf F}\left({\bf U}\right),
    \f}
    where \f${\bf F}\f$ is the spatially discretized non-stiff right-hand-side (treated explicitly),
    \f${\bf G}\f$ is the spatially discretized stiff right-hand-side (treated implicitly) and \b is \b linear, we have that
    \f{equation}{
      \mathcal{F}\left({\bf U}\right) \equiv \dot{\bf U} - {\bf G}\left({\bf U}\right)
      \Rightarrow {\bf J}\left({\bf U}\right) = \left[\alpha{\bf I} - \frac {\partial {\bf G}\left({\bf U}\right)} {\partial {\bf U}} \right].
    \f}
    where \f$\alpha\f$ (#PETScContext::shift) is the shift variable specific to the time integration method.
    So this function computes
    \f{equation}{
      {\bf f} = \alpha {\bf y} - \left[ {\bf G}\left({\bf U}_0+{\bf y}\right)-{\bf G}\left({\bf U}_0\right) \right]
    \f}
    In the code, \f${\bf y}\f$ is \a Y, \f$\bf f\f$ is \a F, and \f${\bf U}_0\f$ is \a #HyPar::m_uref (the reference solution at which the
    nonlinear Jacobian is computed).

    Since \f$\mathcal{F}\f$ is linear,
    \f{equation}{
      {\bf J}{\bf y} = \left[ \mathcal{F}\left({\bf U}_0+{\bf y}\right)-\mathcal{F}\left({\bf U}_0\right) \right]
                     = \mathcal{F}\left({\bf y}\right).
    \f}
    However, the Jacobian is not computed as \f$\mathcal{F}\left({\bf y}\right)\f$ because of the following reason:
    this function is used by the equation solver within the implicit time integrator in PETSc, and \f${\bf y}\f$
    represents the change in the solution, i.e. \f$\Delta {\bf U}\f$, and not the solution \f$\bf U\f$. Thus,
    evaluating \f$\mathcal{F}\left({\bf y}\right)\f$ using #HyPar::HyperbolicFunction, #HyPar::ParabolicFunction,
    and #HyPar::SourceFunction is incorrect since these functions expect \f$\bf U\f$ as the input.

    Note that \f${\bf G}\left({\bf U}\right)\f$ represents all terms that the user has indicated to be
    integrated in time implicitly (#PETScContext::flag_hyperbolic_f, #PETScContext::flag_hyperbolic_df,
    #PETScContext::flag_hyperbolic, #PETScContext::flag_parabolic, and #PETScContext::flag_source).

    \b Notes:
    + For a nonlinear spatial discretization, the right-hand-side \b must be computed without the nonlinearity
      (i.e. with a previously computed or "frozen" discretization operator). This ensures that the Jacobian being
      computed is consistent, and is truly linear.
    + See https://petsc.org/release/docs/manualpages/TS/index.html for documentation on
      PETSc's time integrators.
    + All functions and variables whose names start with Vec, Mat, PC, KSP, SNES, and TS are defined by PETSc. Refer to
      the PETSc documentation (https://petsc.org/release/docs/). Usually, googling with the function
      or variable name yields the specific doc page dealing with that function/variable.
*/
PetscErrorCode PetscJacobianFunctionIMEX_Linear(Mat a_Jacobian, /*!< Jacobian matrix */
                                                Vec a_Y,/*!<  Input vector */
                                                Vec a_F /*!<  Output vector
                                                            (Jacobian times input vector */)
{
  PETScContext* context(nullptr);

  PetscFunctionBegin;

  MatShellGetContext(a_Jacobian,&context);
  SimulationObject* sim = (SimulationObject*) context->m_simobj;
  int nsims = context->m_nsims;

  double normY;
  VecNorm(a_Y,NORM_2,&normY);

  if (normY < 1e-16) {

    /* a_F = 0 */
    VecZeroEntries(a_F);
    /* [J]a_Y = aY - a_F(a_Y) */
    VecAXPBY(a_F,context->m_shift,0,a_Y);

  } else {

    double t = context->m_waqt; /* current stage/step time */

    for (int ns = 0; ns < nsims; ns++) {

      HyPar* solver = &(sim[ns].solver);
      MPIVariables* mpi = &(sim[ns].mpi);
      solver->m_count_i_jac_function++;

      int size = solver->m_npoints_local_wghosts;

      double *u       = solver->m_u;
      double *uref    = solver->m_uref;
      double *rhsref  = solver->m_rhsref;
      double *rhs     = solver->m_rhs;

      /* copy solution from PETSc vector */
      TransferVecFromPETSc(u,a_Y,context,ns,context->m_offsets[ns]);
      _ArrayAYPX_(uref,1.0,u,size*solver->m_nvars);
      /* apply boundary conditions and exchange data over MPI interfaces */
      solver->ApplyBoundaryConditions(solver,mpi,u,NULL,t);
      MPIExchangeBoundariesnD(  solver->m_ndims,
                                solver->m_nvars,
                                solver->m_dim_local,
                                solver->m_ghosts,
                                mpi,
                                u );

      /* Evaluate hyperbolic, parabolic and source terms  and the RHS for U+dU */
      _ArraySetValue_(rhs,size*solver->m_nvars,0.0);
      if ((!strcmp(solver->m_split_hyperbolic_flux,"yes")) && solver->m_flag_fdf_specified) {
        if (context->m_flag_hyperbolic_f == _IMPLICIT_) {
          solver->HyperbolicFunction( solver->m_hyp,u,solver,mpi,t,0,
                                      solver->FdFFunction,solver->UpwindFdF);
          _ArrayAXPY_(solver->m_hyp,-1.0,rhs,size*solver->m_nvars);
        }
        if (context->m_flag_hyperbolic_df == _IMPLICIT_) {
          solver->HyperbolicFunction( solver->m_hyp,u,solver,mpi,t,0,
                                      solver->dFFunction,solver->UpwinddF);
          _ArrayAXPY_(solver->m_hyp,-1.0,rhs,size*solver->m_nvars);
        }
      } else if (!strcmp(solver->m_split_hyperbolic_flux,"yes")) {
        if (context->m_flag_hyperbolic_f == _IMPLICIT_) {
          solver->HyperbolicFunction( solver->m_hyp,u,solver,mpi,t,0,
                                      solver->FFunction,solver->Upwind);
          _ArrayAXPY_(solver->m_hyp,-1.0,rhs,size*solver->m_nvars);
          solver->HyperbolicFunction( solver->m_hyp,u,solver,mpi,t,0,
                                      solver->dFFunction,solver->UpwinddF);
          _ArrayAXPY_(solver->m_hyp, 1.0,rhs,size*solver->m_nvars);
        }
        if (context->m_flag_hyperbolic_df == _IMPLICIT_) {
          solver->HyperbolicFunction( solver->m_hyp,u,solver,mpi,t,0,
                                      solver->dFFunction,solver->UpwinddF);
          _ArrayAXPY_(solver->m_hyp,-1.0,rhs,size*solver->m_nvars);
        }
      } else {
        if (context->m_flag_hyperbolic == _IMPLICIT_) {
          solver->HyperbolicFunction( solver->m_hyp,u,solver,mpi,t,0,
                                      solver->FFunction,solver->Upwind);
          _ArrayAXPY_(solver->m_hyp,-1.0,rhs,size*solver->m_nvars);
        }
      }
      if (context->m_flag_parabolic == _IMPLICIT_) {
        solver->ParabolicFunction (solver->m_par,u,solver,mpi,t);
        _ArrayAXPY_(solver->m_par, 1.0,rhs,size*solver->m_nvars);
      }
      if (context->m_flag_source == _IMPLICIT_) {
        solver->SourceFunction (solver->m_source,u,solver,mpi,t);
        _ArrayAXPY_(solver->m_source, 1.0,rhs,size*solver->m_nvars);
      }

      _ArrayAXPY_(rhsref,-1.0,rhs,size*solver->m_nvars);
      /* Transfer RHS to PETSc vector */
      TransferVecToPETSc(rhs,a_F,context,ns,context->m_offsets[ns]);
    }

    /* [J]a_Y = aY - a_F(a_Y) */
    VecAXPBY(a_F,context->m_shift,-1.0,a_Y);

  }

  PetscFunctionReturn(0);
}

#endif

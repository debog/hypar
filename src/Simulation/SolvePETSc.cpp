/*! @file SolvePETSc.cpp
    @brief Integrate in time using PETSc
    @author Debojyoti Ghosh

    Integrate the spatially discretized system in time using PETSc's TS module.\n
    (https://petsc.org/release/docs/manualpages/TS/index.html)
*/

#ifdef with_petsc

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <basic.h>
#include <io_cpp.h>
#include <petscinterface.h>
#include <simulation_object.h>

#undef __FUNCT__
#define __FUNCT__ "SolvePETSc"

extern "C" int CalculateError (void*,void*); /*!< Calculate the error in the final solution */
extern "C" int OutputSolution (void*,int);   /*!< Write solutions to file */

/*! \brief Integrate in time with PETSc

    This function integrates the semi-discrete ODE (obtained from discretizing
    the PDE in space) using the time-integration module of PETSc 
    (https://petsc.org/release/docs/manualpages/TS/index.html).
    The time-integration context is set up using the parameters specified in 
    the input file. However, they can also be specified using PETSc's command
    line inputs.\n
    \n
    See PETSc's documentation and examples for more details on how to use its
    TS module. All functions and data types whose names start with Vec, Mat, 
    PC, KSP, SNES, and TS are PETSc functions - refer to the PETSc documentation 
    (usually googling with the function name shows the man page for that function 
    on PETSc's website).
*/
int SolvePETSc( void* s, /*!< Array of simulation objects of type #SimulationObject */
                int   nsims,  /*!< number of simulation objects */
                int   rank,   /*!< MPI rank of this process */
                int   nproc   /*!< Number of MPI processes */ )
{
  SimulationObject* sim = (SimulationObject*) s;

  TS              ts;     /* time integration object               */
  Vec             Y,Z;    /* PETSc solution vectors                */
  Mat             A, B;   /* Jacobian and preconditioning matrices */
  TSType          time_scheme;  /* time integration method         */
  TSProblemType   ptype;  /* problem type - nonlinear or linear    */
  int             flag_mat_a = 0, flag_mat_b = 0, iAuxSize = 0, i;

  PetscFunctionBegin;

  /* Register custom time-integration methods, if specified */
  PetscRegisterTIMethods(rank);
  if(!rank) printf("Setting up PETSc time integration... \n");

  /* create and set a PETSc context */
  PETScContext context;

  context.rank = rank;
  context.nproc = nproc;

  context.simobj = sim;
  context.nsims = nsims;

  /* default: everything explicit */
  context.flag_hyperbolic     = _EXPLICIT_; 
  context.flag_hyperbolic_f   = _EXPLICIT_; 
  context.flag_hyperbolic_df  = _EXPLICIT_; 
  context.flag_parabolic      = _EXPLICIT_; 
  context.flag_source         = _EXPLICIT_; 

  context.tic = 0;
  context.flag_is_linear = 0;
  context.globalDOF.clear();
  context.points.clear();

  PetscCreatePointList(&context);

  /* create and initialize PETSc solution vector and other parameters */
  /* PETSc solution vector does not have ghost points */
  VecCreate(MPI_COMM_WORLD,&Y);
  VecSetSizes(Y,context.ndofs,PETSC_DECIDE);
  VecSetUp(Y);

  /* copy initial solution to PETSc's vector */
  for (int ns = 0; ns < nsims; ns++) {
    TransferVecToPETSc( sim[ns].solver.u,
                        Y,
                        &context,
                        ns,
                        context.offsets[ns] );
  }

  /* Create the global DOF mapping for all the grid points */
  PetscGlobalDOF(&context);

  /* Define and initialize the time-integration object */
  TSCreate(MPI_COMM_WORLD,&ts);
  TSSetMaxSteps(ts,sim[0].solver.n_iter);
  TSSetMaxTime(ts,sim[0].solver.dt*sim[0].solver.n_iter);
  TSSetTimeStep(ts,sim[0].solver.dt);
  TSSetTime(ts,0.0);
  TSSetExactFinalTime(ts,TS_EXACTFINALTIME_MATCHSTEP);
  TSSetType(ts,TSBEULER);

  /* set default time step adaptivity to none */
  TSAdapt adapt;
  TSGetAdapt(ts,&adapt);
  TSAdaptSetType(adapt,TSADAPTNONE);

  /* set options from input */
  TSSetFromOptions(ts);

  /* Define the right and left -hand side functions for each time-integration scheme */
  TSGetType(ts,&time_scheme);
  TSGetProblemType(ts,&ptype);
  if (!strcmp(time_scheme,TSARKIMEX)) {

    /* implicit - explicit time integration */

    TSSetRHSFunction(ts,PETSC_NULL,PetscRHSFunctionIMEX,&context);
    TSSetIFunction  (ts,PETSC_NULL,PetscIFunctionIMEX,  &context);

    SNES     snes;
    KSP      ksp;
    PC       pc;
    SNESType snestype;
    TSGetSNES(ts,&snes);
    SNESGetType(snes,&snestype);

    /* Matrix-free representation of the Jacobian */
    flag_mat_a = 1;
    MatCreateShell(MPI_COMM_WORLD,context.ndofs,context.ndofs,PETSC_DETERMINE,
                   PETSC_DETERMINE,&context,&A);
    if ((!strcmp(snestype,SNESKSPONLY)) || (ptype == TS_LINEAR)) {
      /* linear problem */
      context.flag_is_linear = 1;
      MatShellSetOperation(A,MATOP_MULT,(void (*)(void))PetscJacobianFunctionIMEX_Linear);
      SNESSetType(snes,SNESKSPONLY);
    } else {
      /* nonlinear problem */
      context.flag_is_linear = 0;
      context.jfnk_eps = 1e-7;
      PetscOptionsGetReal(NULL,NULL,"-jfnk_epsilon",&context.jfnk_eps,NULL);
      MatShellSetOperation(A,MATOP_MULT,(void (*)(void))PetscJacobianFunctionIMEX_JFNK);
    }
    MatSetUp(A);

    context.flag_use_precon = 0;
    PetscOptionsGetBool(  PETSC_NULL,PETSC_NULL,
                          "-with_pc",
                          (PetscBool*)(&context.flag_use_precon),
                          PETSC_NULL );

    if (context.flag_use_precon) {
      /* check if flux Jacobian of the physical model is defined */
      for (int ns = 0; ns < nsims; ns++) {
        if (!sim[ns].solver.JFunction) {
          if (!rank) {
            fprintf(stderr,"Error in SolvePETSc(): solver->JFunction (point-wise flux Jacobian) must ");
            fprintf(stderr,"be defined for preconditioning.\n");
          }
          PetscFunctionReturn(1);
        }
      }
      /* Set up preconditioner matrix */
      flag_mat_b = 1;
      MatCreateAIJ( MPI_COMM_WORLD,
                    context.ndofs, context.ndofs,
                    PETSC_DETERMINE, PETSC_DETERMINE,
                    (sim[0].solver.ndims*2+1)*sim[0].solver.nvars, NULL,
                    2*sim[0].solver.ndims*sim[0].solver.nvars, NULL,
                    &B );
      MatSetBlockSize(B,sim[0].solver.nvars);
      /* Set the IJacobian function for TS */
      TSSetIJacobian(ts,A,B,PetscIJacobianIMEX,&context);
    } else {
      /* Set the IJacobian function for TS */
      TSSetIJacobian(ts,A,A,PetscIJacobianIMEX,&context);
      /* Set PC (preconditioner) to none */
      SNESGetKSP(snes,&ksp);
      KSPGetPC(ksp,&pc);
      PCSetType(pc,PCNONE);
    }

    /* read the implicit/explicit flags for each of the terms for IMEX schemes */
    /* default -> hyperbolic - explicit, parabolic and source - implicit       */
    PetscBool flag = PETSC_FALSE;

    context.flag_hyperbolic     = _EXPLICIT_; 
    context.flag_hyperbolic_f   = _EXPLICIT_; 
    context.flag_hyperbolic_df  = _IMPLICIT_; 
    context.flag_parabolic      = _IMPLICIT_; 
    context.flag_source         = _IMPLICIT_; 

    if (!strcmp(sim[0].solver.SplitHyperbolicFlux,"yes")) {

      flag = PETSC_FALSE; 
      PetscOptionsGetBool(PETSC_NULL,PETSC_NULL,"-hyperbolic_f_explicit",&flag,PETSC_NULL);
      if (flag == PETSC_TRUE) context.flag_hyperbolic_f = _EXPLICIT_; 
      flag = PETSC_FALSE; 
      PetscOptionsGetBool(PETSC_NULL,PETSC_NULL,"-hyperbolic_f_implicit",&flag,PETSC_NULL);
      if (flag == PETSC_TRUE) context.flag_hyperbolic_f = _IMPLICIT_; 

      flag = PETSC_FALSE; 
      PetscOptionsGetBool(PETSC_NULL,PETSC_NULL,"-hyperbolic_df_explicit",&flag,PETSC_NULL);
      if (flag == PETSC_TRUE) context.flag_hyperbolic_df = _EXPLICIT_; 
      flag = PETSC_FALSE; 
      PetscOptionsGetBool(PETSC_NULL,PETSC_NULL,"-hyperbolic_df_implicit",&flag,PETSC_NULL);
      if (flag == PETSC_TRUE) context.flag_hyperbolic_df = _IMPLICIT_; 

    } else {

      flag = PETSC_FALSE; 
      PetscOptionsGetBool(PETSC_NULL,PETSC_NULL,"-hyperbolic_explicit",&flag,PETSC_NULL);
      if (flag == PETSC_TRUE) context.flag_hyperbolic = _EXPLICIT_; 
      flag = PETSC_FALSE; 
      PetscOptionsGetBool(PETSC_NULL,PETSC_NULL,"-hyperbolic_implicit",&flag,PETSC_NULL);
      if (flag == PETSC_TRUE) context.flag_hyperbolic = _IMPLICIT_; 

    }

    flag = PETSC_FALSE; 
    PetscOptionsGetBool(PETSC_NULL,PETSC_NULL,"-parabolic_explicit",&flag,PETSC_NULL);
    if (flag == PETSC_TRUE) context.flag_parabolic = _EXPLICIT_; 
    flag = PETSC_FALSE; 
    PetscOptionsGetBool(PETSC_NULL,PETSC_NULL,"-parabolic_implicit",&flag,PETSC_NULL);
    if (flag == PETSC_TRUE) context.flag_parabolic = _IMPLICIT_; 

    flag = PETSC_FALSE; 
    PetscOptionsGetBool(PETSC_NULL,PETSC_NULL,"-source_explicit",&flag,PETSC_NULL);
    if (flag == PETSC_TRUE) context.flag_source = _EXPLICIT_; 
    flag = PETSC_FALSE; 
    PetscOptionsGetBool(PETSC_NULL,PETSC_NULL,"-source_implicit",&flag,PETSC_NULL);
    if (flag == PETSC_TRUE) context.flag_source = _IMPLICIT_; 

    flag = PETSC_FALSE;
    PetscOptionsGetBool(PETSC_NULL,PETSC_NULL,"-ts_arkimex_fully_implicit",&flag,PETSC_NULL);
    if (flag == PETSC_TRUE) {
      context.flag_hyperbolic_f   = _IMPLICIT_;
      context.flag_hyperbolic_df  = _IMPLICIT_;
      context.flag_hyperbolic     = _IMPLICIT_;
      context.flag_parabolic      = _IMPLICIT_;
      context.flag_source         = _IMPLICIT_;
    }

    /* print out a summary of the treatment of each term */
    if (!rank) {
      printf("Implicit-Explicit time-integration:-\n");
      if (!strcmp(sim[0].solver.SplitHyperbolicFlux,"yes")) {
        if (context.flag_hyperbolic_f == _EXPLICIT_)  printf("Hyperbolic (f-df) term: Explicit\n");
        else                                          printf("Hyperbolic (f-df) term: Implicit\n");
        if (context.flag_hyperbolic_df == _EXPLICIT_) printf("Hyperbolic (df)   term: Explicit\n");
        else                                          printf("Hyperbolic (df)   term: Implicit\n");
      } else {
        if (context.flag_hyperbolic == _EXPLICIT_)    printf("Hyperbolic        term: Explicit\n");
        else                                          printf("Hyperbolic        term: Implicit\n");
      }
      if (context.flag_parabolic == _EXPLICIT_)       printf("Parabolic         term: Explicit\n");
      else                                            printf("Parabolic         term: Implicit\n");
      if (context.flag_source    == _EXPLICIT_)       printf("Source            term: Explicit\n");
      else                                            printf("Source            term: Implicit\n");
    }

  } else if (     (!strcmp(time_scheme,TSEULER)) 
              ||  (!strcmp(time_scheme,TSRK   )) 
              ||  (!strcmp(time_scheme,TSSSP  )) ) {

    /* Explicit time integration */    
    TSSetRHSFunction(ts,PETSC_NULL,PetscRHSFunctionExpl,&context);

  } else {

    /* Implicit or explicit time integration */

    TSSetRHSFunction(ts,PETSC_NULL,PetscRHSFunctionImpl,&context);

    SNES     snes;
    KSP      ksp;
    PC       pc;
    SNESType snestype;
    TSGetSNES(ts,&snes);
    SNESGetType(snes,&snestype);

    /* Matrix-free representation of the Jacobian */
    flag_mat_a = 1;
    MatCreateShell(MPI_COMM_WORLD,context.ndofs,context.ndofs,PETSC_DETERMINE,
                          PETSC_DETERMINE,&context,&A);
    if ((!strcmp(snestype,SNESKSPONLY)) || (ptype == TS_LINEAR)) {
      /* linear problem */
      context.flag_is_linear = 1;
      MatShellSetOperation(A,MATOP_MULT,(void (*)(void))PetscJacobianFunction_Linear);
      SNESSetType(snes,SNESKSPONLY);
    } else {
      /* nonlinear problem */
      context.flag_is_linear = 0;
      context.jfnk_eps = 1e-7;
      PetscOptionsGetReal(NULL,NULL,"-jfnk_epsilon",&context.jfnk_eps,NULL);
      MatShellSetOperation(A,MATOP_MULT,(void (*)(void))PetscJacobianFunction_JFNK);
    }
    MatSetUp(A);

    context.flag_use_precon = 0;
    PetscOptionsGetBool(  PETSC_NULL,
                          PETSC_NULL,
                          "-with_pc",
                          (PetscBool*)(&context.flag_use_precon),
                          PETSC_NULL );

    /* Since we are using a MatShell to represent the action of the Jacobian 
       on a vector (Jacobian-free approach),
       we need to define the Jacobian through TSSetIJacobian, and not TSSetRHSJacobian, 
       so that we can define the complete operator (the shifted Jacobian aI - J ) */

    if (context.flag_use_precon) {
      /* check if flux Jacobian of the physical model is defined */
      for (int ns = 0; ns < nsims; ns++) {
        if (!sim[ns].solver.JFunction) {
          if (!rank) {
            fprintf(stderr,"Error in SolvePETSc(): solver->JFunction (point-wise flux Jacobian) must ");
            fprintf(stderr,"be defined for preconditioning.\n");
          }
          PetscFunctionReturn(1);
        }
      }
      /* Set up preconditioner matrix */
      flag_mat_b = 1;
      MatCreateAIJ( MPI_COMM_WORLD,
                    context.ndofs, context.ndofs,
                    PETSC_DETERMINE, PETSC_DETERMINE,
                    (sim[0].solver.ndims*2+1)*sim[0].solver.nvars, NULL,
                    2*sim[0].solver.ndims*sim[0].solver.nvars, NULL,
                    &B );
      MatSetBlockSize(B,sim[0].solver.nvars);
      /* Set the RHSJacobian function for TS */
      TSSetIJacobian(ts,A,B,PetscIJacobian,&context);
    } else {
      /* Set the RHSJacobian function for TS */
      TSSetIJacobian(ts,A,A,PetscIJacobian,&context);
      /* Set PC (preconditioner) to none */
      SNESGetKSP(snes,&ksp);
      KSPGetPC(ksp,&pc);
      PCSetType(pc,PCNONE);
    }

  }

  /* Set pre/post-stage and post-timestep function */
  TSSetPreStep (ts,PetscPreTimeStep );
  TSSetPreStage(ts,PetscPreStage    );
  TSSetPostStage(ts,PetscPostStage  );
  TSSetPostStep(ts,PetscPostTimeStep);
  /* Set solution vector for TS */
  TSSetSolution(ts,Y);
  /* Set it all up */
  TSSetUp(ts);
  /* Set application context */
  TSSetApplicationContext(ts,&context);

  if (!rank) {
    if (context.flag_is_linear) printf("SolvePETSc(): Problem type is linear.\n");
    else                        printf("SolvePETSc(): Problem type is nonlinear.\n");
  }

  if (!rank) printf("** Starting PETSc time integration **\n");
  TSSolve(ts,Y);
  if (!rank) printf("** Completed PETSc time integration **\n");

  /* Get the number of time steps */
  for (int ns = 0; ns < nsims; ns++) {
    TSGetStepNumber(ts,&(sim[ns].solver.n_iter));
  }

  /* get and write to file any auxiliary solutions */
  char aux_fname_root[4] = "ts0";
  TSGetSolutionComponents(ts,&iAuxSize,NULL);
  if (iAuxSize) {
    if (iAuxSize > 10) iAuxSize = 10;
    if (!rank) printf("Number of auxiliary solutions from time integration: %d\n",iAuxSize);
    VecDuplicate(Y,&Z);
    for (i=0; i<iAuxSize; i++) {
      TSGetSolutionComponents(ts,&i,&Z);
      for (int ns = 0; ns < nsims; ns++) {
        TransferVecFromPETSc(sim[ns].solver.u,Z,&context,ns,context.offsets[ns]);
        WriteArray( sim[ns].solver.ndims,
                    sim[ns].solver.nvars,
                    sim[ns].solver.dim_global,
                    sim[ns].solver.dim_local,
                    sim[ns].solver.ghosts,
                    sim[ns].solver.x,
                    sim[ns].solver.u,
                    &(sim[ns].solver),
                    &(sim[ns].mpi),
                    aux_fname_root );
      }
      aux_fname_root[2]++;
    }
    VecDestroy(&Z);
  }

  /* if available, get error estimates */
  PetscTimeError(ts);

  /* copy final solution from PETSc's vector */
  for (int ns = 0; ns < nsims; ns++) {
    TransferVecFromPETSc(sim[ns].solver.u,Y,&context,ns,context.offsets[ns]);
  }

  /* clean up */
  if (flag_mat_a) { MatDestroy(&A); }
  if (flag_mat_b) { MatDestroy(&B); }
  TSDestroy(&ts);
  VecDestroy(&Y);

  /* write a final solution file, if last iteration did not write one */
  if (context.tic) { 
    for (int ns = 0; ns < nsims; ns++) {
      HyPar* solver = &(sim[ns].solver);
      MPIVariables* mpi = &(sim[ns].mpi);
      if (solver->PhysicsOutput) {
        solver->PhysicsOutput(solver,mpi);
      }
      CalculateError(solver,mpi);
    }
    OutputSolution(sim, nsims); 
  }
  /* calculate error if exact solution has been provided */
  for (int ns = 0; ns < nsims; ns++) {
    CalculateError(&(sim[ns].solver), &(sim[ns].mpi));
  }

  PetscCleanup(&context);
  PetscFunctionReturn(0);
}

#endif

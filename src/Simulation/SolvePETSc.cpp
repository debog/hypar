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

#ifdef with_librom
#include <librom_interface.h>
#endif

#undef __FUNCT__
#define __FUNCT__ "SolvePETSc"

extern "C" int CalculateError (void*,void*); /*!< Calculate the error in the final solution */
int OutputSolution (void*,int,double);   /*!< Write solutions to file */
extern "C" void ResetFilenameIndex(char*, int); /*!< Reset filename index */
#ifdef with_librom
extern "C" int CalculateROMDiff(void*,void*); /*!< Calculate the diff of PDE and ROM solutions */
int OutputROMSolution(void*,int,double);   /*!< Write ROM solutions to file */
#endif

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

  DM              dm; /* data management object */
  TS              ts; /* time integration object */
  Vec             Y,Z; /* PETSc solution vectors */
  Mat             A, B; /* Jacobian and preconditioning matrices */
  MatFDColoring   fdcoloring; /* coloring for sparse Jacobian computation */
  TSType          time_scheme; /* time integration method */
  TSProblemType   ptype; /* problem type - nonlinear or linear */

  int flag_mat_a = 0,
      flag_mat_b = 0,
      flag_fdcoloring = 0,
      iAuxSize = 0, i;

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
  context.ti_runtime = 0.0;
  context.waqt = 0.0;
  context.dt = sim[0].solver.dt;
  context.stage_times.clear();
  context.stage_index = 0;

#ifdef with_librom
  if (!rank) printf("Setting up libROM interface.\n");
  context.rom_interface = new libROMInterface( sim,
                                               nsims,
                                               rank,
                                               nproc,
                                               sim[0].solver.dt );
  context.rom_mode = ((libROMInterface*)context.rom_interface)->mode();
  context.op_times_arr.clear();
#endif

#ifdef with_librom
  if (      (context.rom_mode == _ROM_MODE_TRAIN_)
        ||  (context.rom_mode == _ROM_MODE_INITIAL_GUESS_ )
        ||  (context.rom_mode == _ROM_MODE_NONE_ ) ) {

    if (context.rom_mode == _ROM_MODE_INITIAL_GUESS_) {
      ((libROMInterface*)context.rom_interface)->loadROM();
      ((libROMInterface*)context.rom_interface)->projectInitialSolution(sim);
    }
#endif
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
    TSSetTime(ts,context.waqt);
    TSSetExactFinalTime(ts,TS_EXACTFINALTIME_MATCHSTEP);
    TSSetType(ts,TSBEULER);

    /* set default time step adaptivity to none */
    TSAdapt adapt;
    TSAdaptType adapt_type = TSADAPTNONE;
    TSGetAdapt(ts,&adapt);
    TSAdaptSetType(adapt,adapt_type);

    /* set options from input */
    TSSetFromOptions(ts);

    /* create DM */
    DMShellCreate(MPI_COMM_WORLD, &dm);
    DMShellSetGlobalVector(dm, Y);
    TSSetDM(ts, dm);

#ifdef with_librom
    TSAdaptGetType(adapt,&adapt_type);
    if (strcmp(adapt_type, TSADAPTNONE)) {
      if (!rank) printf("Warning: libROM interface not yet implemented for adaptive timestepping.\n");
    }
#endif

    /* Define the right and left -hand side functions for each time-integration scheme */
    TSGetType(ts,&time_scheme);
    TSGetProblemType(ts,&ptype);

    if (!strcmp(time_scheme,TSARKIMEX)) {

      /* implicit - explicit time integration */

      TSSetRHSFunction(ts,nullptr,PetscRHSFunctionIMEX,&context);
      TSSetIFunction  (ts,nullptr,PetscIFunctionIMEX,  &context);

      SNES     snes;
      KSP      ksp;
      PC       pc;
      SNESType snestype;
      TSGetSNES(ts,&snes);
      SNESGetType(snes,&snestype);

#ifdef with_librom
      if (context.rom_mode == _ROM_MODE_INITIAL_GUESS_) {
        SNESSetComputeInitialGuess(snes, PetscSetInitialGuessROM, &context);
      }
#endif

      context.flag_use_precon = 0;
      PetscOptionsGetBool(  nullptr,nullptr,
                            "-with_pc",
                            (PetscBool*)(&context.flag_use_precon),
                            nullptr );

      char precon_mat_type_c_st[_MAX_STRING_SIZE_] = "default";
      PetscOptionsGetString(  nullptr,
                              nullptr,
                              "-pc_matrix_type",
                              precon_mat_type_c_st,
                              _MAX_STRING_SIZE_,
                              nullptr );
      context.precon_matrix_type = std::string(precon_mat_type_c_st);

      if (context.flag_use_precon) {

        if (context.precon_matrix_type == "default") {

          /* Matrix-free representation of the Jacobian */
          flag_mat_a = 1;
          MatCreateShell( MPI_COMM_WORLD,
                          context.ndofs,
                          context.ndofs,
                          PETSC_DETERMINE,
                          PETSC_DETERMINE,
                          &context,
                          &A);
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
          /* check if Jacobian of the physical model is defined */
          for (int ns = 0; ns < nsims; ns++) {
            if ((!sim[ns].solver.JFunction) && (!sim[ns].solver.KFunction)) {
              if (!rank) {
                fprintf(stderr,"Error in SolvePETSc(): solver->JFunction  or solver->KFunction ");
                fprintf(stderr,"(point-wise Jacobians for hyperbolic or parabolic terms) must ");
                fprintf(stderr,"be defined for preconditioning.\n");
              }
              PetscFunctionReturn(1);
            }
          }
          /* Set up preconditioner matrix */
          flag_mat_b = 1;
          MatCreateAIJ( MPI_COMM_WORLD,
                        context.ndofs,
                        context.ndofs,
                        PETSC_DETERMINE,
                        PETSC_DETERMINE,
                        (sim[0].solver.ndims*2+1)*sim[0].solver.nvars, NULL,
                        2*sim[0].solver.ndims*sim[0].solver.nvars, NULL,
                        &B );
          MatSetBlockSize(B,sim[0].solver.nvars);
          /* Set the IJacobian function for TS */
          TSSetIJacobian(ts,A,B,PetscIJacobianIMEX,&context);

        } else if (context.precon_matrix_type == "fd") {

          flag_mat_a = 1;
          MatCreateSNESMF(snes,&A);
          flag_mat_b = 1;
          MatCreateAIJ( MPI_COMM_WORLD,
                        context.ndofs,
                        context.ndofs,
                        PETSC_DETERMINE,
                        PETSC_DETERMINE,
                        (sim[0].solver.ndims*2+1)*sim[0].solver.nvars, NULL,
                        2*sim[0].solver.ndims*sim[0].solver.nvars, NULL,
                        &B);
          MatSetOption(B, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
          /* Set the Jacobian function for SNES */
          SNESSetJacobian(snes, A, B, SNESComputeJacobianDefault, NULL);

        } else if (context.precon_matrix_type == "colored_fd") {

          int stencil_width = 1;
          PetscOptionsGetInt( NULL,
                              NULL,
                              "-pc_matrix_colored_fd_stencil_width",
                              &stencil_width,
                              NULL );

          flag_mat_a = 1;
          MatCreateSNESMF(snes,&A);
          flag_mat_b = 1;
          MatCreateAIJ( MPI_COMM_WORLD,
                        context.ndofs,
                        context.ndofs,
                        PETSC_DETERMINE,
                        PETSC_DETERMINE,
                        (sim[0].solver.ndims*2+1)*sim[0].solver.nvars, NULL,
                        2*sim[0].solver.ndims*sim[0].solver.nvars, NULL,
                        &B);
          MatSetOption(B, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
          if (!rank) {
            printf("PETSc:    Setting Jacobian non-zero pattern (stencil width %d).\n",
                    stencil_width );
          }
          PetscJacobianMatNonzeroEntriesImpl(B, stencil_width, &context);

          /* Set the Jacobian function for SNES */
          SNESSetJacobian(snes, A, B, SNESComputeJacobianDefaultColor, NULL);

        } else {

          if (!rank) {
            fprintf(  stderr,"Invalid input for \"-pc_matrix_type\": %s.\n",
                      context.precon_matrix_type.c_str());
          }
          PetscFunctionReturn(0);

        }

        /* set PC side to right */
        SNESGetKSP(snes,&ksp);
        KSPSetPCSide(ksp, PC_RIGHT);

      } else {

        /* Matrix-free representation of the Jacobian */
        flag_mat_a = 1;
        MatCreateShell( MPI_COMM_WORLD,
                        context.ndofs,
                        context.ndofs,
                        PETSC_DETERMINE,
                        PETSC_DETERMINE,
                        &context,
                        &A);
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
        /* Set the RHSJacobian function for TS */
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
        PetscOptionsGetBool(nullptr,nullptr,"-hyperbolic_f_explicit",&flag,nullptr);
        if (flag == PETSC_TRUE) context.flag_hyperbolic_f = _EXPLICIT_;
        flag = PETSC_FALSE;
        PetscOptionsGetBool(nullptr,nullptr,"-hyperbolic_f_implicit",&flag,nullptr);
        if (flag == PETSC_TRUE) context.flag_hyperbolic_f = _IMPLICIT_;

        flag = PETSC_FALSE;
        PetscOptionsGetBool(nullptr,nullptr,"-hyperbolic_df_explicit",&flag,nullptr);
        if (flag == PETSC_TRUE) context.flag_hyperbolic_df = _EXPLICIT_;
        flag = PETSC_FALSE;
        PetscOptionsGetBool(nullptr,nullptr,"-hyperbolic_df_implicit",&flag,nullptr);
        if (flag == PETSC_TRUE) context.flag_hyperbolic_df = _IMPLICIT_;

      } else {

        flag = PETSC_FALSE;
        PetscOptionsGetBool(nullptr,nullptr,"-hyperbolic_explicit",&flag,nullptr);
        if (flag == PETSC_TRUE) context.flag_hyperbolic = _EXPLICIT_;
        flag = PETSC_FALSE;
        PetscOptionsGetBool(nullptr,nullptr,"-hyperbolic_implicit",&flag,nullptr);
        if (flag == PETSC_TRUE) context.flag_hyperbolic = _IMPLICIT_;

      }

      flag = PETSC_FALSE;
      PetscOptionsGetBool(nullptr,nullptr,"-parabolic_explicit",&flag,nullptr);
      if (flag == PETSC_TRUE) context.flag_parabolic = _EXPLICIT_;
      flag = PETSC_FALSE;
      PetscOptionsGetBool(nullptr,nullptr,"-parabolic_implicit",&flag,nullptr);
      if (flag == PETSC_TRUE) context.flag_parabolic = _IMPLICIT_;

      flag = PETSC_FALSE;
      PetscOptionsGetBool(nullptr,nullptr,"-source_explicit",&flag,nullptr);
      if (flag == PETSC_TRUE) context.flag_source = _EXPLICIT_;
      flag = PETSC_FALSE;
      PetscOptionsGetBool(nullptr,nullptr,"-source_implicit",&flag,nullptr);
      if (flag == PETSC_TRUE) context.flag_source = _IMPLICIT_;

      flag = PETSC_FALSE;
      PetscOptionsGetBool(nullptr,nullptr,"-ts_arkimex_fully_implicit",&flag,nullptr);
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
      TSSetRHSFunction(ts,nullptr,PetscRHSFunctionExpl,&context);

    } else if (     (!strcmp(time_scheme,TSCN))
                ||  (!strcmp(time_scheme,TSBEULER )) ) {


      /* Implicit time integration */

      TSSetIFunction(ts,nullptr,PetscIFunctionImpl,&context);

      SNES     snes;
      KSP      ksp;
      PC       pc;
      SNESType snestype;
      TSGetSNES(ts,&snes);
      SNESGetType(snes,&snestype);

#ifdef with_librom
      if (context.rom_mode == _ROM_MODE_INITIAL_GUESS_) {
        SNESSetComputeInitialGuess(snes, PetscSetInitialGuessROM, &context);
      }
#endif

      context.flag_use_precon = 0;
      PetscOptionsGetBool(  nullptr,
                            nullptr,
                            "-with_pc",
                            (PetscBool*)(&context.flag_use_precon),
                            nullptr );

      char precon_mat_type_c_st[_MAX_STRING_SIZE_] = "default";
      PetscOptionsGetString(  nullptr,
                              nullptr,
                              "-pc_matrix_type",
                              precon_mat_type_c_st,
                              _MAX_STRING_SIZE_,
                              nullptr );
      context.precon_matrix_type = std::string(precon_mat_type_c_st);

      if (context.flag_use_precon) {

        if (context.precon_matrix_type == "default") {

          /* Matrix-free representation of the Jacobian */
          flag_mat_a = 1;
          MatCreateShell( MPI_COMM_WORLD,
                          context.ndofs,
                          context.ndofs,
                          PETSC_DETERMINE,
                          PETSC_DETERMINE,
                          &context,
                          &A);
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
          /* check if Jacobian of the physical model is defined */
          for (int ns = 0; ns < nsims; ns++) {
            if ((!sim[ns].solver.JFunction) && (!sim[ns].solver.KFunction)) {
              if (!rank) {
                fprintf(stderr,"Error in SolvePETSc(): solver->JFunction  or solver->KFunction ");
                fprintf(stderr,"(point-wise Jacobians for hyperbolic or parabolic terms) must ");
                fprintf(stderr,"be defined for preconditioning.\n");
              }
              PetscFunctionReturn(1);
            }
          }
          /* Set up preconditioner matrix */
          flag_mat_b = 1;
          MatCreateAIJ( MPI_COMM_WORLD,
                        context.ndofs,
                        context.ndofs,
                        PETSC_DETERMINE,
                        PETSC_DETERMINE,
                        (sim[0].solver.ndims*2+1)*sim[0].solver.nvars, NULL,
                        2*sim[0].solver.ndims*sim[0].solver.nvars, NULL,
                        &B );
          MatSetBlockSize(B,sim[0].solver.nvars);
          /* Set the IJacobian function for TS */
          TSSetIJacobian(ts,A,B,PetscIJacobian,&context);

        } else if (context.precon_matrix_type == "fd") {

          flag_mat_a = 1;
          MatCreateSNESMF(snes,&A);
          flag_mat_b = 1;
          MatCreateAIJ( MPI_COMM_WORLD,
                        context.ndofs,
                        context.ndofs,
                        PETSC_DETERMINE,
                        PETSC_DETERMINE,
                        (sim[0].solver.ndims*2+1)*sim[0].solver.nvars, NULL,
                        2*sim[0].solver.ndims*sim[0].solver.nvars, NULL,
                        &B);
          MatSetOption(B, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
          /* Set the Jacobian function for SNES */
          SNESSetJacobian(snes, A, B, SNESComputeJacobianDefault, NULL);

        } else if (context.precon_matrix_type == "colored_fd") {

          int stencil_width = 1;
          PetscOptionsGetInt( NULL,
                              NULL,
                              "-pc_matrix_colored_fd_stencil_width",
                              &stencil_width,
                              NULL );

          flag_mat_a = 1;
          MatCreateSNESMF(snes,&A);
          flag_mat_b = 1;
          MatCreateAIJ( MPI_COMM_WORLD,
                        context.ndofs,
                        context.ndofs,
                        PETSC_DETERMINE,
                        PETSC_DETERMINE,
                        (sim[0].solver.ndims*2+1)*sim[0].solver.nvars, NULL,
                        2*sim[0].solver.ndims*sim[0].solver.nvars, NULL,
                        &B);
          MatSetOption(B, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
          if (!rank) {
            printf("PETSc:    Setting Jacobian non-zero pattern (stencil width %d).\n",
                    stencil_width );
          }
          PetscJacobianMatNonzeroEntriesImpl(B, stencil_width, &context);

          /* Set the Jacobian function for SNES */
          SNESSetJacobian(snes, A, B, SNESComputeJacobianDefaultColor, NULL);

        } else {

          if (!rank) {
            fprintf(  stderr,"Invalid input for \"-pc_matrix_type\": %s.\n",
                      context.precon_matrix_type.c_str());
          }
          PetscFunctionReturn(0);

        }

        /* set PC side to right */
        SNESGetKSP(snes,&ksp);
        KSPSetPCSide(ksp, PC_RIGHT);

      } else {

        /* Matrix-free representation of the Jacobian */
        flag_mat_a = 1;
        MatCreateShell( MPI_COMM_WORLD,
                        context.ndofs,
                        context.ndofs,
                        PETSC_DETERMINE,
                        PETSC_DETERMINE,
                        &context,
                        &A);
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
        /* Set the RHSJacobian function for TS */
        TSSetIJacobian(ts,A,A,PetscIJacobian,&context);
        /* Set PC (preconditioner) to none */
        SNESGetKSP(snes,&ksp);
        KSPGetPC(ksp,&pc);
        PCSetType(pc,PCNONE);
      }

    } else {

      if (!rank) {
        fprintf(stderr, "Time integration type %s is not yet supported.\n", time_scheme);
      }
      PetscFunctionReturn(0);

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
    context.ti_runtime = 0.0;
    TSSolve(ts,Y);
    if (!rank) {
      printf("** Completed PETSc time integration (Final time: %f), total wctime: %f (seconds) **\n",
              context.waqt, context.ti_runtime );
    }

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
    VecDestroy(&Y);
    if (flag_mat_a) { MatDestroy(&A); }
    if (flag_mat_b) { MatDestroy(&B); }
    if (flag_fdcoloring) { MatFDColoringDestroy(&fdcoloring); }
    TSDestroy(&ts);
    DMDestroy(&dm);

    /* write a final solution file, if last iteration did not write one */
    if (context.tic) {
      for (int ns = 0; ns < nsims; ns++) {
        HyPar* solver = &(sim[ns].solver);
        MPIVariables* mpi = &(sim[ns].mpi);
        if (solver->PhysicsOutput) {
          solver->PhysicsOutput(solver,mpi, context.waqt);
        }
        CalculateError(solver,mpi);
      }
      OutputSolution(sim, nsims, context.waqt);
    }
    /* calculate error if exact solution has been provided */
    for (int ns = 0; ns < nsims; ns++) {
      CalculateError(&(sim[ns].solver), &(sim[ns].mpi));
    }

    PetscCleanup(&context);

#ifdef with_librom
    context.op_times_arr.push_back(context.waqt);

    for (int ns = 0; ns < nsims; ns++) {
      ResetFilenameIndex( sim[ns].solver.filename_index,
                          sim[ns].solver.index_length );
    }

    if (((libROMInterface*)context.rom_interface)->mode() == _ROM_MODE_TRAIN_) {

      ((libROMInterface*)context.rom_interface)->train();
      if (!rank) printf("libROM: total training wallclock time: %f (seconds).\n",
                        ((libROMInterface*)context.rom_interface)->trainWallclockTime() );

      double total_rom_predict_time = 0;
      for (int iter = 0; iter < context.op_times_arr.size(); iter++) {

        double waqt = context.op_times_arr[iter];

        ((libROMInterface*)context.rom_interface)->predict(sim, waqt);
        if (!rank) printf(  "libROM: Predicted solution at time %1.4e using ROM, wallclock time: %f.\n",
                            waqt, ((libROMInterface*)context.rom_interface)->predictWallclockTime() );
        total_rom_predict_time += ((libROMInterface*)context.rom_interface)->predictWallclockTime();

        /* calculate diff between ROM and PDE solutions */
        if (iter == (context.op_times_arr.size()-1)) {
          if (!rank) printf("libROM:   Calculating diff between PDE and ROM solutions.\n");
          for (int ns = 0; ns < nsims; ns++) {
            CalculateROMDiff(  &(sim[ns].solver),
                               &(sim[ns].mpi) );
          }
        }
        /* write the ROM solution to file */
        OutputROMSolution(sim, nsims, waqt);

      }

      if (!rank) {
        printf( "libROM: total prediction/query wallclock time: %f (seconds).\n",
                total_rom_predict_time );
      }

      ((libROMInterface*)context.rom_interface)->saveROM();

    } else {

      for (int ns = 0; ns < nsims; ns++) {
        sim[ns].solver.rom_diff_norms[0]
          = sim[ns].solver.rom_diff_norms[1]
          = sim[ns].solver.rom_diff_norms[2]
          = -1;
      }

    }

  } else if (context.rom_mode == _ROM_MODE_PREDICT_) {

    for (int ns = 0; ns < nsims; ns++) {
      sim[ns].solver.rom_diff_norms[0]
        = sim[ns].solver.rom_diff_norms[1]
        = sim[ns].solver.rom_diff_norms[2]
        = -1;
      strcpy(sim[ns].solver.ConservationCheck,"no");
    }

    ((libROMInterface*)context.rom_interface)->loadROM();
    ((libROMInterface*)context.rom_interface)->projectInitialSolution(sim);

    {
      int start_iter = sim[0].solver.restart_iter;
      int n_iter = sim[0].solver.n_iter;
      double dt = sim[0].solver.dt;

      double cur_time = start_iter * dt;
      context.op_times_arr.push_back(cur_time);

      for (int iter = start_iter; iter < n_iter; iter++) {
        cur_time += dt;
        if (    ( (iter+1)%sim[0].solver.file_op_iter == 0)
            &&  ( (iter+1) < n_iter) ) {
          context.op_times_arr.push_back(cur_time);
        }
      }

      double t_final = n_iter*dt;
      context.op_times_arr.push_back(t_final);
    }

    double total_rom_predict_time = 0;
    for (int iter = 0; iter < context.op_times_arr.size(); iter++) {

      double waqt = context.op_times_arr[iter];

      ((libROMInterface*)context.rom_interface)->predict(sim, waqt);
      if (!rank) printf(  "libROM: Predicted solution at time %1.4e using ROM, wallclock time: %f.\n",
                          waqt, ((libROMInterface*)context.rom_interface)->predictWallclockTime() );
      total_rom_predict_time += ((libROMInterface*)context.rom_interface)->predictWallclockTime();

      /* write the solution to file */
      for (int ns = 0; ns < nsims; ns++) {
        if (sim[ns].solver.PhysicsOutput) {
          sim[ns].solver.PhysicsOutput( &(sim[ns].solver),
                                        &(sim[ns].mpi),
                                        waqt );
        }
      }
      OutputSolution(sim, nsims, waqt);

    }

    /* calculate error if exact solution has been provided */
    for (int ns = 0; ns < nsims; ns++) {
      CalculateError(&(sim[ns].solver),
                     &(sim[ns].mpi) );
    }

    if (!rank) {
      printf( "libROM: total prediction/query wallclock time: %f (seconds).\n",
              total_rom_predict_time );
    }

  }

  delete ((libROMInterface*)context.rom_interface);
#endif

  PetscFunctionReturn(0);
}

#endif

/*! @file NavierStokes3DInitialize.c
    @author Debojyoti Ghosh
    @brief Initialization of the physics-related variables and function pointers for the 3D Navier-Stokes system
*/
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <boundaryconditions.h>
#include <physicalmodels/navierstokes3d.h>
#include <mpivars.h>
#include <hypar.h>

double NavierStokes3DComputeCFL        (void*,void*,double,double);

int NavierStokes3DFlux              (double*,double*,int,void*,double);
int NavierStokes3DStiffFlux         (double*,double*,int,void*,double);
int NavierStokes3DNonStiffFlux      (double*,double*,int,void*,double);
int NavierStokes3DRoeAverage        (double*,double*,double*,void*);
int NavierStokes3DParabolicFunction (double*,double*,void*,void*,double);
int NavierStokes3DSource            (double*,double*,void*,void*,double);

int NavierStokes3DJacobian          (double*,double*,void*,int,int,int);
int NavierStokes3DStiffJacobian     (double*,double*,void*,int,int,int);

int NavierStokes3DLeftEigenvectors  (double*,double*,void*,int);
int NavierStokes3DRightEigenvectors (double*,double*,void*,int);

int NavierStokes3DUpwindRoe         (double*,double*,double*,double*,double*,double*,int,void*,double);
int NavierStokes3DUpwindRF          (double*,double*,double*,double*,double*,double*,int,void*,double);
int NavierStokes3DUpwindLLF         (double*,double*,double*,double*,double*,double*,int,void*,double);
int NavierStokes3DUpwindRusanov     (double*,double*,double*,double*,double*,double*,int,void*,double);

int NavierStokes3DUpwinddFRoe       (double*,double*,double*,double*,double*,double*,int,void*,double);
int NavierStokes3DUpwindFdFRoe      (double*,double*,double*,double*,double*,double*,int,void*,double);

int NavierStokes3DUpwindRusanovModified   (double*,double*,double*,double*,double*,double*,int,void*,double);
int NavierStokes3DUpwinddFRusanovModified (double*,double*,double*,double*,double*,double*,int,void*,double);
int NavierStokes3DUpwindFdFRusanovModified(double*,double*,double*,double*,double*,double*,int,void*,double);

int NavierStokes3DGravityField      (void*,void*);
int NavierStokes3DModifiedSolution  (double*,double*,int,void*,void*,double);

int NavierStokes3DIBAdiabatic  (void*,void*,double*,double);
int NavierStokes3DIBIsothermal (void*,void*,double*,double);

int NavierStokes3DPreStep           (double*,void*,void*,double);
int NavierStokes3DIBForces          (void*,void*,double);

#if defined(HAVE_CUDA)
int gpuNavierStokes3DInitialize          (void*,void*);
int gpuNavierStokes3DFlux                (double*,double*,int,void*,double);
int gpuNavierStokes3DParabolicFunction   (double*,double*,void*,void*,double);
int gpuNavierStokes3DSource              (double*,double*,void*,void*,double);
int gpuNavierStokes3DModifiedSolution    (double*,double*,int,void*,void*,double);
int gpuNavierStokes3DUpwindRusanov       (double*,double*,double*,double*,double*,double*,int,void*,double);
int gpuNavierStokes3DUpwindRoe           (double*,double*,double*,double*,double*,double*,int,void*,double);
int gpuNavierStokes3DPreStep             (double*,void*,void*,double);
#endif

/*! Initialize the 3D Navier-Stokes (#NavierStokes3D) module:
    Sets the default parameters, read in and set physics-related parameters,
    and set the physics-related function pointers in #HyPar.

    This file reads the file "physics.inp" that must have the following format:

        begin
            <keyword>   <value>
            <keyword>   <value>
            <keyword>   <value>
            ...
            <keyword>   <value>
        end

    where the list of keywords are:

    Keyword name       | Type                 | Variable                                                                | Default value
    ------------------ | -------------------- | ----------------------------------------------------------------------- | ------------------------
    gamma              | double               | #NavierStokes3D::gamma                                                  | 1.4
    Pr                 | double               | #NavierStokes3D::Pr                                                     | 0.72
    Re                 | double               | #NavierStokes3D::Re                                                     | -1
    Minf               | double               | #NavierStokes3D::Minf                                                   | 1.0
    gravity            | double,double,double | #NavierStokes3D::grav_x,#NavierStokes3D::grav_y,#NavierStokes3D::grav_z | 0.0,0.0,0.0
    rho_ref            | double               | #NavierStokes3D::rho0                                                   | 1.0
    p_ref              | double               | #NavierStokes3D::p0                                                     | 1.0
    HB                 | int                  | #NavierStokes3D::HB                                                     | 1
    R                  | double               | #NavierStokes3D::R                                                      | 1.0
    upwinding          | char[]               | #NavierStokes3D::upw_choice                                             | "roe" (#_ROE_)
    ib_wall_type       | char[]               | #NavierStokes3D::ib_wall_type                                           | "adiabatic" (#_IB_ADIABATIC_)
    ib_ramp_time       | double               | #NavierStokes3D::t_ib_ramp                                              | -1
    ib_ramp_width      | double               | #NavierStokes3D::t_ib_width                                             | 0.05
    ib_ramp_type       | char[]               | #NavierStokes3D::ib_ramp_type                                           | "linear" (#_IB_RAMP_LINEAR_)

    + if "ib_wall_type" (#NavierStokes3D::ib_wall_type) is specified as "isothermal",
      it should be followed by the wall temperature (##NavierStokes3D::T_ib_wall), i.e,

        begin
            ...
            ib_wall_type  isothermal 1.0
            ...
        end

    + If "HB" (#NavierStokes3D::HB) is specified as 3, it should be followed by the the
      Brunt-Vaisala frequency (#NavierStokes3D::N_bv), i.e.

        begin
            ...
            HB      3 0.01
            ...
        end

    \b Note: "physics.inp" is \b optional; if absent, default values will be used.
*/
int NavierStokes3DInitialize( void *s, /*!< Solver object of type #HyPar */
                              void *m  /*!< MPI object of type #MPIVariables */
                            )
{
  HyPar           *solver  = (HyPar*)          s;
  MPIVariables    *mpi     = (MPIVariables*)   m;
  NavierStokes3D  *physics = (NavierStokes3D*) solver->physics;
  int             ferr     = 0;

  static int count = 0;

  if (solver->nvars != _MODEL_NVARS_) {
    fprintf(stderr,"Error in NavierStokes3DInitialize(): nvars has to be %d.\n",_MODEL_NVARS_);
    return(1);
  }
  if (solver->ndims != _MODEL_NDIMS_) {
    fprintf(stderr,"Error in NavierStokes3DInitialize(): ndims has to be %d.\n",_MODEL_NDIMS_);
    return(1);
  }

  /* default values */
  physics->gamma     = 1.4;
  physics->Pr        = 0.72;
  physics->Re        = -1;
  physics->Minf      = 1.0;
  physics->C1        = 1.458e-6;
  physics->C2        = 110.4;
  physics->grav_x    = 0.0;
  physics->grav_y    = 0.0;
  physics->grav_z    = 0.0;
  physics->rho0      = 1.0;
  physics->p0        = 1.0;
  physics->HB        = 1;
  physics->R         = 1.0;
  physics->N_bv      = 0.0;
  physics->T_ib_wall = -DBL_MAX;
  physics->t_ib_ramp = -1.0;
  physics->t_ib_width= 0.05;
  physics->ib_T_tol  = 5;
  strcpy(physics->upw_choice,"roe");
  strcpy(physics->ib_write_surface_data,"yes");
  strcpy(physics->ib_wall_type,"adiabatic");
  strcpy(physics->ib_ramp_type,"linear");

  /* reading physical model specific inputs - all processes */
  if (!mpi->rank) {
    FILE *in;
    if (!count) printf("Reading physical model inputs from file \"physics.inp\".\n");
    in = fopen("physics.inp","r");
    if (!in) printf("Warning: File \"physics.inp\" not found. Using default values.\n");
    else {
      char word[_MAX_STRING_SIZE_];
      ferr = fscanf(in,"%s",word);
      if (ferr != 1) {
        fprintf(stderr, "Read error while reading physics.inp in NavierStokes3DInitialize().\n");
        return 1;
      }
      if (!strcmp(word, "begin")){
        while (strcmp(word, "end")){
          ferr = fscanf(in,"%s",word);
          if (ferr != 1) {
            fprintf(stderr, "Read error while reading physics.inp in NavierStokes3DInitialize().\n");
            return 1;
          }
          if (!strcmp(word, "gamma")) {
            ferr = fscanf(in,"%lf",&physics->gamma);
            if (ferr != 1) {
              fprintf(stderr, "Read error while reading physics.inp in NavierStokes3DInitialize().\n");
              return 1;
            }
          } else if (!strcmp(word,"upwinding")) {
            ferr = fscanf(in,"%s",physics->upw_choice);
            if (ferr != 1) {
              fprintf(stderr, "Read error while reading physics.inp in NavierStokes3DInitialize().\n");
              return 1;
            }
          } else if (!strcmp(word,"Pr")) {
            ferr = fscanf(in,"%lf",&physics->Pr);
            if (ferr != 1) {
              fprintf(stderr, "Read error while reading physics.inp in NavierStokes3DInitialize().\n");
              return 1;
            }
          } else if (!strcmp(word,"Re")) {
            ferr = fscanf(in,"%lf",&physics->Re);
            if (ferr != 1) {
              fprintf(stderr, "Read error while reading physics.inp in NavierStokes3DInitialize().\n");
              return 1;
            }
          } else if (!strcmp(word,"Minf")) {
            ferr = fscanf(in,"%lf",&physics->Minf);
            if (ferr != 1) {
              fprintf(stderr, "Read error while reading physics.inp in NavierStokes3DInitialize().\n");
              return 1;
            }
          } else if (!strcmp(word,"gravity")) {
            ferr = fscanf(in,"%lf",&physics->grav_x);
            if (ferr != 1) {
              fprintf(stderr, "Read error while reading physics.inp in NavierStokes3DInitialize().\n");
              return 1;
            }
            ferr = fscanf(in,"%lf",&physics->grav_y);
            if (ferr != 1) {
              fprintf(stderr, "Read error while reading physics.inp in NavierStokes3DInitialize().\n");
              return 1;
            }
            ferr = fscanf(in,"%lf",&physics->grav_z);
            if (ferr != 1) {
              fprintf(stderr, "Read error while reading physics.inp in NavierStokes3DInitialize().\n");
              return 1;
            }
          } else if (!strcmp(word,"rho_ref")) {
            ferr = fscanf(in,"%lf",&physics->rho0);
            if (ferr != 1) {
              fprintf(stderr, "Read error while reading physics.inp in NavierStokes3DInitialize().\n");
              return 1;
            }
          } else if (!strcmp(word,"p_ref")) {
            ferr = fscanf(in,"%lf",&physics->p0);
            if (ferr != 1) {
              fprintf(stderr, "Read error while reading physics.inp in NavierStokes3DInitialize().\n");
              return 1;
            }
          } else if (!strcmp(word,"HB")) {
            ferr = fscanf(in,"%d",&physics->HB);
            if (ferr != 1) {
              fprintf(stderr, "Read error while reading physics.inp in NavierStokes3DInitialize().\n");
              return 1;
            }
            if (physics->HB==3) {
              ferr = fscanf(in,"%lf",&physics->N_bv);
              if (ferr != 1) {
                fprintf(stderr, "Read error while reading physics.inp in NavierStokes3DInitialize().\n");
                return 1;
              }
            }
          } else if (!strcmp(word,"R")) {
            ferr = fscanf(in,"%lf",&physics->R);
            if (ferr != 1) {
              fprintf(stderr, "Read error while reading physics.inp in NavierStokes3DInitialize().\n");
              return 1;
            }
          } else if (!strcmp(word,"ib_surface_data")) {
            ferr = fscanf(in,"%s",physics->ib_write_surface_data);
            if (ferr != 1) {
              fprintf(stderr, "Read error while reading physics.inp in NavierStokes3DInitialize().\n");
              return 1;
            }
          } else if (!strcmp(word,"ib_wall_type")) {
            ferr = fscanf(in,"%s",physics->ib_wall_type);
            if (ferr != 1) {
              fprintf(stderr, "Read error while reading physics.inp in NavierStokes3DInitialize().\n");
              return 1;
            }
            if (!strcmp(physics->ib_wall_type,_IB_ISOTHERMAL_)) {
              ferr = fscanf(in,"%lf",&physics->T_ib_wall);
              if (ferr != 1) {
                fprintf(stderr, "Read error while reading physics.inp in NavierStokes3DInitialize().\n");
                return 1;
              }
            }
            if (!solver->flag_ib) {
              printf("Warning: in NavierStokes3DInitialize().\n");
              printf("Warning: no immersed body present; specification of ib_wall_type unnecessary.\n");
            }
          } else if (!strcmp(word,"ib_ramp_time")) {
            ferr = fscanf(in,"%lf",&physics->t_ib_ramp);
            if (ferr != 1) {
              fprintf(stderr, "Read error while reading physics.inp in NavierStokes3DInitialize().\n");
              return 1;
            }
            if (!solver->flag_ib) {
              printf("Warning: in NavierStokes3DInitialize().\n");
              printf("Warning: no immersed body present; specification of ib_ramp_time unnecessary.\n");
            }
          } else if (!strcmp(word,"ib_ramp_width")) {
            ferr = fscanf(in,"%lf",&physics->t_ib_width);
            if (ferr != 1) {
              fprintf(stderr, "Read error while reading physics.inp in NavierStokes3DInitialize().\n");
              return 1;
            }
            if (!solver->flag_ib) {
              printf("Warning: in NavierStokes3DInitialize().\n");
              printf("Warning: no immersed body present; specification of ib_ramp_width unnecessary.\n");
            }
          } else if (!strcmp(word,"ib_ramp_type")) {
            ferr = fscanf(in,"%s",physics->ib_ramp_type);
            if (ferr != 1) {
              fprintf(stderr, "Read error while reading physics.inp in NavierStokes3DInitialize().\n");
              return 1;
            }
            if (!solver->flag_ib) {
              printf("Warning: in NavierStokes3DInitialize().\n");
              printf("Warning: no immersed body present; specification of ib_ramp_type unnecessary.\n");
            }
          } else if (!strcmp(word,"ib_T_tolerance")) {
            ferr = fscanf(in,"%lf",&physics->ib_T_tol);
            if (ferr != 1) {
              fprintf(stderr, "Read error while reading physics.inp in NavierStokes3DInitialize().\n");
              return 1;
            }
            if (!solver->flag_ib) {
              printf("Warning: in NavierStokes3DInitialize().\n");
              printf("Warning: no immersed body present; specification of ib_T_tolerance unnecessary.\n");
            }
          } else if (strcmp(word,"end")) {
            char useless[_MAX_STRING_SIZE_];
            ferr = fscanf(in,"%s",useless); if (ferr != 1) return(ferr);
            printf("Warning: keyword %s in file \"physics.inp\" with value %s not ",word,useless);
            printf("recognized or extraneous. Ignoring.\n");
          }
        }
      } else {
        fprintf(stderr,"Error: Illegal format in file \"physics.inp\".\n");
        return(1);
      }
    }
    fclose(in);
  }

  IERR MPIBroadcast_character (physics->upw_choice            ,_MAX_STRING_SIZE_,0,&mpi->world); CHECKERR(ierr);
  IERR MPIBroadcast_character (physics->ib_write_surface_data ,_MAX_STRING_SIZE_,0,&mpi->world); CHECKERR(ierr);
  IERR MPIBroadcast_character (physics->ib_wall_type          ,_MAX_STRING_SIZE_,0,&mpi->world); CHECKERR(ierr);
  IERR MPIBroadcast_character (physics->ib_ramp_type          ,_MAX_STRING_SIZE_,0,&mpi->world); CHECKERR(ierr);
  IERR MPIBroadcast_double    (&physics->gamma                ,1                ,0,&mpi->world); CHECKERR(ierr);
  IERR MPIBroadcast_double    (&physics->Pr                   ,1                ,0,&mpi->world); CHECKERR(ierr);
  IERR MPIBroadcast_double    (&physics->Re                   ,1                ,0,&mpi->world); CHECKERR(ierr);
  IERR MPIBroadcast_double    (&physics->Minf                 ,1                ,0,&mpi->world); CHECKERR(ierr);
  IERR MPIBroadcast_double    (&physics->grav_x               ,1                ,0,&mpi->world); CHECKERR(ierr);
  IERR MPIBroadcast_double    (&physics->grav_y               ,1                ,0,&mpi->world); CHECKERR(ierr);
  IERR MPIBroadcast_double    (&physics->grav_z               ,1                ,0,&mpi->world); CHECKERR(ierr);
  IERR MPIBroadcast_double    (&physics->rho0                 ,1                ,0,&mpi->world); CHECKERR(ierr);
  IERR MPIBroadcast_double    (&physics->p0                   ,1                ,0,&mpi->world); CHECKERR(ierr);
  IERR MPIBroadcast_double    (&physics->R                    ,1                ,0,&mpi->world); CHECKERR(ierr);
  IERR MPIBroadcast_double    (&physics->N_bv                 ,1                ,0,&mpi->world); CHECKERR(ierr);
  IERR MPIBroadcast_double    (&physics->T_ib_wall            ,1                ,0,&mpi->world); CHECKERR(ierr);
  IERR MPIBroadcast_double    (&physics->t_ib_ramp            ,1                ,0,&mpi->world); CHECKERR(ierr);
  IERR MPIBroadcast_double    (&physics->t_ib_width           ,1                ,0,&mpi->world); CHECKERR(ierr);
  IERR MPIBroadcast_double    (&physics->ib_T_tol             ,1                ,0,&mpi->world); CHECKERR(ierr);
  IERR MPIBroadcast_integer   (&physics->HB                   ,1                ,0,&mpi->world); CHECKERR(ierr);

  /* if file output is disabled in HyPar, respect that */
  if (!strcmp(solver->op_file_format,"none")) {
    if (!strcmp(physics->ib_write_surface_data,"yes")) {
      if (!mpi->rank) {
        printf("Warning from NavierStokes3DInitialize(): solver->op_file_format is set to \"none\", thus ");
        printf("setting physics->ib_write_surface_data to \"no\" (no solution files will be written).\n");
      }
    }
    strcpy(physics->ib_write_surface_data,"no");
  }

  /* Scaling Re by M_inf */
  physics->Re /= physics->Minf;

  /* check that a well-balanced upwinding scheme is being used for cases with gravity */
  if (   ((physics->grav_x != 0.0) || (physics->grav_y != 0.0) || (physics->grav_z != 0.0) )
      && (strcmp(physics->upw_choice,_RUSANOV_)) ) {
    if (!mpi->rank) {
      fprintf(stderr,"Error in NavierStokes3DInitialize: %s upwinding is needed for flows ",_RUSANOV_);
      fprintf(stderr,"with gravitational forces.\n");
    }
    return(1);
  }
  /* check that solver has the correct choice of diffusion formulation, if viscous flow */
  if (strcmp(solver->spatial_type_par,_NC_2STAGE_) && (physics->Re > 0)) {
    if (!mpi->rank) {
      fprintf(stderr,"Error in NavierStokes3DInitialize(): Parabolic term spatial discretization must be \"%s\"\n",_NC_2STAGE_);
    }
    return(1);
  }

  /* initializing physical model-specific functions */
#if defined(HAVE_CUDA)
  if (solver->use_gpu) {
    solver->PreStep     = gpuNavierStokes3DPreStep;
    solver->FFunction   = gpuNavierStokes3DFlux;
    solver->SFunction   = gpuNavierStokes3DSource;
    solver->UFunction   = gpuNavierStokes3DModifiedSolution;
  } else {
#endif
    solver->PreStep               = NavierStokes3DPreStep;
    solver->ComputeCFL            = NavierStokes3DComputeCFL;
    solver->FFunction             = NavierStokes3DFlux;
    solver->SFunction             = NavierStokes3DSource;
    solver->UFunction             = NavierStokes3DModifiedSolution;
    solver->AveragingFunction     = NavierStokes3DRoeAverage;
    solver->GetLeftEigenvectors   = NavierStokes3DLeftEigenvectors;
    solver->GetRightEigenvectors  = NavierStokes3DRightEigenvectors;
#if defined(HAVE_CUDA)
  }
#endif

  if (solver->flag_ib) {
    if (!strcmp(physics->ib_wall_type,_IB_ADIABATIC_)) {
      solver->IBFunction = NavierStokes3DIBAdiabatic;
    } else if (!strcmp(physics->ib_wall_type,_IB_ISOTHERMAL_)) {
      solver->IBFunction = NavierStokes3DIBIsothermal;
    } else {
      fprintf(stderr, "Error in NavierStokes3DInitialize()\n");
      fprintf(stderr, "  invalid value for IB wall type (%s).\n",
              physics->ib_wall_type );
    }
    if (!strcmp(physics->ib_write_surface_data,"yes")) {
      solver->PhysicsOutput = NavierStokes3DIBForces;
    }
  }

#if defined(HAVE_CUDA)
  if (solver->use_gpu) {

    if (!strcmp(solver->SplitHyperbolicFlux,"yes")) {
      if (!mpi->rank) {
        fprintf(stderr,"Error in NavierStokes3DInitialize(): Not available on GPU!");
      }
      return 1;
    } else {
      solver->JFunction = NavierStokes3DJacobian;
      if (!strcmp(physics->upw_choice,_ROE_    )) {
        solver->Upwind = gpuNavierStokes3DUpwindRoe;
      } else if (!strcmp(physics->upw_choice,_RUSANOV_)) {
        solver->Upwind = gpuNavierStokes3DUpwindRusanov;
      } else {
        if (!mpi->rank) {
          fprintf(stderr,"Error in NavierStokes3DInitialize(): %s is not a valid upwinding scheme on GPU. ",
                  physics->upw_choice);
          fprintf(stderr,"Choices are %s and %s.\n",_ROE_,_RUSANOV_);
        }
        return(1);
      }
    }

  } else {
#endif

    if (!strcmp(solver->SplitHyperbolicFlux,"yes")) {
      solver->FdFFunction = NavierStokes3DNonStiffFlux;
      solver->dFFunction  = NavierStokes3DStiffFlux;
      solver->JFunction   = NavierStokes3DStiffJacobian;
      if (!strcmp(physics->upw_choice,_ROE_)) {
        solver->Upwind    = NavierStokes3DUpwindRoe;
        solver->UpwinddF  = NavierStokes3DUpwinddFRoe;
        solver->UpwindFdF = NavierStokes3DUpwindFdFRoe;
      } else if (!strcmp(physics->upw_choice,_RUSANOV_)) {
        solver->Upwind    = NavierStokes3DUpwindRusanovModified;
        solver->UpwinddF  = NavierStokes3DUpwinddFRusanovModified;
        solver->UpwindFdF = NavierStokes3DUpwindFdFRusanovModified;
      } else {
        if (!mpi->rank) {
          fprintf(stderr,"Error in NavierStokes3DInitialize(): %s is not a valid upwinding scheme ",
                  physics->upw_choice);
          fprintf(stderr,"for use with split hyperbolic flux form. Use %s or %s.\n",
                  _ROE_,_RUSANOV_);
        }
        return(1);
      }
    } else {
      solver->JFunction      = NavierStokes3DJacobian;
      if      (!strcmp(physics->upw_choice,_ROE_    )) solver->Upwind = NavierStokes3DUpwindRoe;
      else if (!strcmp(physics->upw_choice,_RF_     )) solver->Upwind = NavierStokes3DUpwindRF;
      else if (!strcmp(physics->upw_choice,_LLF_    )) solver->Upwind = NavierStokes3DUpwindLLF;
      else if (!strcmp(physics->upw_choice,_RUSANOV_)) solver->Upwind = NavierStokes3DUpwindRusanov;
      else {
        if (!mpi->rank) {
          fprintf(stderr,"Error in NavierStokes3DInitialize(): %s is not a valid upwinding scheme. ",
                  physics->upw_choice);
          fprintf(stderr,"Choices are %s, %s, %s, and %s.\n",_ROE_,_RF_,_LLF_,_RUSANOV_);
        }
        return(1);
      }
    }

#if defined(HAVE_CUDA)
  }
#endif

  /* set the value of gamma in all the boundary objects */
  int n;
  DomainBoundary  *boundary = (DomainBoundary*) solver->boundary;
  for (n = 0; n < solver->nBoundaryZones; n++)  boundary[n].gamma = physics->gamma;

  /* finally, hijack the main solver's dissipation function pointer
   * to this model's own function, since it's difficult to express
   * the dissipation terms in the general form                      */
#if defined(HAVE_CUDA)
  if (solver->use_gpu) {
    solver->ParabolicFunction = gpuNavierStokes3DParabolicFunction;
  } else {
#endif
    solver->ParabolicFunction = NavierStokes3DParabolicFunction;
#if defined(HAVE_CUDA)
  }
#endif

  /* allocate array to hold the gravity field */
  int *dim    = solver->dim_local;
  int ghosts  = solver->ghosts;
  int d, size = 1; for (d=0; d<_MODEL_NDIMS_; d++) size *= (dim[d] + 2*ghosts);
  physics->grav_field_f = (double*) calloc (size, sizeof(double));
  physics->grav_field_g = (double*) calloc (size, sizeof(double));
  /* allocate arrays to hold the fast Jacobian for split form of the hyperbolic flux */
  physics->fast_jac     = (double*) calloc (_MODEL_NDIMS_*size*_MODEL_NVARS_*_MODEL_NVARS_,sizeof(double));
  physics->solution     = (double*) calloc (size*_MODEL_NVARS_,sizeof(double));
  /* initialize the gravity fields */
  IERR NavierStokes3DGravityField(solver,mpi); CHECKERR(ierr);

#if defined(HAVE_CUDA)
  if (solver->use_gpu) gpuNavierStokes3DInitialize(s,m);
#endif

  count++;
  return(0);
}

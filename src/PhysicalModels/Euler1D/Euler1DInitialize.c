/*! @file Euler1DInitialize.c
    @author Debojyoti Ghosh
    @brief Initialize the 1D Euler equations module.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <physicalmodels/euler1d.h>
#include <mpivars.h>
#include <hypar.h>

double Euler1DComputeCFL (void*,void*,double,double);
int    Euler1DFlux       (double*,double*,int,void*,double);
int    Euler1DStiffFlux  (double*,double*,int,void*,double);
int    Euler1DSource     (double*,double*,void*,void*,double);

int    Euler1DUpwindRoe     (double*,double*,double*,double*,double*,double*,int,void*,double);
int    Euler1DUpwinddFRoe   (double*,double*,double*,double*,double*,double*,int,void*,double);

int    Euler1DUpwindRF      (double*,double*,double*,double*,double*,double*,int,void*,double);
int    Euler1DUpwinddFRF    (double*,double*,double*,double*,double*,double*,int,void*,double);

int    Euler1DUpwindLLF     (double*,double*,double*,double*,double*,double*,int,void*,double);
int    Euler1DUpwinddFLLF   (double*,double*,double*,double*,double*,double*,int,void*,double);

int    Euler1DUpwindSWFS    (double*,double*,double*,double*,double*,double*,int,void*,double);
int    Euler1DUpwindRusanov (double*,double*,double*,double*,double*,double*,int,void*,double);

int    Euler1DJacobian      (double*,double*,void*,int,int,int);
int    Euler1DStiffJacobian (double*,double*,void*,int,int,int);

int    Euler1DRoeAverage        (double*,double*,double*,void*);
int    Euler1DLeftEigenvectors  (double*,double*,void*,int);
int    Euler1DRightEigenvectors (double*,double*,void*,int);

int    Euler1DGravityField      (void*,void*);
int    Euler1DSourceUpwindLLF   (double*,double*,double*,double*,int,void*,double);
int    Euler1DSourceUpwindRoe   (double*,double*,double*,double*,int,void*,double);

int    Euler1DModifiedSolution  (double*,double*,int,void*,void*,double);
int    Euler1DPreStep           (double*,void*,void*,double);

/*! Function to initialize the 1D inviscid Euler equations (#Euler1D) module:
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

    Keyword name       | Type         | Variable                      | Default value
    ------------------ | ------------ | ----------------------------- | ------------------------
    gamma              | double       | #Euler1D::gamma               | 1.4
    gravity            | double       | #Euler1D::grav                | 0.0
    grav_type          | int          | #Euler1D::grav_type           | 0
    upwinding          | char[]       | #Euler1D::upw_choice          | "roe" (#_ROE_)

    \b Note: "physics.inp" is \b optional; if absent, default values will be used.
*/
int Euler1DInitialize(
                      void *s, /*!< Solver object of type #HyPar */
                      void *m  /*!< Object of type #MPIVariables containing MPI-related info */
                     )
{
  HyPar         *solver  = (HyPar*)         s;
  MPIVariables  *mpi     = (MPIVariables*)  m;
  Euler1D       *physics = (Euler1D*)       solver->physics;
  int           ferr, d;

  static int count = 0;

  if (solver->nvars != _MODEL_NVARS_) {
    fprintf(stderr,"Error in Euler1DInitialize(): nvars has to be %d.\n",_MODEL_NVARS_);
    return(1);
  }
  if (solver->ndims != _MODEL_NDIMS_) {
    fprintf(stderr,"Error in Euler1DInitialize(): ndims has to be %d.\n",_MODEL_NDIMS_);
    return(1);
  }

  /* default values */
  physics->gamma      = 1.4;
  physics->grav       = 0.0;
  physics->grav_type  = 0;
  strcpy(physics->upw_choice,"roe");

  /* reading physical model specific inputs */
  if (!mpi->rank) {
    FILE *in;
    if (!count) printf("Reading physical model inputs from file \"physics.inp\".\n");
    in = fopen("physics.inp","r");
    if (!in) printf("Warning: File \"physics.inp\" not found. Using default values.\n");
    else {
      char word[_MAX_STRING_SIZE_];
      ferr = fscanf(in,"%s",word); if (ferr != 1) return(1);
      if (!strcmp(word, "begin")){
        while (strcmp(word, "end")){
          ferr = fscanf(in,"%s",word); if (ferr != 1) return(1);
          if (!strcmp(word, "gamma")) {
            ferr = fscanf(in,"%lf",&physics->gamma);
            if (ferr != 1) return(1);
          } else if (!strcmp(word, "gravity")) {
            ferr = fscanf(in,"%lf",&physics->grav);
            if (ferr != 1) return(1);
          } else if (!strcmp(word, "gravity_type")) {
            ferr = fscanf(in,"%d",&physics->grav_type);
            if (ferr != 1) return(1);
          } else if (!strcmp(word,"upwinding")) {
            ferr = fscanf(in,"%s",physics->upw_choice);
            if (ferr != 1) return(1);
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

#ifndef serial
  IERR MPIBroadcast_double    (&physics->gamma    ,1,0,&mpi->world);                  CHECKERR(ierr);
  IERR MPIBroadcast_double    (&physics->grav     ,1,0,&mpi->world);                  CHECKERR(ierr);
  IERR MPIBroadcast_integer   (&physics->grav_type,1,0,&mpi->world);                  CHECKERR(ierr);
  IERR MPIBroadcast_character (physics->upw_choice,_MAX_STRING_SIZE_,0,&mpi->world);  CHECKERR(ierr);
#endif

  if ((physics->grav != 0.0) && (strcmp(physics->upw_choice,_LLF_)) && (strcmp(physics->upw_choice,_ROE_))) {
    if (!mpi->rank) {
      fprintf(stderr,"Error in Euler1DInitialize: %s or %s upwinding is needed for flows ",_LLF_,_ROE_);
      fprintf(stderr,"with gravitational forces.\n");
    }
    return(1);
  }

  /* initializing physical model-specific functions */
  solver->PreStep            = Euler1DPreStep;
  solver->ComputeCFL         = Euler1DComputeCFL;
  solver->FFunction          = Euler1DFlux;
  solver->SFunction          = Euler1DSource;
  solver->UFunction          = Euler1DModifiedSolution;
  if      (!strcmp(physics->upw_choice,_ROE_    )) solver->Upwind = Euler1DUpwindRoe;
  else if (!strcmp(physics->upw_choice,_RF_     )) solver->Upwind = Euler1DUpwindRF;
  else if (!strcmp(physics->upw_choice,_LLF_    )) solver->Upwind = Euler1DUpwindLLF;
  else if (!strcmp(physics->upw_choice,_SWFS_   )) solver->Upwind = Euler1DUpwindSWFS;
  else if (!strcmp(physics->upw_choice,_RUSANOV_)) solver->Upwind = Euler1DUpwindRusanov;
  else {
    if (!mpi->rank) fprintf(stderr,"Error in Euler1DInitialize(): %s is not a valid upwinding scheme.\n",
                            physics->upw_choice);
    return(1);
  }
  if (!strcmp(solver->SplitHyperbolicFlux,"yes")) {
    solver->dFFunction = Euler1DStiffFlux;
    solver->JFunction  = Euler1DStiffJacobian;
    if      (!strcmp(physics->upw_choice,_ROE_ )) solver->UpwinddF = Euler1DUpwinddFRoe;
    else if (!strcmp(physics->upw_choice,_RF_  )) solver->UpwinddF = Euler1DUpwinddFRF;
    else if (!strcmp(physics->upw_choice,_LLF_ )) solver->UpwinddF = Euler1DUpwinddFLLF;
    else {
      if (!mpi->rank) {
        fprintf(stderr,"Error in Euler1DInitialize(): %s is not a valid upwinding scheme ",
                physics->upw_choice);
        fprintf(stderr,"when split form of the hyperbolic flux is used. Use %s, %s or %s.\n",
                _ROE_,_RF_,_LLF_);
      }
      return(1);
    }
  } else {
    solver->dFFunction = NULL;
    solver->UpwinddF   = NULL;
    solver->JFunction  = Euler1DJacobian;
  }
  solver->AveragingFunction     = Euler1DRoeAverage;
  solver->GetLeftEigenvectors   = Euler1DLeftEigenvectors;
  solver->GetRightEigenvectors  = Euler1DRightEigenvectors;

  if      (!strcmp(physics->upw_choice,_LLF_ )) physics->SourceUpwind = Euler1DSourceUpwindLLF;
  else if (!strcmp(physics->upw_choice,_ROE_ )) physics->SourceUpwind = Euler1DSourceUpwindRoe;

  /* allocate array to hold the gravity field */
  int *dim    = solver->dim_local;
  int ghosts  = solver->ghosts;
  int size = 1; for (d=0; d<_MODEL_NDIMS_; d++) size *= (dim[d] + 2*ghosts);
  physics->grav_field = (double*) calloc (size, sizeof(double));
  physics->fast_jac   = (double*) calloc (size*_MODEL_NVARS_*_MODEL_NVARS_, sizeof(double));
  physics->solution   = (double*) calloc (size*_MODEL_NVARS_, sizeof(double));
  IERR Euler1DGravityField(solver,mpi); CHECKERR(ierr);

  count++;
  return(0);
}

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mathfunctions.h>
#include <boundaryconditions.h>
#include <physicalmodels/numa2d-cons.h>
#include <mpivars.h>
#include <hypar.h>

double Numa2DConsComputeCFL (void*,void*,double,double);
int    Numa2DConsFlux       (double*,double*,int,void*,double);
int    Numa2DConsStiffFlux  (double*,double*,int,void*,double);
int    Numa2DConsSource     (double*,double*,void*,double);
int    Numa2DConsUpwindRF   (double*,double*,double*,double*,double*,double*,int,void*,double);
int    Numa2DConsRusanov    (double*,double*,double*,double*,double*,double*,int,void*,double);

void   Numa2DConsCalculateStandardAtmosphere_1(void*,double,double*,double*,double*,double*);
void   Numa2DConsCalculateStandardAtmosphere_2(void*,double,double*,double*,double*,double*);

int Numa2DConsInitialize(void *s,void *m)
{
  HyPar           *solver  = (HyPar*)         s;
  MPIVariables    *mpi     = (MPIVariables*)  m; 
  Numa2DCons      *physics = (Numa2DCons*)    solver->physics;
  int             ferr     = 0;

  if (solver->nvars != _MODEL_NVARS_) {
    fprintf(stderr,"Error in Numa2DConsInitialize(): nvars has to be %d.\n",_MODEL_NVARS_);
    return(1);
  }
  if (solver->ndims != _MODEL_NDIMS_) {
    fprintf(stderr,"Error in Numa2DConsInitialize(): ndims has to be %d.\n",_MODEL_NDIMS_);
    return(1);
  }

  /* default values */
  physics->gamma  = 1.4; 
  physics->R      = 287.058;        /* J kg^{-1} K^{-1} */
  physics->Omega  = 7.2921150E-05;  /* rad s^{-1}       */
  physics->g      = 9.8;            /* m s^{-2}         */

  physics->Pref   = 101327.0;       /* N m^{-2}         */
  physics->Tref   = 288.15;         /* Kelvin           */

  strcpy(physics->upwind,_RUSANOV_UPWINDING_);

  /* default choice of initial atmosphere */
  physics->init_atmos = 1;

  /* reading physical model specific inputs - rank 0 */
  if (!mpi->rank) {
    FILE *in;
    printf("Reading physical model inputs from file \"physics.inp\".\n");
    in = fopen("physics.inp","r");
    if (!in) printf("Warning: File \"physics.inp\" not found. Using default values.\n");
    else {
      char word[_MAX_STRING_SIZE_];
      ferr = fscanf(in,"%s",word); if (ferr != 1) return(1);
      if (!strcmp(word, "begin")){
	      while (strcmp(word, "end")){
		      ferr = fscanf(in,"%s",word); if (ferr != 1) return(1);
          if (!strcmp(word, "gamma")) { 
            ferr = fscanf(in,"%lf",&physics->gamma); if (ferr != 1) return(1);
          } else if (!strcmp(word,"R")) {
            ferr = fscanf(in,"%lf",&physics->R); if (ferr != 1) return(1);
          } else if (!strcmp(word,"g")) {
            ferr = fscanf(in,"%lf",&physics->g); if (ferr != 1) return(1);
          } else if (!strcmp(word,"Omega")) {
            ferr = fscanf(in,"%lf",&physics->Omega); if (ferr != 1) return(1);
          } else if (!strcmp(word,"Pref")) {
            ferr = fscanf(in,"%lf",&physics->Pref); if (ferr != 1) return(1);
          } else if (!strcmp(word,"Tref")) {
            ferr = fscanf(in,"%lf",&physics->Tref); if (ferr != 1) return(1);
          } else if (!strcmp(word,"init_atmos")) {
            ferr = fscanf(in,"%d",&physics->init_atmos); if (ferr != 1) return(1);
          } else if (!strcmp(word,"upwinding")) {
            ferr = fscanf(in,"%s",physics->upwind); if (ferr != 1) return(1);
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
  IERR MPIBroadcast_double    (&physics->gamma      ,1                ,0,&mpi->world); CHECKERR(ierr);
  IERR MPIBroadcast_double    (&physics->R          ,1                ,0,&mpi->world); CHECKERR(ierr);
  IERR MPIBroadcast_double    (&physics->g          ,1                ,0,&mpi->world); CHECKERR(ierr);
  IERR MPIBroadcast_double    (&physics->Omega      ,1                ,0,&mpi->world); CHECKERR(ierr);
  IERR MPIBroadcast_double    (&physics->Pref       ,1                ,0,&mpi->world); CHECKERR(ierr);
  IERR MPIBroadcast_double    (&physics->Tref       ,1                ,0,&mpi->world); CHECKERR(ierr);
  IERR MPIBroadcast_integer   (&physics->init_atmos ,1                ,0,&mpi->world); CHECKERR(ierr);
  IERR MPIBroadcast_character ( physics->upwind     ,_MAX_STRING_SIZE_,0,&mpi->world); CHECKERR(ierr);
#endif

  /* calculate the mean hydrostatic atmosphere as a function of altitude */
  if (physics->init_atmos == 1) {
    physics->StandardAtmosphere = Numa2DConsCalculateStandardAtmosphere_1;
  } else if (physics->init_atmos == 2) {
    physics->StandardAtmosphere = Numa2DConsCalculateStandardAtmosphere_2;
  } else {
    if (!mpi->rank) {
      fprintf(stderr,"Error in Numa2DConsInitialize(): invalid choice of initial atmosphere (init_atmos).\n");
      return(1);
    }
  }
  CHECKERR(ierr);

  /* initializing physical model-specific functions */
  if (!strcmp(solver->SplitHyperbolicFlux,"yes")) 
    solver->dFFunction    = Numa2DConsStiffFlux;
  else solver->dFFunction = NULL;
  solver->FFunction       = Numa2DConsFlux;
  solver->ComputeCFL      = Numa2DConsComputeCFL;
  solver->SFunction       = Numa2DConsSource;
  if (!strcmp(physics->upwind,_RUSANOV_UPWINDING_)) solver->Upwind = Numa2DConsRusanov;
  else {
    if (!mpi->rank) fprintf(stderr,"Error in Numa2DConsInitialize(): Invalid choice of upwinding scheme.\n");
    return(1);
  }

  /* set the value of gamma in all the boundary objects */
  int n;
  DomainBoundary  *boundary = (DomainBoundary*) solver->boundary;
  for (n = 0; n < solver->nBoundaryZones; n++)  boundary[n].gamma = physics->gamma;

  return(0);
}

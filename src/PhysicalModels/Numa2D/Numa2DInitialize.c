#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mathfunctions.h>
#include <boundaryconditions.h>
#include <physicalmodels/numa2d.h>
#include <mpivars.h>
#include <hypar.h>

double Numa2DComputeCFL         (void*,void*,double,double);
int    Numa2DFlux               (double*,double*,int,void*,double);
int    Numa2DStiffFlux          (double*,double*,int,void*,double);
int    Numa2DSource             (double*,double*,void*,void*,double);
int    Numa2DParabolicFunction  (double*,double*,void*,void*,double);

int    Numa2DRusanovFlux      (double*,double*,double*,double*,double*,double*,int,void*,double);
int    Numa2DRusanovLinearFlux(double*,double*,double*,double*,double*,double*,int,void*,double);

void   Numa2DCalculateStandardAtmosphere_1(void*,double,double*,double*,double*,double*);
void   Numa2DCalculateStandardAtmosphere_2(void*,double,double*,double*,double*,double*);

int Numa2DInitialize(void *a_s,void *a_m)
{
  HyPar           *solver  = (HyPar*)         a_s;
  MPIVariables    *mpi     = (MPIVariables*)  a_m;
  Numa2D          *physics = (Numa2D*)        solver->m_physics;
  int             ferr     = 0;

  static int count = 0;

  if (solver->m_nvars != _MODEL_NVARS_) {
    fprintf(stderr,"Error in Numa2DInitialize(): nvars has to be %d.\n",_MODEL_NVARS_);
    return(1);
  }
  if (solver->m_ndims != _MODEL_NDIMS_) {
    fprintf(stderr,"Error in Numa2DInitialize(): ndims has to be %d.\n",_MODEL_NDIMS_);
    return(1);
  }

  /* default values */
  physics->m_gamma  = 1.4;
  physics->m_R      = 287.058;        /* J kg^{-1} K^{-1} */
  physics->m_g      = 9.8;            /* a_m a_s^{-2}         */
  physics->m_mu     = 0.0;            /* a_m^2 a_s^{-1}       */

  physics->Pref   = 101327.0;       /* N a_m^{-2}         */
  physics->Tref   = 288.15;         /* Kelvin           */

  strcpy(physics->m_upwind,_RUSANOV_UPWINDING_);

  /* default choice of initial atmosphere */
  physics->m_init_atmos = 1;

  /* reading physical model specific inputs - rank 0 */
  if (!mpi->m_rank) {
    FILE *in;
    if (!count) printf("Reading physical model inputs from file \"physics.inp\".\n");
    in = fopen("physics.inp","r");
    if (!in) printf("Warning: File \"physics.inp\" not found. Using default values.\n");
    else {
      char word[_MAX_STRING_SIZE_];
      ferr = fscanf(in,"%a_s",word); if (ferr != 1) return(1);
      if (!strcmp(word, "begin")){
        while (strcmp(word, "end")){
          ferr = fscanf(in,"%a_s",word); if (ferr != 1) return(1);
          if (!strcmp(word, "gamma")) {
            ferr = fscanf(in,"%lf",&physics->m_gamma); if (ferr != 1) return(1);
          } else if (!strcmp(word,"R")) {
            ferr = fscanf(in,"%lf",&physics->m_R); if (ferr != 1) return(1);
          } else if (!strcmp(word,"g")) {
            ferr = fscanf(in,"%lf",&physics->m_g); if (ferr != 1) return(1);
          } else if (!strcmp(word,"mu")) {
            ferr = fscanf(in,"%lf",&physics->m_mu); if (ferr != 1) return(1);
          } else if (!strcmp(word,"Pref")) {
            ferr = fscanf(in,"%lf",&physics->Pref); if (ferr != 1) return(1);
          } else if (!strcmp(word,"Tref")) {
            ferr = fscanf(in,"%lf",&physics->Tref); if (ferr != 1) return(1);
          } else if (!strcmp(word,"init_atmos")) {
            ferr = fscanf(in,"%d",&physics->m_init_atmos); if (ferr != 1) return(1);
          } else if (!strcmp(word,"upwinding")) {
            ferr = fscanf(in,"%a_s",physics->m_upwind); if (ferr != 1) return(1);
          } else if (strcmp(word,"end")) {
            char useless[_MAX_STRING_SIZE_];
            ferr = fscanf(in,"%a_s",useless); if (ferr != 1) return(ferr);
            printf("Warning: keyword %a_s in file \"physics.inp\" with value %a_s not ",word,useless);
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
  IERR MPIBroadcast_double   (&physics->m_gamma     ,1                ,0,&mpi->m_world);CHECKERR(ierr);
  IERR MPIBroadcast_double   (&physics->m_R         ,1                ,0,&mpi->m_world);CHECKERR(ierr);
  IERR MPIBroadcast_double   (&physics->m_g         ,1                ,0,&mpi->m_world);CHECKERR(ierr);
  IERR MPIBroadcast_double   (&physics->m_mu        ,1                ,0,&mpi->m_world);CHECKERR(ierr);
  IERR MPIBroadcast_double   (&physics->Pref      ,1                ,0,&mpi->m_world);CHECKERR(ierr);
  IERR MPIBroadcast_double   (&physics->Tref      ,1                ,0,&mpi->m_world);CHECKERR(ierr);
  IERR MPIBroadcast_integer  (&physics->m_init_atmos,1                ,0,&mpi->m_world);CHECKERR(ierr);
  IERR MPIBroadcast_character( physics->m_upwind    ,_MAX_STRING_SIZE_,0,&mpi->m_world);CHECKERR(ierr);
#endif

  /* calculate the mean hydrostatic atmosphere as a function of altitude */
  if (physics->m_init_atmos == 1) {
    physics->StandardAtmosphere = Numa2DCalculateStandardAtmosphere_1;
  } else if (physics->m_init_atmos == 2) {
    physics->StandardAtmosphere = Numa2DCalculateStandardAtmosphere_2;
  } else {
    if (!mpi->m_rank) {
      fprintf(stderr,"Error in Numa2DInitialize(): invalid choice of initial atmosphere (init_atmos).\n");
      return(1);
    }
  }
  CHECKERR(ierr);

  /* initializing physical model-specific functions */
  if (!strcmp(solver->m_split_hyperbolic_flux,"yes"))
    solver->dFFunction    = Numa2DStiffFlux;
  else solver->dFFunction = NULL;
  solver->FFunction       = Numa2DFlux;
  solver->ComputeCFL      = Numa2DComputeCFL;
  solver->SFunction       = Numa2DSource;
  if (!strcmp(physics->m_upwind,_RUSANOV_UPWINDING_)) {
    solver->Upwind        = Numa2DRusanovFlux;
    if (!strcmp(solver->m_split_hyperbolic_flux,"yes"))
      solver->UpwinddF    = Numa2DRusanovLinearFlux;
    else solver->UpwinddF = NULL;
  } else {
    if (!mpi->m_rank) fprintf(stderr,"Error in Numa2DInitialize(): Invalid choice of upwinding scheme.\n");
    return(1);
  }

  /* set the value of gamma in all the boundary objects */
  int n;
  DomainBoundary  *boundary = (DomainBoundary*) solver->m_boundary;
  for (n = 0; n < solver->m_n_boundary_zones; n++)  boundary[n].m_gamma = physics->m_gamma;

  /* finally, hijack the main solver'a_s dissipation function pointer
   * to this model'a_s own function, since it'a_s difficult to express
   * the dissipation terms in the general form                      */
  solver->ParabolicFunction = Numa2DParabolicFunction;

  /* check that solver has the correct choice of diffusion formulation */
  if (strcmp(solver->m_spatial_type_par,_NC_2STAGE_)) {
    if (!mpi->m_rank) {
      fprintf(stderr,"Error in Numa2DInitialize(): Parabolic term spatial discretization must be \"%a_s\"\n",_NC_2STAGE_);
    }
    return(1);
  }

  count++;
  return(0);
}

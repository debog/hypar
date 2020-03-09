/*! @file LinearADRInitialize.c
    @author Debojyoti Ghosh
    @brief Initialize the linear advection-diffusion-reaction module
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <bandedmatrix.h>
#include <physicalmodels/linearadr.h>
#include <mpivars.h>
#include <hypar.h>

int    LinearADRAdvectionField    (void*,void*);
double LinearADRComputeCFL        (void*,void*,double,double);
double LinearADRComputeDiffNumber (void*,void*,double,double);
int    LinearADRAdvection         (double*,double*,int,void*,double);
int    LinearADRDiffusionG        (double*,double*,int,void*,double);
int    LinearADRDiffusionH        (double*,double*,int,int,void*,double);
int    LinearADRReaction          ();
int    LinearADRUpwind            (double*,double*,double*,double*,
                                   double*,double*,int,void*,double);
int    LinearADRJacobian          (double*,double*,void*,int,int);
int    LinearADRWriteAdvField     (void*,void*);

/*! Initialize the linear advection-diffusion-reaction physics module - 
    allocate and set physics-related parameters, read physics-related inputs
    from file, and set the physics-related function pointers in #HyPar
*/
int LinearADRInitialize(
                        void *s, /*!< Solver object of type #HyPar */
                        void *m  /*!< Object of type #MPIVariables containing MPI-related info */
                       )
{
  HyPar         *solver  = (HyPar*)         s;
  MPIVariables  *mpi     = (MPIVariables*)  m; 
  LinearADR     *physics = (LinearADR*)     solver->physics;
  int           i,ferr;

  /* default values */
  physics->constant_advection = -1;
  strcpy(physics->adv_filename,"none");
  physics->a = NULL;
  physics->adv_arr_size = -1;

  physics->d = (double*) calloc (solver->ndims*solver->nvars,sizeof(double));
  _ArraySetValue_(physics->d,solver->ndims*solver->nvars,0.0);

  /* reading physical model specific inputs - all processes */
  if (!mpi->rank) {

    FILE *in;
    printf("Reading physical model inputs from file \"physics.inp\".\n");
    in = fopen("physics.inp","r");

    if (!in) {

      fprintf(stderr,"Error: File \"physics.inp\" not found.\n");
      return(1);

    } else {

      char word[_MAX_STRING_SIZE_];
      ferr = fscanf(in,"%s",word); if (ferr != 1) return(1);

      if (!strcmp(word, "begin")) {
	      while (strcmp(word, "end")) {

		      ferr = fscanf(in,"%s",word); if (ferr != 1) return(1);

          if (!strcmp(word, "advection_filename")) {
            if (physics->constant_advection != -1) {
              fprintf(stderr,"Error in LinearADRInitialize():\n");
              fprintf(stderr,"Maybe advection_filename and advection are both specified.\n");
              return 1;
            } else {
              /* read filename for spatially-varying advection field */
              physics->constant_advection = 0;
              physics->adv_arr_size = solver->npoints_local_wghosts*solver->ndims*solver->nvars;
              physics->a = (double*) calloc ( physics->adv_arr_size, sizeof(double) );
              _ArraySetValue_(physics->a, physics->adv_arr_size, 0.0);
              ferr = fscanf(in,"%s",physics->adv_filename); if (ferr != 1) return(ferr);
            }
          } else if (!strcmp(word, "advection")) {
            if (physics->constant_advection != -1) {
              fprintf(stderr,"Error in LinearADRInitialize():\n");
              fprintf(stderr,"Maybe advection_filename and advection are both specified.\n");
              return 1;
            } else {
              /* read advection coefficients */
              physics->constant_advection = 1;
              physics->adv_arr_size = solver->ndims*solver->nvars;
              physics->a = (double*) calloc ( physics->adv_arr_size, sizeof(double) );
              for (i=0; i<solver->ndims*solver->nvars; i++) ferr = fscanf(in,"%lf",&physics->a[i]);
              if (ferr != 1) return(1);
            }
          } else if (!strcmp(word, "diffusion")) {
            /* read diffusion coefficients */
            for (i=0; i<solver->ndims*solver->nvars; i++) ferr = fscanf(in,"%lf",&physics->d[i]);
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
  MPIBroadcast_integer(&physics->constant_advection,1,0,&mpi->world);
  MPIBroadcast_integer(&physics->adv_arr_size,1,0,&mpi->world);
#endif

  if (mpi->rank) {
    physics->a = (double*) calloc (physics->adv_arr_size,sizeof(double));
    _ArraySetValue_(physics->a, physics->adv_arr_size, 0.0);
  }

  if (physics->constant_advection == 1) {
#ifndef serial
    MPIBroadcast_double(physics->a,physics->adv_arr_size,0,&mpi->world);
#endif
  } else if (physics->constant_advection == 0) {
#ifndef serial
    MPIBroadcast_character(physics->adv_filename, _MAX_STRING_SIZE_,0,&mpi->world);
#endif
    if (!mpi->rank) {
      printf("Reading advection field from %s.\n", physics->adv_filename);
    }
    int retval = LinearADRAdvectionField(solver, mpi);
    if (retval) {
      if (!mpi->rank) {
        fprintf(stderr, "Error in LinearADRInitialize():\n");
        fprintf(stderr, "LinearADRAdvectionField() returned with error.\n");
      }
      return retval;
    }
  }

#ifndef serial
  MPIBroadcast_double(physics->d,solver->ndims*solver->nvars,0,&mpi->world);
#endif

  if (!strcmp(solver->SplitHyperbolicFlux,"yes")) {
    if (!mpi->rank) {
      fprintf(stderr,"Error in LinearADRInitialize: This physical model does not have a splitting ");
      fprintf(stderr,"of the hyperbolic term defined.\n");
    }
    return(1);
  }

  /* initializing physical model-specific functions */
  solver->ComputeCFL         = LinearADRComputeCFL;
  solver->ComputeDiffNumber  = LinearADRComputeDiffNumber;
  solver->FFunction          = LinearADRAdvection;
  solver->GFunction          = LinearADRDiffusionG;
  solver->HFunction          = LinearADRDiffusionH;
  solver->SFunction          = LinearADRReaction;
  solver->JFunction          = LinearADRJacobian;
  solver->Upwind             = LinearADRUpwind;

  if (physics->constant_advection == 0) {
    solver->PhysicsOutput = LinearADRWriteAdvField;
  }

  return(0);
}

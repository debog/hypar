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

int    LinearADRAdvectionField    (void*,void*,int,int,int*);
double LinearADRComputeCFL        (void*,void*,double,double);
double LinearADRComputeDiffNumber (void*,void*,double,double);
int    LinearADRAdvection         (double*,double*,int,void*,double);
int    LinearADRDiffusionG        (double*,double*,int,void*,double);
int    LinearADRDiffusionH        (double*,double*,int,int,void*,double);
int    LinearADRReaction          ();
int    LinearADRUpwind            (double*,double*,double*,double*,
                                   double*,double*,int,void*,double);
int    LinearADRCenteredFlux      (double*,double*,double*,double*,
                                   double*,double*,int,void*,double);
int    LinearADRWriteAdvField     (void*,void*,double);

int    LinearADRAdvectionJacobian (double*,double*,void*,int,int,int);
int    LinearADRDiffusionJacobian (double*,double*,void*,int,int);

/*! Initialize the linear advection-diffusion-reaction physics module -
    allocate and set physics-related parameters, read physics-related inputs
    from file, and set the physics-related function pointers in #HyPar

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
    advection_filename | char[]       | #LinearADR::adv_filename      | "none"
    advection          | double[]     | #LinearADR::a                 | 0
    diffusion          | double[]     | #LinearADR::d                 | 0
    centered_flux      | char[]       | #LinearADR::centered_flux     | "no"

    \b Note:
    + "physics.inp" is \b optional; if absent, default values will be used.
    + Please do *not* specify both "advection" and "advection_filename"!
*/
int LinearADRInitialize(void *a_s, /*!< Solver object of type #HyPar */
                        void *a_m  /*!< Object of type #MPIVariables containing MPI-related info */ )
{
  HyPar         *solver  = (HyPar*)         a_s;
  MPIVariables  *mpi     = (MPIVariables*)  a_m;
  LinearADR     *physics = (LinearADR*)     solver->m_physics;
  int           i,ferr;

  static int count = 0;

  /* default values */
  physics->m_constant_advection = -1;
  strcpy(physics->m_adv_filename,"none");
  physics->m_a = NULL;
  physics->m_adv_arr_size = -1;
  strcpy(physics->m_centered_flux,"no");

  physics->m_d = (double*) calloc (solver->m_ndims*solver->m_nvars,sizeof(double));
  _ArraySetValue_(physics->m_d,solver->m_ndims*solver->m_nvars,0.0);

  /* reading physical model specific inputs - all processes */
  if (!mpi->m_rank) {

    FILE *in;
    if (!count) printf("Reading physical model inputs from file \"physics.inp\".\n");
    in = fopen("physics.inp","r");

    if (!in) {

      fprintf(stderr,"Error: File \"physics.inp\" not found.\n");
      return(1);

    } else {

      char word[_MAX_STRING_SIZE_];
      ferr = fscanf(in,"%a_s",word); if (ferr != 1) return(1);
      if (!strcmp(word, "begin")) {
        while (strcmp(word, "end")) {

          ferr = fscanf(in,"%a_s",word); if (ferr != 1) return(1);

          if (!strcmp(word, "advection_filename")) {
            if (physics->m_constant_advection != -1) {
              fprintf(stderr,"Error in LinearADRInitialize():\n");
              fprintf(stderr,"Maybe advection_filename and advection are both specified.\n");
              return 1;
            } else {
              /* read filename for spatially-varying advection field */
              physics->m_constant_advection = 0;
              physics->m_adv_arr_size = solver->m_npoints_local_wghosts*solver->m_ndims*solver->m_nvars;
              physics->m_a = (double*) calloc ( physics->m_adv_arr_size, sizeof(double) );
              _ArraySetValue_(physics->m_a, physics->m_adv_arr_size, 0.0);
              ferr = fscanf(in,"%a_s",physics->m_adv_filename); if (ferr != 1) return(ferr);
            }
          } else if (!strcmp(word, "advection")) {
            if (physics->m_constant_advection != -1) {
              fprintf(stderr,"Error in LinearADRInitialize():\n");
              fprintf(stderr,"Maybe advection_filename and advection are both specified.\n");
              return 1;
            } else {
              /* read advection coefficients */
              physics->m_constant_advection = 1;
              physics->m_adv_arr_size = solver->m_ndims*solver->m_nvars;
              physics->m_a = (double*) calloc ( physics->m_adv_arr_size, sizeof(double) );
              for (i=0; i<solver->m_ndims*solver->m_nvars; i++) ferr = fscanf(in,"%lf",&physics->m_a[i]);
              if (ferr != 1) return(1);
            }
          } else if (!strcmp(word, "diffusion")) {
            /* read diffusion coefficients */
            for (i=0; i<solver->m_ndims*solver->m_nvars; i++) ferr = fscanf(in,"%lf",&physics->m_d[i]);
            if (ferr != 1) return(1);
          } else if (!strcmp(word, "centered_flux")) {
            ferr = fscanf(in, "%a_s", physics->m_centered_flux);
            if (ferr != 1) return(ferr);
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
  MPIBroadcast_integer(&physics->m_constant_advection,1,0,&mpi->m_world);
  MPIBroadcast_integer(&physics->m_adv_arr_size,1,0,&mpi->m_world);
#endif

  if (mpi->m_rank) {
    physics->m_a = (double*) calloc (physics->m_adv_arr_size,sizeof(double));
    _ArraySetValue_(physics->m_a, physics->m_adv_arr_size, 0.0);
  }

  if (physics->m_constant_advection == 1) {
#ifndef serial
    MPIBroadcast_double(physics->m_a,physics->m_adv_arr_size,0,&mpi->m_world);
#endif
  } else if (physics->m_constant_advection == 0) {
#ifndef serial
    MPIBroadcast_character(physics->m_adv_filename, _MAX_STRING_SIZE_,0,&mpi->m_world);
#endif
  }

#ifndef serial
  MPIBroadcast_double(physics->m_d,solver->m_ndims*solver->m_nvars,0,&mpi->m_world);
  MPIBroadcast_character(physics->m_centered_flux, _MAX_STRING_SIZE_,0,&mpi->m_world);
#endif

  if (!strcmp(solver->m_split_hyperbolic_flux,"yes")) {
    if (!mpi->m_rank) {
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
  solver->JFunction          = LinearADRAdvectionJacobian;
  solver->KFunction          = LinearADRDiffusionJacobian;

  if (!strcmp(physics->m_centered_flux,"no")) {
    solver->Upwind = LinearADRUpwind;
  } else {
    solver->Upwind = LinearADRCenteredFlux;
  }

  if (physics->m_constant_advection == 0) {
    solver->PhysicsInput = LinearADRAdvectionField;
    solver->PhysicsOutput = LinearADRWriteAdvField;
  } else {
    solver->PhysicsInput = NULL;
    solver->PhysicsOutput = NULL;
  }

  count++;
  return(0);
}

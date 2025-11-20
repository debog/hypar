/*! @file TemplateModelInitialize.c
    @author [YOUR NAME]
    @brief Initialize the Template Model physics module
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <physicalmodels/template_model.h>
#include <mpivars.h>
#include <hypar.h>

/* Forward declarations of functions to be registered with solver */
double TemplateModelComputeCFL        (void*,void*,double,double);
int    TemplateModelAdvection         (double*,double*,int,void*,double);
int    TemplateModelUpwind            (double*,double*,double*,double*,
                                       double*,double*,int,void*,double);

/* Optional function declarations - uncomment and implement as needed */
/* double TemplateModelComputeDiffNumber (void*,void*,double,double); */
/* int    TemplateModelDiffusion         (double*,double*,int,void*,double); */
/* int    TemplateModelSource            (double*,double*,void*,void*,double); */
/* int    TemplateModelPreStep           (double*,void*,void*,double); */
/* int    TemplateModelPostStep          (double*,void*,void*,double); */

/*!
 * Initialize the Template Model physics module:
 * - Allocate and set physics-related parameters
 * - Read physics-specific inputs from physics.inp
 * - Set physics-related function pointers in #HyPar
*/
int TemplateModelInitialize(void *s, /*!< Solver object of type #HyPar */
                            void *m  /*!< MPI object of type #MPIVariables */
                           )
{
  HyPar          *solver  = (HyPar*)          s;
  MPIVariables   *mpi     = (MPIVariables*)   m;
  TemplateModel  *physics = (TemplateModel*)  solver->m_physics;
  int            ferr, i;

  static int count = 0;

  /* Set default parameter values */
  physics->param1 = 1.0;
  physics->param2 = 0.0;
  physics->field_array = NULL;
  physics->field_size = 0;
  physics->use_feature = 0;
  strcpy(physics->option_string, "default");
  physics->reference_length = 1.0;
  physics->reference_velocity = 1.0;
  physics->reference_time = 1.0;

  /* Reading physical model specific inputs - root process only */
  if (!mpi->rank) {
    FILE *in;
    if (!count) printf("Reading physical model inputs from file \"physics.inp\".\n");
    in = fopen("physics.inp","r");
    if (in) {
      char word[_MAX_STRING_SIZE_];
      ferr = fscanf(in,"%s",word); if (ferr != 1) return(1);
      if (!strcmp(word, "begin")) {
        while (strcmp(word, "end")) {
          ferr = fscanf(in,"%s",word); if (ferr != 1) return(1);

          if (!strcmp(word, "param1")) {
            /* [REPLACE] Read double parameter */
            ferr = fscanf(in,"%lf",&physics->param1);
            if (ferr != 1) return(1);

          } else if (!strcmp(word, "param2")) {
            /* [REPLACE] Read double parameter */
            ferr = fscanf(in,"%lf",&physics->param2);
            if (ferr != 1) return(1);

          } else if (!strcmp(word, "use_feature")) {
            /* [REPLACE] Read integer flag */
            ferr = fscanf(in,"%d",&physics->use_feature);
            if (ferr != 1) return(1);

          } else if (!strcmp(word, "option_string")) {
            /* [REPLACE] Read string option */
            ferr = fscanf(in,"%s",physics->option_string);
            if (ferr != 1) return(1);

          } else if (strcmp(word,"end")) {
            /* Unrecognized keyword */
            char useless[_MAX_STRING_SIZE_];
            ferr = fscanf(in,"%s",useless);
            if (ferr != 1) return(ferr);
            printf("Warning: keyword %s in file \"physics.inp\" with value %s not ",
                   word, useless);
            printf("recognized or extraneous. Ignoring.\n");
          }
        }
      } else {
        fprintf(stderr,"Error: Illegal format in file \"physics.inp\".\n");
        fprintf(stderr,"       File should begin with \"begin\" keyword.\n");
        return(1);
      }
      fclose(in);
    } else {
      if (!count) {
        fprintf(stderr,"Warning: File \"physics.inp\" not found. ");
        fprintf(stderr,"Using default parameter values.\n");
      }
    }
  }

#ifndef serial
  /* Broadcast parameters to all processes */
  IERR MPIBroadcast_double(&physics->param1,1,0,&mpi->world); CHECKERR(ierr);
  IERR MPIBroadcast_double(&physics->param2,1,0,&mpi->world); CHECKERR(ierr);
  IERR MPIBroadcast_integer(&physics->use_feature,1,0,&mpi->world); CHECKERR(ierr);
  IERR MPIBroadcast_character(physics->option_string,_MAX_STRING_SIZE_,0,&mpi->world);
  CHECKERR(ierr);
#endif

  /* Print parameters to screen */
  if (!mpi->rank) {
    printf("Template Model Parameters:\n");
    printf("  param1        = %lf\n", physics->param1);
    printf("  param2        = %lf\n", physics->param2);
    printf("  use_feature   = %d\n", physics->use_feature);
    printf("  option_string = %s\n", physics->option_string);
  }

  /* Check for incompatible solver options */
  if (!strcmp(solver->m_split_hyperbolic_flux,"yes")) {
    if (!mpi->rank) {
      fprintf(stderr,"Error in TemplateModelInitialize: ");
      fprintf(stderr,"This physical model does not have flux splitting defined.\n");
      fprintf(stderr,"Set SplitHyperbolicFlux to \"no\" in solver.inp.\n");
    }
    return(1);
  }

  /* [REPLACE] Add any other validation checks */
  if (physics->param1 <= 0.0) {
    if (!mpi->rank) {
      fprintf(stderr,"Error: param1 must be positive.\n");
    }
    return(1);
  }

  /* [REPLACE] Allocate model-specific arrays if needed */
  /* Example: allocate field array based on grid size
  int size = 1;
  for (i=0; i<solver->m_ndims; i++) {
    size *= (solver->m_dim_local[i] + 2*solver->m_ghosts);
  }
  physics->field_size = size;
  physics->field_array = (double*) calloc(size*solver->m_nvars, sizeof(double));
  */

  /* Register required function pointers with solver */
  solver->ComputeCFL  = TemplateModelComputeCFL;
  solver->FFunction   = TemplateModelAdvection;
  solver->Upwind      = TemplateModelUpwind;

  /* Register optional function pointers - uncomment as needed */
  /* solver->ComputeDiffNumber  = TemplateModelComputeDiffNumber; */
  /* solver->GFunction          = TemplateModelDiffusion; */
  /* solver->SFunction          = TemplateModelSource; */
  /* solver->PreStep            = TemplateModelPreStep; */
  /* solver->PostStep           = TemplateModelPostStep; */

  count++;
  return(0);
}

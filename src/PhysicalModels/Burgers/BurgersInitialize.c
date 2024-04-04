/*! @file BurgersInitialize.c
    @author John Loffeld
    @brief Initialize the Burgers module
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <bandedmatrix.h>
#include <physicalmodels/burgers.h>
#include <mpivars.h>
#include <hypar.h>

double BurgersComputeCFL (void*,void*,double,double);
int    BurgersAdvection  (double*,double*,int,void*,double);
int    BurgersUpwind     (double*,double*,double*,double*,
                          double*,double*,int,void*,double);

/*! Initialize the nonlinear Burgers physics module -
    allocate and set physics-related parameters, read physics-related inputs
    from file, and set the physics-related function pointers in #HyPar
*/
int BurgersInitialize(void *s, /*!< Solver object of type #HyPar */
                      void *m  /*!< Object of type #MPIVariables containing MPI-related info */
                     )
{
  HyPar         *solver  = (HyPar*)        s;
  MPIVariables  *mpi     = (MPIVariables*) m;
  Burgers       *physics = (Burgers*)      solver->physics;
  int           i, ferr;

  static int count = 0;

  /* reading physical model specific inputs - all processes */
  if (!mpi->rank) {
    FILE *in;
    in = fopen("physics.inp","r");
    if (in) {
      if (!count) printf("Reading physical model inputs from file \"physics.inp\".\n");
      char word[_MAX_STRING_SIZE_];
      ferr = fscanf(in,"%s",word); if (ferr != 1) return(1);
      if (!strcmp(word, "begin")){
        while (strcmp(word, "end")){
          ferr = fscanf(in,"%s",word); if (ferr != 1) return(1);
          if (strcmp(word,"end")) {
            char useless[_MAX_STRING_SIZE_];
            ferr = fscanf(in,"%s",useless); if (ferr != 1) return(ferr);
            printf("Warning: keyword %s in file \"physics.inp\" with value %s not ",
                    word, useless);
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

  if (!strcmp(solver->SplitHyperbolicFlux,"yes")) {
    if (!mpi->rank) {
      fprintf(stderr,"Error in BurgersInitialize: This physical model does not have a splitting ");
      fprintf(stderr,"of the hyperbolic term defined.\n");
    }
    return(1);
  }

#ifndef serial
#endif

  /* initializing physical model-specific functions */
  solver->ComputeCFL = BurgersComputeCFL;
  solver->FFunction  = BurgersAdvection;
  solver->Upwind     = BurgersUpwind;

  count++;
  return(0);
}

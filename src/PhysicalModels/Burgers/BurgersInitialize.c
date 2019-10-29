/*! @file BurgersInitialize.c
    @author Debojyoti Ghosh
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

int    BurgersAdvection         (double*,double*,int,void*,double);
int    BurgersUpwind            (double*,double*,double*,double*,
                                   double*,double*,int,void*,double);

int BurgersInitialize(
                        void *s, /*!< Solver object of type #HyPar */
                        void *m  /*!< Object of type #MPIVariables containing MPI-related info */
                       )
{
  HyPar         *solver  = (HyPar*)        s;
  MPIVariables  *mpi     = (MPIVariables*) m; 
  Burgers       *physics = (Burgers*)      solver->physics;
  int           i,ferr;

  if (!strcmp(solver->SplitHyperbolicFlux,"yes")) {
    if (!mpi->rank) {
      fprintf(stderr,"Error in BurgersInitialize: This physical model does not have a splitting ");
      fprintf(stderr,"of the hyperbolic term defined.\n");
    }
    return(1);
  }

  /* initializing physical model-specific functions */
  solver->FFunction          = BurgersAdvection;
  solver->Upwind             = BurgersUpwind;

  return(0);
}

/*! @file MUSCLInitialize.c
    @brief Initialize the 2nd or 3rd order MUSCL scheme
    @author Debojyoti Ghosh
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <limiters.h>
#include <interpolation.h>
#include <mpivars.h>
#include <hypar.h>

/*!
    Initialize the 2nd or 3rd order MUSCL scheme.
*/
int MUSCLInitialize(
                      void *s,  /*!< Solver object of type #HyPar */
                      void *m   /*!< MPI object of type #MPIVariables */
                   )
{
  HyPar           *solver = (HyPar*) s;
  MPIVariables    *mpi    = (MPIVariables*) m;
  MUSCLParameters *muscl  = (MUSCLParameters*) solver->interp;

  /* default values */
  muscl->eps = 1e-3;
  strcpy(muscl->limiter_type, _LIM_GM_);

  if (!mpi->rank) {
    FILE *in;
    int ferr;
    in = fopen("muscl.inp","r");
    if (!in) printf("Warning: File muscl.inp not found. Using default parameters for muscl2/muscl3 scheme.\n");
    else {
      printf("Reading MUSCL parameters from muscl.inp.\n");
      char word[_MAX_STRING_SIZE_];
      ferr = fscanf(in,"%s",word); if (ferr != 1) return(1);
      if (!strcmp(word, "begin")){
        while (strcmp(word, "end")){
          ferr = fscanf(in,"%s",word); if (ferr != 1) return(1);
          if      (!strcmp(word,"epsilon")) { ferr = fscanf(in,"%lf",&muscl->eps);         if (ferr != 1) return(1); }
          else if (!strcmp(word,"limiter")) { ferr = fscanf(in,"%s" ,muscl->limiter_type); if (ferr != 1) return(1); }
          else if (strcmp(word,"end")) {
            char useless[_MAX_STRING_SIZE_];
            ferr = fscanf(in,"%s",useless); if (ferr != 1) return(ferr);
            printf("Warning: keyword %s in file \"muscl.inp\" with value %s not ",word,useless);
            printf("recognized or extraneous. Ignoring.\n");
          }
        }
      } else {
        fprintf(stderr,"Error: Illegal format in file \"muscl.inp\".\n");
        return(1);
      }
      fclose(in);
    }
  }

  MPIBroadcast_double(&muscl->eps, 1, 0, &mpi->world);
  MPIBroadcast_character(muscl->limiter_type, _MAX_STRING_SIZE_, 0, &mpi->world);

  if (!strcmp(muscl->limiter_type, _LIM_GM_)) {
    muscl->LimiterFunction = LimiterGeneralizedMinMod;
  } else if (!strcmp(muscl->limiter_type, _LIM_MM_)) {
    muscl->LimiterFunction = LimiterMinMod;
  } else if (!strcmp(muscl->limiter_type, _LIM_VANLEER_)) {
    muscl->LimiterFunction = LimiterVanLeer;
  } else if (!strcmp(muscl->limiter_type, _LIM_SUPERBEE_)) {
    muscl->LimiterFunction = LimiterSuperBee;
  } else {
    if (!mpi->rank) {
      fprintf(stderr, "Warning: %s is an invalid limiter type. Using default (Generalized MinMod).\n",
              muscl->limiter_type);
    }
    muscl->LimiterFunction = LimiterGeneralizedMinMod;
  }

  return(0);
}

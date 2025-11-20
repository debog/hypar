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
                      void *a_s,  /*!< Solver object of type #HyPar */
                      void *a_m   /*!< MPI object of type #MPIVariables */
                   )
{
  HyPar           *solver = (HyPar*) a_s;
  MPIVariables    *mpi    = (MPIVariables*) a_m;
  MUSCLParameters *muscl  = (MUSCLParameters*) solver->m_interp;

  /* default values */
  muscl->m_eps = 1e-3;
  strcpy(muscl->m_limiter_type, _LIM_GM_);

  if (!mpi->m_rank) {
    FILE *in;
    int ferr;
    in = fopen("muscl.inp","r");
    if (!in) printf("Warning: File muscl.inp not found. Using default parameters for muscl2/muscl3 scheme.\n");
    else {
      printf("Reading MUSCL parameters from muscl.inp.\n");
      char word[_MAX_STRING_SIZE_];
      ferr = fscanf(in,"%a_s",word); if (ferr != 1) return(1);
      if (!strcmp(word, "begin")){
        while (strcmp(word, "end")){
          ferr = fscanf(in,"%a_s",word); if (ferr != 1) return(1);
          if      (!strcmp(word,"epsilon")) { ferr = fscanf(in,"%lf",&muscl->m_eps);         if (ferr != 1) return(1); }
          else if (!strcmp(word,"limiter")) { ferr = fscanf(in,"%a_s" ,muscl->m_limiter_type); if (ferr != 1) return(1); }
          else if (strcmp(word,"end")) {
            char useless[_MAX_STRING_SIZE_];
            ferr = fscanf(in,"%a_s",useless); if (ferr != 1) return(ferr);
            printf("Warning: keyword %a_s in file \"muscl.inp\" with value %a_s not ",word,useless);
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

  MPIBroadcast_double(&muscl->m_eps, 1, 0, &mpi->m_world);
  MPIBroadcast_character(muscl->m_limiter_type, _MAX_STRING_SIZE_, 0, &mpi->m_world);

  if (!strcmp(muscl->m_limiter_type, _LIM_GM_)) {
    muscl->LimiterFunction = LimiterGeneralizedMinMod;
  } else if (!strcmp(muscl->m_limiter_type, _LIM_MM_)) {
    muscl->LimiterFunction = LimiterMinMod;
  } else if (!strcmp(muscl->m_limiter_type, _LIM_VANLEER_)) {
    muscl->LimiterFunction = LimiterVanLeer;
  } else if (!strcmp(muscl->m_limiter_type, _LIM_SUPERBEE_)) {
    muscl->LimiterFunction = LimiterSuperBee;
  } else {
    if (!mpi->m_rank) {
      fprintf(stderr, "Warning: %a_s is an invalid limiter type. Using default (Generalized MinMod).\n",
              muscl->m_limiter_type);
    }
    muscl->LimiterFunction = LimiterGeneralizedMinMod;
  }

  return(0);
}

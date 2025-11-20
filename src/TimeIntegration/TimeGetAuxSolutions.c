/*! @file TimeGetAuxSolutions.c
    @brief Returns any "auxiliary" solutions
    @author Debojyoti Ghosh
*/
#include <stdio.h>
#include <string.h>
#include <hypar.h>
#include <timeintegration.h>

/*!
  Return auxiliary solution: Some time integrators may have the concept of
  auxiliary solutions that they evolve along with the main solution #HyPar::m_u
  (these may be used for error estimation, for example). This function returns
  a pointer to such an auxiliary solution. Note that the auxiliary solution has
  the same dimensions and array layout as the main solution.
  + Call with the final argument \a n less than 0 to get the total number of
    auxiliary solutions (this value will be stored in the argument \a N at
    exit).
  + Call with the final argument \a n less than or equal to zero to get the
    \a n-th auxiliary solution (the argument \a uaux will point to this at
    exit).

  Note that auxiliary solutions are numbered in the C convention: 0,1,...,N-1.

  Time integration methods which use auxiliary solutions currently implemented:
  + General Linear Methods with Global Error Estimators (GLM-GEE - #_GLM_GEE_)
    (see TimeGLMGEE(), TimeGLMGEEInitialize() ).
*/
int TimeGetAuxSolutions(
                          int     *a_N,      /*!< Number of auxiliary solutions */
                          double  **a_uaux,  /*!< Pointer to the array holding the auxiliary solution */
                          void    *a_s,      /*!< Solver object of type #HyPar */
                          int     a_nu,      /*!< Index of the auxiliary solution to return */
                          int     a_ns       /*!< Index of the simulation domain of which the auxiliary solution to return */
                       )
{
  HyPar           *solver = (HyPar*) a_s;
  TimeIntegration *TS     = (TimeIntegration*) solver->m_time_integrator;

  if (a_nu >= 0) {
    if (!strcmp(solver->m_time_scheme,_GLM_GEE_)) {
      GLMGEEParameters *params = (GLMGEEParameters*) solver->m_msti;
      *a_uaux = (TS->m_U[params->r+a_nu] + TS->m_u_offsets[a_ns]);
    }
  } else {
    if (!TS) *a_N = 0;
    else {
      if (!strcmp(solver->m_time_scheme,_GLM_GEE_)) {
        GLMGEEParameters *params = (GLMGEEParameters*) solver->m_msti;
        *a_N = params->r - 1;
      } else *a_N = 0;
    }
  }

  return(0);
}

/*! @file timeintegration.h
    @author Debojyoti Ghosh
    @brief Contains function declarations for time integration
*/

#ifndef _TIME_INTEGRATION_H_
#define _TIME_INTEGRATION_H_

#include <timeintegration_struct.h>

/*! Initialize the explicit Runge-Kutta time-integration method */
int TimeExplicitRKInitialize(char*,char*,void*,void*);
/*! Clean up variables related to the explicit Runge-Kutta time-integration method */
int TimeExplicitRKCleanup   (void*);

/*! Initialize the General Linear Methods with Global Error Estimators (GLM-GEE) */
int TimeGLMGEEInitialize(char*,char*,void*,void*);
/*! Clean up variables relted to General Linear Methods with Global Error Estimators (GLM-GEE) */
int TimeGLMGEECleanup   (void*);

/*! Initialize the time integration */
int TimeInitialize      (void*,int, int, int, void*);
/*! Clean up variables related to time integration */
int TimeCleanup         (void*);
/*! Function called at the beginning of a time step */
int TimePreStep         (void*);
/*! Take one step in time */
int TimeStep            (void*);
/*! Function called at the end of a time step */
int TimePostStep        (void*);
/*! Print time integration related information */
int TimePrintStep       (void*);
/*! Compute/estimate error in solution */
int TimeError           (void*,void*,double*);
/*! Function to get auxiliary solutions if available (for example, in GLM-GEE methods) */
int TimeGetAuxSolutions (int*,double**,void*,int,int);

/*! Take a step in time using the Forward Euler method */
int TimeForwardEuler  (void*);
/*! Take a step in time using the explicit Runge-Kutta method */
int TimeRK            (void*);
/*! Take a step in time using the General Linear Methods with Global Error Estimators */
int TimeGLMGEE        (void*);

#endif

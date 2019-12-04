/*! @file timeintegration_cpp.h
    @author Debojyoti Ghosh
    @brief Contains C++ function declarations for time integration
*/

#ifndef _TIME_INTEGRATION_CPP_H_
#define _TIME_INTEGRATION_CPP_H_

#include <timeintegration_struct.h>

/*! Initialize the time integration */
extern "C" int TimeInitialize      (void*,int, int, int, void*);
/*! Clean up variables related to time integration */
extern "C" int TimeCleanup         (void*);
/*! Function called at the beginning of a time step */
extern "C" int TimePreStep         (void*);
/*! Take one step in time */
extern "C" int TimeStep            (void*);
/*! Function called at the end of a time step */
extern "C" int TimePostStep        (void*);
/*! Print time integration related information */
extern "C" int TimePrintStep       (void*);
/*! Compute/estimate error in solution */
extern "C" int TimeError           (void*,void*,double*);
/*! Function to get auxiliary solutions if available (for example, in GLM-GEE methods) */
extern "C" int TimeGetAuxSolutions (int*,double**,void*,int,int);

#endif

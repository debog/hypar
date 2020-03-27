/*! @file limiters.h
    @brief Definitions for limiter functions used in MUSCL-type reconstruction schemes.
    @author Debojyoti Ghosh
*/

#ifndef _LIMITERS_H_
#define _LIMITERS_H_

/*! Generalized MinMod Limiter: LimiterGeneralizedMinMod() */
#define _LIM_GM_ "gmm"
/*! MinMod Limiter: LimiterMinMod() */
#define _LIM_MM_ "minmod"
/*! van Leer Limiter: LimiterVanLeer() */
#define _LIM_VANLEER_ "vanleer"
/*! Superbee Limiter: LimiterSuperBee() */
#define _LIM_SUPERBEE_ "superbee"

/*! Generalized MinMod Limiter */
double LimiterGeneralizedMinMod(double);
/*! MinMod Limiter */
double LimiterMinMod(double);
/*! van Leer Limiter */
double LimiterVanLeer(double);
/*! Superbee Limiter */
double LimiterSuperBee(double);

#endif

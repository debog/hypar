/*! @file common_cpp.h
    @brief Some common functions used here and there (C++ declarations)
    @author Debojyoti Ghosh
*/

#ifndef _COMMON_CPP_H_
#define _COMMON_CPP_H_

/*! Get a string corresponding to an integer, i.e. 41 gives "00041" if
    \a width is 5, or "41" if \a width is 2, or "1" if \a width is 1.
*/
extern "C" void GetStringFromInteger(int, char*, int);

#endif

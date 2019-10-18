/*! @file common.h
    @brief Some common functions used here and there
    @author Debojyoti Ghosh
*/


/*!
  Get a string corresponding to an integer, i.e. 41 gives "00041" if
  \a width is 5, or "41" if \a width is 2, or "1" if \a width is 1.
*/
void GetStringFromInteger(int, char*, int);

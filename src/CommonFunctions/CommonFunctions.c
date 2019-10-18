/*! @file CommonFunctions.c
    @brief Some common functions used here and there
    @author Debojyoti Ghosh
*/

#include <common.h>


/*!
  Get a string corresponding to an integer, i.e. 41 gives "00041" if
  \a width is 5, or "41" if \a width is 2, or "1" if \a width is 1.
*/
void GetStringFromInteger(
                            int   a,    /*!< the integer to convert to a string */
                            char  *A,   /*!< the string */
                            int   width /*!< desired width of the string */
                         )
{
  int i;
  for (i=0; i<width; i++) {
    char digit = (char) (a%10 + '0'); 
    a /= 10;
    A[width-1-i] = digit;
  }
  return;
}


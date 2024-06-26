/*! @file CommonFunctions.c
    @brief Some common functions used here and there
    @author Debojyoti Ghosh
*/

#include <stdio.h>
#include <basic.h>
#include <common.h>
#include <math.h>
#include <stddef.h>


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
  if (width > _MAX_STRING_SIZE_-1) {
    fprintf(stderr,"Error in GetStringFromInteger(): requested width is larger than _MAX_STRING_SIZE_.\n");
  }
  int i;
  for (i=0; i<width; i++) {
    char digit = (char) (a%10 + '0');
    a /= 10;
    A[width-1-i] = digit;
  }
  A[width] = 0;
  return;
}

void takeLog(double* array, int array_size)
{
    for (int i = 0; i < array_size; ++i) {
        if (array[i] > 0) {
            array[i] = log(array[i]);
        } else {
            fprintf(stderr, "Error: Logarithm of non-positive value (%f) at index %zu\n", array[i], i);
            // Handle error appropriately, for example by setting to NaN, zero, or skipping
            array[i] = NAN; // Set to NaN
        }
    }
  return;
}

void takeExp(double* array, int array_size)
{
    for (int i = 0; i < array_size; ++i) {
      array[i] = exp(array[i]);
    }
  return;
}
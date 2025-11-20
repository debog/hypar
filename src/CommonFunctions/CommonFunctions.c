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
                            int   a_a,    /*!< the integer to convert to a string */
                            char  *a_A,   /*!< the string */
                            int   a_width /*!< desired width of the string */
                         )
{
  if (a_width > _MAX_STRING_SIZE_-1) {
    fprintf(stderr,"Error in GetStringFromInteger(): requested a_width is larger than _MAX_STRING_SIZE_.\n");
  }
  int i;
  for (i=0; i<a_width; i++) {
    char digit = (char) (a_a%10 + '0');
    a_a /= 10;
    a_A[a_width-1-i] = digit;
  }
  a_A[a_width] = 0;
  return;
}

/*!
  Take the natural logarithm of each element of the array
*/
void takeLog( double* const a_array,      /*!< the array to take log of */
              const int     a_array_size  /*!< size of the array */ )
{
    for (size_t i = 0; i < a_array_size; ++i) {
        if (a_array[i] > 0) {
            a_array[i] = log(a_array[i]);
        } else {
            fprintf(stderr, "Error: Logarithm of non-positive value (%f) at index %zu\n", a_array[i], i);
            // Handle error appropriately, for example by setting to NaN, zero, or skipping
            a_array[i] = NAN; // Set to NaN
        }
    }
  return;
}

/*!
  Take the exponential of each element of the array
*/
void takeExp( double* const a_array,      /*!< the array to take log of */
              const int     a_array_size  /*!< size of the array */ )
{
    for (size_t i = 0; i < a_array_size; ++i) {
      a_array[i] = exp(a_array[i]);
    }
  return;
}

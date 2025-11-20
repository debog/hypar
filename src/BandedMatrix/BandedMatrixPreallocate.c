/*! @file BandedMatrixPreallocate.c
    @author Debojyoti Ghosh
    @brief Preallocate memory for a banded matrix object.
*/

#include <stdlib.h>
#include <bandedmatrix.h>

/*! Preallocate memory for a banded matrix object. */
int BandedMatrixPreallocate(
                              void  *a_A,           /*!< Banded matrix object of the type BandedMatrix */
                              int   a_nbands,       /*!< Number of bands */
                              int   a_nrows_local,  /*!< Local number of rows */
                              int   a_BlockSize     /*!< Block size */
                           )
{
  BandedMatrix *B = (BandedMatrix*) a_A;

  B->m_nbands       = a_nbands;
  B->m_nrows_local  = a_nrows_local;
  B->m_BlockSize    = a_BlockSize;

  B->m_ncol = (int*) calloc (a_nrows_local*a_nbands, sizeof(int));
  B->m_nrow = (int*) calloc (a_nrows_local, sizeof(int));
  B->m_data = (double*) calloc (a_nrows_local*a_nbands*a_BlockSize*a_BlockSize, sizeof(double));

  return(0);
}

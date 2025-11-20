/*! @file BandedMatrixPreallocate.c
    @author Debojyoti Ghosh
    @brief Preallocate memory for a banded matrix object.
*/

#include <stdlib.h>
#include <bandedmatrix.h>

/*! Preallocate memory for a banded matrix object. */
int BandedMatrixPreallocate(
                              void  *A,           /*!< Banded matrix object of the type BandedMatrix */
                              int   nbands,       /*!< Number of bands */
                              int   nrows_local,  /*!< Local number of rows */
                              int   BlockSize     /*!< Block size */
                           )
{
  BandedMatrix *B = (BandedMatrix*) A;

  B->m_nbands       = nbands;
  B->m_nrows_local  = nrows_local;
  B->m_BlockSize    = BlockSize;

  B->m_ncol = (int*) calloc (nrows_local*nbands, sizeof(int));
  B->m_nrow = (int*) calloc (nrows_local, sizeof(int));
  B->m_data = (double*) calloc (nrows_local*nbands*BlockSize*BlockSize, sizeof(double));

  return(0);
}

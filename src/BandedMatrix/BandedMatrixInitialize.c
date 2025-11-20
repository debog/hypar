/*! @file BandedMatrixInitialize.c
    @author Debojyoti Ghosh
    @brief Initialize a banded matrix object
*/

#include <stdlib.h>
#include <bandedmatrix.h>

/*! Initialize a newly-created banded matrix object. */
int BandedMatrixInitialize(void *a_A /*!< Banded matrix object of type BandedMatrix*/)
{
  BandedMatrix *B = (BandedMatrix*) a_A;

  B->m_nbands       = 0;
  B->m_nrows_local  = 0;
  B->m_BlockSize    = 0;

  B->m_ncol = NULL;
  B->m_nrow = NULL;
  B->m_data = NULL;

  return(0);
}

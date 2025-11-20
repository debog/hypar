/*! @file BandedMatrixDestroy.c
    @author Debojyoti Ghosh
    @brief Destroy a banded matrix object
*/

#include <stdlib.h>
#include <bandedmatrix.h>

/*! Free up allocations inside a banded matrix object */
int BandedMatrixDestroy(void *A /*!< Banded matrix object of type #BandedMatrix */)
{
  BandedMatrix *B = (BandedMatrix*) A;

  if (B->m_ncol) free(B->m_ncol);
  if (B->m_nrow) free(B->m_nrow);
  if (B->m_data) free(B->m_data);

  return(0);
}

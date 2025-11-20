/*! @file CompactSchemeCleanup.c
    @brief Cleans up allocations specific to compact finite-difference methods
    @author Debojyoti Ghosh
*/

#include <stdlib.h>
#include <interpolation.h>
#include <mpivars.h>
#include <hypar.h>

/*!
    Cleans up all allocations related to the compact finite difference methods.
*/
int CompactSchemeCleanup(void *s /*!< CompactScheme object of type #CompactScheme */ )
{
  CompactScheme *compact   = (CompactScheme*) s;

  if (compact->m_A) free(compact->m_A);
  if (compact->m_B) free(compact->m_B);
  if (compact->m_C) free(compact->m_C);
  if (compact->m_R) free(compact->m_R);

  if (compact->m_sendbuf) free(compact->m_sendbuf);
  if (compact->m_recvbuf) free(compact->m_recvbuf);

  return(0);
}

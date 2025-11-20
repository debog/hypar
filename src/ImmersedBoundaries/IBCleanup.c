/*! @file IBCleanup.c
    @brief Clean up immersed boundaries-related allocations.
    @author Debojyoti Ghosh
*/

#include <stdlib.h>
#include <immersedboundaries.h>

/*!
*/
int IBCleanup(void *a_s /*!< Object of type #ImmersedBoundary */)
{
  ImmersedBoundary *ib = (ImmersedBoundary*) a_s;
  if (!ib) return(0);

  free(ib->m_body->m_surface);
  free(ib->m_body);

  if (ib->m_n_boundary_nodes > 0) free(ib->m_boundary);
  if (ib->m_nfacets_local > 0) free(ib->m_fmap);

  return(0);
}

/*! @file IBComputeBoundingBox.c
    @author Debojyoti Ghosh
    @brief Compute bounding box for an immersed body
*/

#include <immersedboundaries.h>

/*! Compute the bounding box for a given body. */
int IBComputeBoundingBox(Body3D *b /*!< The body */)
{
  b->m_xmin = b->m_xmax = b->m_surface[0].m_x1;
  b->m_ymin = b->m_ymax = b->m_surface[0].m_y1;
  b->m_zmin = b->m_zmax = b->m_surface[0].m_z1;

  int n;
  for (n = 0; n < b->m_nfacets; n++) {
    if (b->m_surface[n].m_x1 < b->m_xmin) b->m_xmin = b->m_surface[n].m_x1;
    if (b->m_surface[n].m_x2 < b->m_xmin) b->m_xmin = b->m_surface[n].m_x2;
    if (b->m_surface[n].m_x3 < b->m_xmin) b->m_xmin = b->m_surface[n].m_x3;

    if (b->m_surface[n].m_y1 < b->m_ymin) b->m_ymin = b->m_surface[n].m_y1;
    if (b->m_surface[n].m_y2 < b->m_ymin) b->m_ymin = b->m_surface[n].m_y2;
    if (b->m_surface[n].m_y3 < b->m_ymin) b->m_ymin = b->m_surface[n].m_y3;

    if (b->m_surface[n].m_z1 < b->m_zmin) b->m_zmin = b->m_surface[n].m_z1;
    if (b->m_surface[n].m_z2 < b->m_zmin) b->m_zmin = b->m_surface[n].m_z2;
    if (b->m_surface[n].m_z3 < b->m_zmin) b->m_zmin = b->m_surface[n].m_z3;

    if (b->m_surface[n].m_x1 > b->m_xmax) b->m_xmax = b->m_surface[n].m_x1;
    if (b->m_surface[n].m_x2 > b->m_xmax) b->m_xmax = b->m_surface[n].m_x2;
    if (b->m_surface[n].m_x3 > b->m_xmax) b->m_xmax = b->m_surface[n].m_x3;

    if (b->m_surface[n].m_y1 > b->m_ymax) b->m_ymax = b->m_surface[n].m_y1;
    if (b->m_surface[n].m_y2 > b->m_ymax) b->m_ymax = b->m_surface[n].m_y2;
    if (b->m_surface[n].m_y3 > b->m_ymax) b->m_ymax = b->m_surface[n].m_y3;

    if (b->m_surface[n].m_z1 > b->m_zmax) b->m_zmax = b->m_surface[n].m_z1;
    if (b->m_surface[n].m_z2 > b->m_zmax) b->m_zmax = b->m_surface[n].m_z2;
    if (b->m_surface[n].m_z3 > b->m_zmax) b->m_zmax = b->m_surface[n].m_z3;
  }
  return(0);
}

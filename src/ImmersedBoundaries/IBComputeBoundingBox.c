/*! @file IBComputeBoundingBox.c
    @author Debojyoti Ghosh
    @brief Compute bounding box for an immersed body
*/

#include <immersedboundaries.h>

/*! Compute the bounding box for a given body. */
int IBComputeBoundingBox(Body3D *a_b /*!< The body */)
{
  a_b->m_xmin = a_b->m_xmax = a_b->m_surface[0].m_x1;
  a_b->m_ymin = a_b->m_ymax = a_b->m_surface[0].m_y1;
  a_b->m_zmin = a_b->m_zmax = a_b->m_surface[0].m_z1;

  int n;
  for (n = 0; n < a_b->m_nfacets; n++) {
    if (a_b->m_surface[n].m_x1 < a_b->m_xmin) a_b->m_xmin = a_b->m_surface[n].m_x1;
    if (a_b->m_surface[n].m_x2 < a_b->m_xmin) a_b->m_xmin = a_b->m_surface[n].m_x2;
    if (a_b->m_surface[n].m_x3 < a_b->m_xmin) a_b->m_xmin = a_b->m_surface[n].m_x3;

    if (a_b->m_surface[n].m_y1 < a_b->m_ymin) a_b->m_ymin = a_b->m_surface[n].m_y1;
    if (a_b->m_surface[n].m_y2 < a_b->m_ymin) a_b->m_ymin = a_b->m_surface[n].m_y2;
    if (a_b->m_surface[n].m_y3 < a_b->m_ymin) a_b->m_ymin = a_b->m_surface[n].m_y3;

    if (a_b->m_surface[n].m_z1 < a_b->m_zmin) a_b->m_zmin = a_b->m_surface[n].m_z1;
    if (a_b->m_surface[n].m_z2 < a_b->m_zmin) a_b->m_zmin = a_b->m_surface[n].m_z2;
    if (a_b->m_surface[n].m_z3 < a_b->m_zmin) a_b->m_zmin = a_b->m_surface[n].m_z3;

    if (a_b->m_surface[n].m_x1 > a_b->m_xmax) a_b->m_xmax = a_b->m_surface[n].m_x1;
    if (a_b->m_surface[n].m_x2 > a_b->m_xmax) a_b->m_xmax = a_b->m_surface[n].m_x2;
    if (a_b->m_surface[n].m_x3 > a_b->m_xmax) a_b->m_xmax = a_b->m_surface[n].m_x3;

    if (a_b->m_surface[n].m_y1 > a_b->m_ymax) a_b->m_ymax = a_b->m_surface[n].m_y1;
    if (a_b->m_surface[n].m_y2 > a_b->m_ymax) a_b->m_ymax = a_b->m_surface[n].m_y2;
    if (a_b->m_surface[n].m_y3 > a_b->m_ymax) a_b->m_ymax = a_b->m_surface[n].m_y3;

    if (a_b->m_surface[n].m_z1 > a_b->m_zmax) a_b->m_zmax = a_b->m_surface[n].m_z1;
    if (a_b->m_surface[n].m_z2 > a_b->m_zmax) a_b->m_zmax = a_b->m_surface[n].m_z2;
    if (a_b->m_surface[n].m_z3 > a_b->m_zmax) a_b->m_zmax = a_b->m_surface[n].m_z3;
  }
  return(0);
}

/*! @file IBNearestFacetNormal.c
    @brief Find the nearest facet for immersed boundary points.
    @author Debojyoti Ghosh
*/

#include <stdio.h>
#include <basic.h>
#include <mathfunctions.h>
#include <mpivars.h>
#include <immersedboundaries.h>

/*! For each immersed boundary point, find the nearest facet (#Facet3D) of the immersed
    body (#ImmersedBoundary::body). The "nearest" facet is the one which is closest to
    the boundary point in terms of the distance along the normal defined for that facet.
    + The function will first try to find the nearest facet for which a line starting
      from the boundary point along the direction defind by the facet normal passes
      through that facet.
    + Failing the above criterion, the function will find the nearest facet.

    \b Note: This function is sensitive to the fact that the normals defined for the
    immersed body are the \b outward normals, i.e., pointing away from the body. If
    this function returns an error, make sure this is true in the STL file the body
    is defined in.
*/
int IBNearestFacetNormal(
                          void    *a_ib, /*!< Immersed boundary object of type #ImmersedBoundary */
                          void    *a_m,  /*!< MPI object of type #MPIVariables */
                          double  *a_X,  /*!< Array of (local) spatial coordinates */
                          double  a_large_distance, /*!< A large distance */
                          int     *a_dim_l, /*!< Integer array of local grid size in each spatial dimension */
                          int     a_ghosts /*!< Number of ghost points */
                        )
{
  ImmersedBoundary  *IB       = (ImmersedBoundary*) a_ib;
  MPIVariables      *mpi      = (MPIVariables*) a_m;
  Body3D            *body     = IB->m_body;
  Facet3D           *surface  = body->m_surface;
  IBNode            *boundary = IB->m_boundary;

  double  eps = IB->m_tolerance;
  int     nb  = IB->m_n_boundary_nodes,
          nf  = body->m_nfacets;

  int i, j, k, dg, n;
  for (dg = 0; dg < nb; dg++) {
    i = boundary[dg].m_i;
    j = boundary[dg].m_j;
    k = boundary[dg].m_k;

    double xp, yp, zp;
    _GetCoordinate_(0,i,a_dim_l,a_ghosts,a_X,xp);
    _GetCoordinate_(1,j,a_dim_l,a_ghosts,a_X,yp);
    _GetCoordinate_(2,k,a_dim_l,a_ghosts,a_X,zp);
    boundary[dg].m_x = xp;
    boundary[dg].m_y = yp;
    boundary[dg].m_z = zp;

    double  dist_min = a_large_distance;
    int     n_min    = -1;

    for (n = 0; n < nf; n++) {
      double x1, x2, x3;
      double y1, y2, y3;
      double z1, z2, z3;
      x1 = surface[n].m_x1;  x2 = surface[n].m_x2;  x3 = surface[n].m_x3;
      y1 = surface[n].m_y1;  y2 = surface[n].m_y2;  y3 = surface[n].m_y3;
      z1 = surface[n].m_z1;  z2 = surface[n].m_z2;  z3 = surface[n].m_z3;
      double dist =   surface[n].m_nx*(xp-surface[n].m_x1)
                    + surface[n].m_ny*(yp-surface[n].m_y1)
                    + surface[n].m_nz*(zp-surface[n].m_z1);
      if (dist > 0)  continue;
      if (absolute(dist) < dist_min) {
        short   is_it_in = 0;
        double  x_int, y_int, z_int;
        x_int = xp - dist * surface[n].m_nx;
        y_int = yp - dist * surface[n].m_ny;
        z_int = zp - dist * surface[n].m_nz;
        if (absolute(surface[n].m_nx) > eps) {
          double den = (z2-z3)*(y1-y3)-(y2-y3)*(z1-z3);
          double l1, l2, l3;
          l1 = ((y2-y3)*(z3-z_int)-(z2-z3)*(y3-y_int)) / den;
          l2 = ((z1-z3)*(y3-y_int)-(y1-y3)*(z3-z_int)) / den;
          l3 = 1 - l1 - l2;
          if ((l1 > -eps) && (l2 > -eps) && (l3 > -eps))  is_it_in = 1;
        } else if (absolute(surface[n].m_ny) > eps) {
          double den = (x2-x3)*(z1-z3)-(z2-z3)*(x1-x3);
          double l1, l2, l3;
          l1 = ((z2-z3)*(x3-x_int)-(x2-x3)*(z3-z_int)) / den;
          l2 = ((x1-x3)*(z3-z_int)-(z1-z3)*(x3-x_int)) / den;
          l3 = 1 - l1 - l2;
          if ((l1 > -eps) && (l2 > -eps) && (l3 > -eps))  is_it_in = 1;
        } else {
          double den = (y2-y3)*(x1-x3)-(x2-x3)*(y1-y3);
          double l1, l2, l3;
          l1 = ((x2-x3)*(y3-y_int)-(y2-y3)*(x3-x_int)) / den;
          l2 = ((y1-y3)*(x3-x_int)-(x1-x3)*(y3-y_int)) / den;
          l3 = 1 - l1 - l2;
          if ((l1 > -eps) && (l2 > -eps) && (l3 > -eps))  is_it_in = 1;
        }
        if (is_it_in) {
          dist_min = absolute(dist);
          n_min = n;
        }
      }
    }
    if (n_min == -1) {
      for (n = 0; n < nf; n++) {
        double dist =   surface[n].m_nx*(xp-surface[n].m_x1)
                      + surface[n].m_ny*(yp-surface[n].m_y1)
                      + surface[n].m_nz*(zp-surface[n].m_z1);
        if (dist > eps)  continue;
        else {
          if (absolute(dist) < dist_min) {
            dist_min = absolute(dist);
            n_min = n;
          }
        }
      }
    }

    if (n_min == -1)  {
      fprintf(stderr,"Error in IBNearestFacetNormal(): no nearest normal found for boundary node (%d,%d,%d) ",i,j,k);
      fprintf(stderr,"on rank %d.\n",mpi->m_rank);
      return(1);
    } else boundary[dg].m_face = &surface[n_min];
  }

  return(0);
}

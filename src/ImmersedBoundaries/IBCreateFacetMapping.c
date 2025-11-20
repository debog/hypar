/*! @file IBCreateFacetMapping.c
    @brief Create facet mapping
    @author Debojyoti Ghosh
*/

#include <stdlib.h>
#include <string.h>
#include <arrayfunctions.h>
#include <mathfunctions.h>
#include <immersedboundaries.h>
#include <mpivars.h>

/*! is x inside the interval [a,b]? */
static inline int isInside(
                            double x, /*!< the value to check for */
                            double a, /*!< small end of the interval */
                            double b  /*!< big end of the interval */
                          )
{
  return ((x >= a) && (x <= b));
}

/*! Given a point in the 3D space (xc, yc, zc), this function finds the
    indices of the 8 grid points that define the grid cell the given
    point is in, as well as the trilinear interpolation coefficients
    for each of the surrounding grid points.
*/
static int interpNodesCoeffs(
                              void    *a_m,       /*!< MPI object of type #MPIVariables */
                              double  a_xc,       /*!< x-coordinate of the point */
                              double  a_yc,       /*!< y-coordinate of the point */
                              double  a_zc,       /*!< z-coordinate of the point */
                              double  *a_x,       /*!< array of x-coordinates of the grid */
                              double  *a_y,       /*!< array of y-coordinates of the grid */
                              double  *a_z,       /*!< array of z-coordinates of the grid */
                              int     *a_dim,     /*!< local dimensions of the grid */
                              int     a_ghosts,   /*!< number of ghost points */
                              char    *a_mode,    /*!< "mode", i.e., #ImmersedBoundary::mode */
                              int     *a_ii,      /*!< i-index of the surrounding node at the high end
                                                     (i.e. smallest i such that a_x[i] > a_xc) */
                              int     *a_jj,      /*!< j-index of the surrounding node at the high end
                                                     (i.e. smallest j such that a_y[j] > a_yc) */
                              int     *a_kk,      /*!< k-index of the surrounding node at the high end
                                                     (i.e. smallest k such that a_z[k] > a_zc) */
                              int     *a_inodes,  /*!< array to store the indices of the surrounding nodes */
                              double  *a_icoeffs  /*!< array to store the interpolation coefficients of the surrounding nodes */
                            )
{
  MPIVariables *mpi = (MPIVariables*) a_m;

  int i, j, k, ic, jc, kc;
  ic = jc = kc = -1;

  double  xmin = 0.5 * (a_x[a_ghosts-1]         + a_x[a_ghosts]),
          xmax = 0.5 * (a_x[a_dim[0]+a_ghosts-1]  + a_x[a_dim[0]+a_ghosts]),
          ymin = 0.5 * (a_y[a_ghosts-1]         + a_y[a_ghosts]),
          ymax = 0.5 * (a_y[a_dim[1]+a_ghosts-1]  + a_y[a_dim[1]+a_ghosts]),
          zmin = 0.5 * (a_z[a_ghosts-1]         + a_z[a_ghosts]),
          zmax = 0.5 * (a_z[a_dim[2]+a_ghosts-1]  + a_z[a_dim[2]+a_ghosts]);

  for (i = 0; i < a_dim[0]+2*a_ghosts-1; i++) {
    if (isInside(a_xc,a_x[i],a_x[i+1])) {
      ic = i;
      break;
    }
  }
  if      (ic <= a_ghosts-1)        ic = a_ghosts;
  else if (ic >= a_dim[0]+a_ghosts-1) ic = a_dim[0]+a_ghosts-2;

  for (j = 0; j < a_dim[1]+2*a_ghosts-1; j++) {
    if (isInside(a_yc,a_y[j],a_y[j+1])) {
      jc = j;
      break;
    }
  }
  if      (jc <= a_ghosts-1)        jc = a_ghosts;
  else if (jc >= a_dim[1]+a_ghosts-1) jc = a_dim[1]+a_ghosts-2;

  for (k = 0; k < a_dim[2]+2*a_ghosts-1; k++) {
    if (isInside(a_zc,a_z[k],a_z[k+1])) {
      kc = k;
      break;
    }
  }
  if      (kc <= a_ghosts-1)        kc = a_ghosts;
  else if (kc >= a_dim[2]+a_ghosts-1) kc = a_dim[2]+a_ghosts-2;

  if      (!strcmp(a_mode,_IB_XY_))  { kc = a_ghosts; a_zc = 0.5*(zmin+zmax); }
  else if (!strcmp(a_mode,_IB_XZ_))  { jc = a_ghosts; a_yc = 0.5*(ymin+ymax); }
  else if (!strcmp(a_mode,_IB_YZ_))  { ic = a_ghosts; a_xc = 0.5*(xmin+xmax); }

  if (ic == -1) {
    fprintf(stderr,"Error in interpNodesCoeffs() (in ImmersedBoundaries/IBCreateFacetMapping.c) on rank %d: ic = -1.\n", mpi->m_rank);
    return(1);
  }
  if (jc == -1) {
    fprintf(stderr,"Error in interpNodesCoeffs() (in ImmersedBoundaries/IBCreateFacetMapping.c) on rank %d: jc = -1.\n", mpi->m_rank);
    return(1);
  }
  if (kc == -1) {
    fprintf(stderr,"Error in interpNodesCoeffs() (in ImmersedBoundaries/IBCreateFacetMapping.c) on rank %d: kc = -1.\n", mpi->m_rank);
    return(1);
  }
  ic++;
  jc++;
  kc++;

  if (a_ii) *a_ii = ic;
  if (a_jj) *a_jj = jc;
  if (a_kk) *a_kk = kc;

  int pc[_IB_NNODES_], index[_IB_NDIMS_];
  index[0]=ic-1-a_ghosts; index[1]=jc-1-a_ghosts; index[2]=kc-1-a_ghosts; _ArrayIndex1D_(_IB_NDIMS_,a_dim,index,a_ghosts,pc[0]);
  index[0]=ic-a_ghosts  ; index[1]=jc-1-a_ghosts; index[2]=kc-1-a_ghosts; _ArrayIndex1D_(_IB_NDIMS_,a_dim,index,a_ghosts,pc[1]);
  index[0]=ic-1-a_ghosts; index[1]=jc-a_ghosts  ; index[2]=kc-1-a_ghosts; _ArrayIndex1D_(_IB_NDIMS_,a_dim,index,a_ghosts,pc[2]);
  index[0]=ic-a_ghosts  ; index[1]=jc-a_ghosts  ; index[2]=kc-1-a_ghosts; _ArrayIndex1D_(_IB_NDIMS_,a_dim,index,a_ghosts,pc[3]);
  index[0]=ic-1-a_ghosts; index[1]=jc-1-a_ghosts; index[2]=kc-a_ghosts  ; _ArrayIndex1D_(_IB_NDIMS_,a_dim,index,a_ghosts,pc[4]);
  index[0]=ic-a_ghosts  ; index[1]=jc-1-a_ghosts; index[2]=kc-a_ghosts  ; _ArrayIndex1D_(_IB_NDIMS_,a_dim,index,a_ghosts,pc[5]);
  index[0]=ic-1-a_ghosts; index[1]=jc-a_ghosts  ; index[2]=kc-a_ghosts  ; _ArrayIndex1D_(_IB_NDIMS_,a_dim,index,a_ghosts,pc[6]);
  index[0]=ic-a_ghosts  ; index[1]=jc-a_ghosts  ; index[2]=kc-a_ghosts  ; _ArrayIndex1D_(_IB_NDIMS_,a_dim,index,a_ghosts,pc[7]);
  _ArrayCopy1D_(pc,a_inodes,_IB_NNODES_);

  double coeffs[_IB_NNODES_];
  TrilinearInterpCoeffs(a_x[ic-1],a_x[ic],a_y[jc-1],a_y[jc],a_z[kc-1],a_z[kc],a_xc,a_yc,a_zc,&coeffs[0]);
  _ArrayCopy1D_(coeffs,a_icoeffs,_IB_NNODES_);

  return(0);
}

/*!
  This function creates a "facet map", i.e., on each MPI rank, it does the following:
  + Makes a list of facets (defining the immersed body surface) that lie within the
    local computational domain of this MPI rank ("local facets").
  + For each local facet, finds and stores the indices of the grid points that
    surround it, as well as the trilinear interpolation coefficients.
  + For each local facet, finds a "near-surface" point, i.e., a point near the surface
    ("near" in terms of the local grid spacing) along the outward surface normal (i.e.,
    outside the body), and finds and stores the indices of the grid points that
    surround it, as well as the trilinear interpolation coefficients.

Note: each MPI rank has a copy of the entire immersed body, i.e., all the facets.
*/
int IBCreateFacetMapping(
                          void    *a_ib,    /*!< Immersed boundary object of type #ImmersedBoundary */
                          void    *a_m,     /*!< MPI object of type #MPIVariables */
                          double  *a_X,     /*!< Array of local spatial coordinates */
                          int     *a_dim,   /*!< Local dimensions */
                          int     a_ghosts  /*!< Number of ghost points */
                        )
{
  ImmersedBoundary  *IB     = (ImmersedBoundary*) a_ib;
  MPIVariables      *mpi    = (MPIVariables*) a_m;
  Body3D            *body   = IB->m_body;
  int               nfacets = body->m_nfacets, n, count, ierr;
  Facet3D           *facets = body->m_surface;

  double  *x = a_X,
          *y = (x + a_dim[0] + 2*a_ghosts),
          *z = (y + a_dim[1] + 2*a_ghosts);

  double  xmin = 0.5 * (x[a_ghosts-1]         + x[a_ghosts]),
          xmax = 0.5 * (x[a_dim[0]+a_ghosts-1]  + x[a_dim[0]+a_ghosts]),
          ymin = 0.5 * (y[a_ghosts-1]         + y[a_ghosts]),
          ymax = 0.5 * (y[a_dim[1]+a_ghosts-1]  + y[a_dim[1]+a_ghosts]),
          zmin = 0.5 * (z[a_ghosts-1]         + z[a_ghosts]),
          zmax = 0.5 * (z[a_dim[2]+a_ghosts-1]  + z[a_dim[2]+a_ghosts]);

  count = 0;
  for (n = 0; n < nfacets; n++) {

    /* find facet centroid */
    double xc, yc, zc;
    xc = (facets[n].m_x1 + facets[n].m_x2 + facets[n].m_x3) / 3.0;
    yc = (facets[n].m_y1 + facets[n].m_y2 + facets[n].m_y3) / 3.0;
    zc = (facets[n].m_z1 + facets[n].m_z2 + facets[n].m_z3) / 3.0;

    if (!strcmp(IB->m_mode,_IB_3D_)) {
      if (isInside(xc,xmin,xmax) && isInside(yc,ymin,ymax) && isInside(zc,zmin,zmax)) count++;
    } else if (!strcmp(IB->m_mode,_IB_XY_)) {
      if (isInside(xc,xmin,xmax) && isInside(yc,ymin,ymax)) count++;
    } else if (!strcmp(IB->m_mode,_IB_XZ_)) {
      if (isInside(xc,xmin,xmax) && isInside(zc,zmin,zmax)) count++;
    } else if (!strcmp(IB->m_mode,_IB_YZ_)) {
      if (isInside(yc,ymin,ymax) && isInside(zc,zmin,zmax)) count++;
    }
  }

  int nfacets_local = count;
  if (nfacets_local > 0) {

    FacetMap *fmap = (FacetMap*) calloc (nfacets_local, sizeof(FacetMap));
    count = 0;

    for (n = 0; n < nfacets; n++) {

      double xc, yc, zc;
      xc = (facets[n].m_x1 + facets[n].m_x2 + facets[n].m_x3) / 3.0;
      yc = (facets[n].m_y1 + facets[n].m_y2 + facets[n].m_y3) / 3.0;
      zc = (facets[n].m_z1 + facets[n].m_z2 + facets[n].m_z3) / 3.0;

      int flag = 0;
      if (!strcmp(IB->m_mode,_IB_3D_)) {
        if (isInside(xc,xmin,xmax) && isInside(yc,ymin,ymax) && isInside(zc,zmin,zmax)) flag = 1;
      } else if (!strcmp(IB->m_mode,_IB_XY_)) {
        if (isInside(xc,xmin,xmax) && isInside(yc,ymin,ymax)) flag = 1;
      } else if (!strcmp(IB->m_mode,_IB_XZ_)) {
        if (isInside(xc,xmin,xmax) && isInside(zc,zmin,zmax)) flag = 1;
      } else if (!strcmp(IB->m_mode,_IB_YZ_)) {
        if (isInside(yc,ymin,ymax) && isInside(zc,zmin,zmax)) flag = 1;
      }

      if (flag == 1) {
        fmap[count].m_facet = facets + n;
        fmap[count].m_index = n;

        fmap[count].m_xc = xc;
        fmap[count].m_yc = yc;
        fmap[count].m_zc = zc;

        int ic,jc, kc;

        ierr = interpNodesCoeffs( mpi,
                                  xc, yc, zc,
                                  x, y, z,
                                  a_dim,
                                  a_ghosts,
                                  IB->m_mode,
                                  &ic, &jc, &kc,
                                  fmap[count].m_interp_nodes,
                                  fmap[count].m_interp_coeffs );
        if (ierr) {
          fprintf(stderr, "Error in IBCreateFacetMapping(): \n");
          fprintf(stderr, "  interpNodesCoeffs() returned with error code %d on rank %d.\n",
                  ierr, mpi->m_rank );
          return(ierr);
        }

        double dx = x[ic] - x[ic-1];
        double dy = y[jc] - y[jc-1];
        double dz = z[kc] - z[kc-1];
        double ds;
        if      (!strcmp(IB->m_mode,_IB_XY_)) ds = min(dx,dy);
        else if (!strcmp(IB->m_mode,_IB_XZ_)) ds = min(dx,dz);
        else if (!strcmp(IB->m_mode,_IB_YZ_)) ds = min(dy,dz);
        else                                ds = min3(dx,dy,dz);

        double nx = fmap[count].m_facet->m_nx;
        double ny = fmap[count].m_facet->m_ny;
        double nz = fmap[count].m_facet->m_nz;

        if (nx == 0.0) nx += IB->m_delta*ds;
        if (ny == 0.0) ny += IB->m_delta*ds;
        if (nz == 0.0) nz += IB->m_delta*ds;

        double xns = xc + sign(nx)*ds;
        double yns = yc + sign(ny)*ds;
        double zns = zc + sign(nz)*ds;

        fmap[count].m_dx = xns - xc;
        fmap[count].m_dy = yns - yc;
        fmap[count].m_dz = zns - zc;

        ierr = interpNodesCoeffs( mpi,
                                  xns, yns, zns,
                                  x, y, z,
                                  a_dim,
                                  a_ghosts,
                                  IB->m_mode,
                                  NULL,NULL,NULL,
                                  fmap[count].m_interp_nodes_ns,
                                  fmap[count].m_interp_coeffs_ns );
        if (ierr) {
          fprintf(stderr, "Error in IBCreateFacetMapping(): \n");
          fprintf(stderr, "  interpNodesCoeffs() returned with error code %d on rank %d.\n",
                  ierr, mpi->m_rank );
          return(ierr);
        }

        count++;
      }
    }

    IB->m_nfacets_local = nfacets_local;
    IB->m_fmap = fmap;

  } else {

    IB->m_nfacets_local = 0;
    IB->m_fmap = NULL;

  }

  return(0);
}

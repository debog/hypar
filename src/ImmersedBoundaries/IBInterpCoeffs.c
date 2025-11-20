/*! @file IBInterpCoeffs.c
    @brief Compute interpolation nodes and coefficients for immersed boundary points.
    @author Debojyoti Ghosh
*/

#include <basic.h>
#include <arrayfunctions.h>
#include <mathfunctions.h>
#include <mpivars.h>
#include <immersedboundaries.h>

/*!
  Compute the interpolation nodes and coefficients for immersed boundary points: For each
  immersed boundary point, do the following:
  + From the immersed boundary point, extend a probe in the direction defined by the outward
    normal of the "nearest" facet (computed in IBNearestFacetNormal()), till the probe tip
    reaches a point in space such that all surrounding (#_IB_NNODES_) grid points are "interior" points,
    i.e., outside the immersed body (they are "interior" to the computational domain).
  + Store the indices of the surrounding grid points, as well as the trilinear interpolation
    coefficients to interpolate a variable from the surrounding points to the probe tip.
*/
int IBInterpCoeffs(
                    void    *a_ib,    /*!< Immersed boundary object of type #ImmersedBoundary */
                    void    *a_m,     /*!< MPI object of type #MPIVariables */
                    double  *a_X,     /*!< Array of (local) spatial coordinates */
                    int     *a_dim_l, /*!< Integer array of local grid size in each spatial dimension */
                    int     a_ghosts, /*!< Number of ghost points */
                    double *a_blank /*!< Blanking array: for grid points within the
                                         body, this value will be set to 0 */
                  )
{
  ImmersedBoundary  *IB       = (ImmersedBoundary*) a_ib;
  MPIVariables      *mpi      = (MPIVariables*) a_m;
  IBNode            *boundary = IB->m_boundary;

  double  eps         = IB->m_tolerance;
  int     maxiter     = IB->m_itr_max,
          n_boundary  = IB->m_n_boundary_nodes;

  int imax        = a_dim_l[0],
      jmax        = a_dim_l[1],
      kmax        = a_dim_l[2];

  int is = mpi->m_is[0],
      js = mpi->m_is[1],
      ks = mpi->m_is[2];

  static int index[_IB_NDIMS_];

  int dg;
  for (dg = 0; dg < n_boundary; dg++) {
    int    i, j, k, p;
    double xb, yb, zb;
    double nx, ny, nz;
    double xx, yy, zz;
    double dx, dy, dz;
    double ds, dist;
    double xtip, ytip, ztip;

    index[0] = i = boundary[dg].m_i;
    index[1] = j = boundary[dg].m_j;
    index[2] = k = boundary[dg].m_k;
    _ArrayIndex1D_(_IB_NDIMS_,a_dim_l,index,a_ghosts,p);

    xb = boundary[dg].m_x;
    yb = boundary[dg].m_y;
    zb = boundary[dg].m_z;

    nx = boundary[dg].m_face->m_nx;
    ny = boundary[dg].m_face->m_ny;
    nz = boundary[dg].m_face->m_nz;
    xx = boundary[dg].m_face->m_x1;
    yy = boundary[dg].m_face->m_y1;
    zz = boundary[dg].m_face->m_z1;

    dist = nx*(xx-xb) + ny*(yy-yb) + nz*(zz-zb);

    double x1, x2, y1, y2, z1, z2;
    _GetCoordinate_(0,(i+1),a_dim_l,a_ghosts,a_X,x1);
    _GetCoordinate_(0,(i-1),a_dim_l,a_ghosts,a_X,x2);
    _GetCoordinate_(1,(j+1),a_dim_l,a_ghosts,a_X,y1);
    _GetCoordinate_(1,(j-1),a_dim_l,a_ghosts,a_X,y2);
    _GetCoordinate_(2,(k+1),a_dim_l,a_ghosts,a_X,z1);
    _GetCoordinate_(2,(k-1),a_dim_l,a_ghosts,a_X,z2);
    dx = 0.5 * (x1 - x2);
    dy = 0.5 * (y1 - y2);
    dz = 0.5 * (z1 - z2);
    ds = min3(dx, dy, dz);

    xtip = xb + dist*nx;
    ytip = yb + dist*ny;
    ztip = zb + dist*nz;

    int is_it_in = 0;
    int iter = 0;
    int itip, jtip, ktip;
    while(!is_it_in && (iter < maxiter)) {
      iter++;
      itip = i;
      jtip = j;
      ktip = k;

      if (xtip > xb)  {
        double xx;
        _GetCoordinate_(0,itip,a_dim_l,a_ghosts,a_X,xx);
        while ((xx < xtip) && (itip < imax+a_ghosts-1)) {
          itip++;
          _GetCoordinate_(0,itip,a_dim_l,a_ghosts,a_X,xx);
        }
      }  else {
        double xx;
        _GetCoordinate_(0,(itip-1),a_dim_l,a_ghosts,a_X,xx);
        while ((xx > xtip) && (itip > -a_ghosts)) {
          itip--;
          _GetCoordinate_(0,(itip-1),a_dim_l,a_ghosts,a_X,xx);
        }
      }

      if (ytip > yb) {
        double yy;
        _GetCoordinate_(1,jtip,a_dim_l,a_ghosts,a_X,yy);
        while ((yy < ytip) && (jtip < jmax+a_ghosts-1)) {
          jtip++;
          _GetCoordinate_(1,jtip,a_dim_l,a_ghosts,a_X,yy);
        }
      } else {
        double yy;
        _GetCoordinate_(1,(jtip-1),a_dim_l,a_ghosts,a_X,yy);
        while ((yy > ytip) && (jtip > -a_ghosts)) {
          jtip--;
          _GetCoordinate_(1,(jtip-1),a_dim_l,a_ghosts,a_X,yy);
        }
      }

      if (ztip > zb) {
        double zz;
        _GetCoordinate_(2,ktip,a_dim_l,a_ghosts,a_X,zz);
        while ((zz < ztip) && (ktip < kmax+a_ghosts-1))  {
          ktip++;
          _GetCoordinate_(2,ktip,a_dim_l,a_ghosts,a_X,zz);
        }
      } else {
        double zz;
        _GetCoordinate_(2,(ktip-1),a_dim_l,a_ghosts,a_X,zz);
        while ((zz > ztip) && (ktip > -a_ghosts))  {
          ktip--;
          _GetCoordinate_(2,(ktip-1),a_dim_l,a_ghosts,a_X,zz);
        }
      }

      int ptip[_IB_NNODES_];
      index[0] = itip  ; index[1] = jtip  ; index[2] = ktip  ; _ArrayIndex1D_(_IB_NDIMS_,a_dim_l,index,a_ghosts,ptip[0]);
      index[0] = itip-1; index[1] = jtip  ; index[2] = ktip  ; _ArrayIndex1D_(_IB_NDIMS_,a_dim_l,index,a_ghosts,ptip[1]);
      index[0] = itip  ; index[1] = jtip-1; index[2] = ktip  ; _ArrayIndex1D_(_IB_NDIMS_,a_dim_l,index,a_ghosts,ptip[2]);
      index[0] = itip  ; index[1] = jtip  ; index[2] = ktip-1; _ArrayIndex1D_(_IB_NDIMS_,a_dim_l,index,a_ghosts,ptip[3]);
      index[0] = itip-1; index[1] = jtip-1; index[2] = ktip  ; _ArrayIndex1D_(_IB_NDIMS_,a_dim_l,index,a_ghosts,ptip[4]);
      index[0] = itip  ; index[1] = jtip-1; index[2] = ktip-1; _ArrayIndex1D_(_IB_NDIMS_,a_dim_l,index,a_ghosts,ptip[5]);
      index[0] = itip-1; index[1] = jtip  ; index[2] = ktip-1; _ArrayIndex1D_(_IB_NDIMS_,a_dim_l,index,a_ghosts,ptip[6]);
      index[0] = itip-1; index[1] = jtip-1; index[2] = ktip-1; _ArrayIndex1D_(_IB_NDIMS_,a_dim_l,index,a_ghosts,ptip[7]);

      int nflow = 0;
      nflow += a_blank[ptip[0]];
      nflow += a_blank[ptip[1]];
      nflow += a_blank[ptip[2]];
      nflow += a_blank[ptip[3]];
      nflow += a_blank[ptip[4]];
      nflow += a_blank[ptip[5]];
      nflow += a_blank[ptip[6]];
      nflow += a_blank[ptip[7]];
      if (nflow == _IB_NNODES_) {
        is_it_in = 1;
      } else if (nflow < _IB_NNODES_) {
        is_it_in = 0;
        xtip += nx*absolute(ds);
        ytip += ny*absolute(ds);
        ztip += nz*absolute(ds);
      } else {
        fprintf(stderr,"Error in IBInterpCoeffs() (Bug in code) - counting interior points surrounding probe tip \n");
        fprintf(stderr,"on rank %d.\n", mpi->m_rank);
        fprintf(stderr,"Value of nflow is %d but can only be positive and <= %d.\n",nflow,_IB_NNODES_);
        return(1);
      }
    }

    if (!is_it_in) {
      fprintf(stderr,"Error in IBInterpCoeffs() on rank %d - interior point not found for immersed boundary point (%d,%d,%d)!\n",
              mpi->m_rank, i, j, k);
      return(1);
    }

    double tlx[2],tly[2],tlz[2];
    _GetCoordinate_(0,(itip-1),a_dim_l,a_ghosts,a_X,tlx[0]);
    _GetCoordinate_(0,(itip  ),a_dim_l,a_ghosts,a_X,tlx[1]);
    _GetCoordinate_(1,(jtip-1),a_dim_l,a_ghosts,a_X,tly[0]);
    _GetCoordinate_(1,(jtip  ),a_dim_l,a_ghosts,a_X,tly[1]);
    _GetCoordinate_(2,(ktip-1),a_dim_l,a_ghosts,a_X,tlz[0]);
    _GetCoordinate_(2,(ktip  ),a_dim_l,a_ghosts,a_X,tlz[1]);

    int ptip[_IB_NNODES_];
    index[0]=itip-1; index[1]=jtip-1; index[2]=ktip-1; _ArrayIndex1D_(_IB_NDIMS_,a_dim_l,index,a_ghosts,ptip[0]);
    index[0]=itip  ; index[1]=jtip-1; index[2]=ktip-1; _ArrayIndex1D_(_IB_NDIMS_,a_dim_l,index,a_ghosts,ptip[1]);
    index[0]=itip-1; index[1]=jtip  ; index[2]=ktip-1; _ArrayIndex1D_(_IB_NDIMS_,a_dim_l,index,a_ghosts,ptip[2]);
    index[0]=itip  ; index[1]=jtip  ; index[2]=ktip-1; _ArrayIndex1D_(_IB_NDIMS_,a_dim_l,index,a_ghosts,ptip[3]);
    index[0]=itip-1; index[1]=jtip-1; index[2]=ktip  ; _ArrayIndex1D_(_IB_NDIMS_,a_dim_l,index,a_ghosts,ptip[4]);
    index[0]=itip  ; index[1]=jtip-1; index[2]=ktip  ; _ArrayIndex1D_(_IB_NDIMS_,a_dim_l,index,a_ghosts,ptip[5]);
    index[0]=itip-1; index[1]=jtip  ; index[2]=ktip  ; _ArrayIndex1D_(_IB_NDIMS_,a_dim_l,index,a_ghosts,ptip[6]);
    index[0]=itip  ; index[1]=jtip  ; index[2]=ktip  ; _ArrayIndex1D_(_IB_NDIMS_,a_dim_l,index,a_ghosts,ptip[7]);
    _ArrayCopy1D_(ptip,boundary[dg].m_interp_nodes,_IB_NNODES_);

    double coeffs[_IB_NNODES_];
    TrilinearInterpCoeffs(tlx[0],tlx[1],tly[0],tly[1],tlz[0],tlz[1],xtip,ytip,ztip,&coeffs[0]);
    _ArrayCopy1D_(coeffs,boundary[dg].m_interp_coeffs,_IB_NNODES_);

    double tipdist = absolute(nx*(xx-xtip) + ny*(yy-ytip) + nz*(zz-ztip));
    boundary[dg].m_interp_node_distance = tipdist;
    boundary[dg].m_surface_distance = absolute(dist);
    if (tipdist < eps) {
      fprintf(stderr,"Warning in IBInterpCoeffs() on rank %d - how can probe tip be on surface? Tipdist = %e\n",
              mpi->m_rank,tipdist);
    }

  }
  return(0);
}

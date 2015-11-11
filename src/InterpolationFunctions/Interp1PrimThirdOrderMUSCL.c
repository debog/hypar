/*! @file Interp1PrimThirdOrderMUSCL.c
    @author Debojyoti Ghosh
    @brief 3rd order MUSCL scheme with Koren's limiter (component-wise application to vectors)
*/

#include <stdio.h>
#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mathfunctions.h>
#include <interpolation.h>
#include <mpivars.h>
#include <hypar.h>

#ifdef with_omp
#include <omp.h>
#endif

#undef  _MINIMUM_GHOSTS_
/*! \def _MINIMUM_GHOSTS_
 * Minimum number of ghost points required for this interpolation 
 * method.
*/
#define _MINIMUM_GHOSTS_ 2

/*! @brief 3rd order MUSCL scheme with Koren's limiter (component-wise) on a uniform grid

    Computes the interpolated values of the first primitive of a function \f${\bf f}\left({\bf u}\right)\f$
    at the interfaces from the cell-centered values of the function using the 3rd order MUSCL scheme with Koren's limiter on a 
    uniform grid. The first primitive is defined as a function \f${\bf h}\left({\bf u}\right)\f$ that satisfies:
    \f{equation}{
      {\bf f}\left({\bf u}\left(x\right)\right) = \frac{1}{\Delta x} \int_{x-\Delta x/2}^{x+\Delta x/2} {\bf h}\left({\bf u}\left(\zeta\right)\right)d\zeta,
    \f}
    where \f$x\f$ is the spatial coordinate along the dimension of the interpolation. This function computes numerical approximation 
    \f$\hat{\bf f}_{j+1/2} \approx {\bf h}_{j+1/2}\f$ as:
    using the 3rd order MUSCL scheme with Koren's limiter as follows:
    \f{equation}{
      \hat{\bf f}_{j+1/2} = {\bf f}_{j-1} + \phi \left[\frac{1}{3}\left({\bf f}_j-{\bf f}_{j-1}\right) + \frac{1}{6}\left({\bf f}_{j-1}-{\bf f}_{j-2}\right)\right]
    \f}
    where
    \f{equation}{
      \phi = \frac {3\left({\bf f}_j-{\bf f}_{j-1}\right)\left({\bf f}_{j-1}-{\bf f}_{j-2}\right) + \epsilon} {2\left[\left({\bf f}_j-{\bf f}_{j-1}\right)-\left({\bf f}_{j-1}-{\bf f}_{j-2}\right)\right]^2 + 3\left({\bf f}_j-{\bf f}_{j-1}\right)\left({\bf f}_{j-1}-{\bf f}_{j-2}\right) + \epsilon}.
    \f}
    and \f$\epsilon\f$ is a small constant (typically \f$10^{-3}\f$).

    \b Implementation \b Notes:
    + The scalar interpolation method is applied to the vector function in a component-wise manner.
    + The method described above corresponds to a left-biased interpolation. The corresponding right-biased
      interpolation can be obtained by reflecting the equations about interface j+1/2.
    + The function computes the interpolant for the entire grid in one call. It loops over all the grid lines along the interpolation direction
      and carries out the 1D interpolation along these grid lines.
    + Location of cell-centers and cell interfaces along the spatial dimension of the interpolation is shown in the following figure:
      @image html chap1_1Ddomain.png
      @image latex chap1_1Ddomain.eps width=0.9\textwidth

    \b Function \b arguments:

    Argument  | Type      | Explanation             
    --------- | --------- | ---------------------------------------------
    fI        | double*   | Array to hold the computed interpolant at the grid interfaces. This array must have the same layout as the solution, but with \b no \b ghost \b points. Its size should be the same as u in all dimensions, except dir (the dimension along which to interpolate) along which it should be larger by 1 (number of interfaces is 1 more than the number of interior cell centers).
    fC        | double*   | Array with the cell-centered values of the flux function \f${\bf f}\left({\bf u}\right)\f$. This array must have the same layout and size as the solution, \b with \b ghost \b points. 
    u         | double*   | The solution array \f${\bf u}\f$ (with ghost points). If the interpolation is characteristic based, this is needed to compute the eigendecomposition. For a multidimensional problem, the layout is as follows: u is a contiguous 1D array of size (nvars*dim[0]*dim[1]*...*dim[D-1]) corresponding to the multi-dimensional solution, with the following ordering - nvars, dim[0], dim[1], ..., dim[D-1], where nvars is the number of solution components (#HyPar::nvars), dim is the local size (#HyPar::dim_local), D is the number of spatial dimensions.
    x         | double*   | The grid array (with ghost points). This is used only by non-uniform-grid interpolation methods. For multidimensional problems, the layout is as follows: x is a contiguous 1D array of size (dim[0]+dim[1]+...+dim[D-1]), with the spatial coordinates along dim[0] stored from 0,...,dim[0]-1, the spatial coordinates along dim[1] stored along dim[0],...,dim[0]+dim[1]-1, and so forth.
    upw       | int       | Upwinding direction: if positive, a left-biased interpolant will be computed; if negative, a right-biased interpolant will be computed. If the interpolation method is central, then this has no effect.
    dir       | int       | Spatial dimension along which to interpolate (eg: 0 for 1D; 0 or 1 for 2D; 0,1 or 2 for 3D)
    s         | void*     | Solver object of type #HyPar: the following variables are needed - #HyPar::ghosts, #HyPar::ndims, #HyPar::nvars, #HyPar::dim_local.
    m         | void*     | MPI object of type #MPIVariables: this is needed only by compact interpolation method that need to solve a global implicit system across MPI ranks.
    uflag     | int       | A flag indicating if the function being interpolated \f${\bf f}\f$ is the solution itself \f${\bf u}\f$ (if 1, \f${\bf f}\left({\bf u}\right) \equiv {\bf u}\f$).

    Reference:
    + van Leer, B., Towards the Ultimate Conservative Difference Scheme. 2: Monotonicity and Conservation Combined in a Second-Order Scheme, 
      J. of Comput. Phys., 14 (4), 1974, pp.361-370, http://dx.doi.org/10.1016/0021-9991(74)90019-9
    + Koren, B., A Robust Upwind Discretization Method for Advection, Diffusion and Source Terms, Centrum voor Wiskunde en Informatica, Amsterdam, 1993
 */
int Interp1PrimThirdOrderMUSCL(
                                double *fI,  /*!< Array of interpolated function values at the interfaces */
                                double *fC,  /*!< Array of cell-centered values of the function \f${\bf f}\left({\bf u}\right)\f$ */
                                double *u,   /*!< Array of cell-centered values of the solution \f${\bf u}\f$ */
                                double *x,   /*!< Grid coordinates */
                                int    upw,  /*!< Upwind direction (left or right biased) */
                                int    dir,  /*!< Spatial dimension along which to interpolation */
                                void   *s,   /*!< Object of type #HyPar containing solver-related variables */
                                void   *m,   /*!< Object of type #MPIVariables containing MPI-related variables */
                                int    uflag /*!< Flag to indicate if \f$f(u) \equiv u\f$, i.e, if the solution is being reconstructed */
                              )
{
  HyPar           *solver = (HyPar*)          s;
  MUSCLParameters *muscl  = (MUSCLParameters*) solver->interp;

  int ghosts = solver->ghosts;
  int ndims  = solver->ndims;
  int nvars  = solver->nvars;
  int *dim   = solver->dim_local;

  /* define some constants */
  double one_third = 1.0/3.0;
  double one_sixth = 1.0/6.0;

  /* create index and bounds for the outer loop, i.e., to loop over all 1D lines along
     dimension "dir"                                                                    */
  int indexC[ndims], indexI[ndims], index_outer[ndims], bounds_outer[ndims], bounds_inter[ndims];
  _ArrayCopy1D_(dim,bounds_outer,ndims); bounds_outer[dir] =  1;
  _ArrayCopy1D_(dim,bounds_inter,ndims); bounds_inter[dir] += 1;
  int N_outer; _ArrayProduct1D_(bounds_outer,ndims,N_outer);

  int i;
  if (upw > 0) {
#pragma omp parallel for schedule(auto) default(shared) private(i,index_outer,indexC,indexI)
    for (i=0; i<N_outer; i++) {
      _ArrayIndexnD_(ndims,i,bounds_outer,index_outer,0);
      _ArrayCopy1D_(index_outer,indexC,ndims);
      _ArrayCopy1D_(index_outer,indexI,ndims);
      for (indexI[dir] = 0; indexI[dir] < dim[dir]+1; indexI[dir]++) {
        int p; _ArrayIndex1D_(ndims,bounds_inter,indexI,0,p);
        int qm1,qm2,qp1,v;
        indexC[dir] = indexI[dir]-2; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qm2);
        indexC[dir] = indexI[dir]-1; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qm1);
        indexC[dir] = indexI[dir]  ; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qp1);
        for (v=0; v<nvars; v++)  {
          /* Defining stencil points */
          double m2, m1, p1;
          m2 = fC[qm2*nvars+v];
          m1 = fC[qm1*nvars+v];
          p1 = fC[qp1*nvars+v];

          double fdiff = p1 - m1;
          double bdiff = m1 - m2;
          double limit =  (3*fdiff*bdiff + muscl->eps) 
                        / (2*(fdiff-bdiff)*(fdiff-bdiff) + 3*fdiff*bdiff + muscl->eps);

          fI[p*nvars+v] = m1 +  limit * (one_third*fdiff + one_sixth*bdiff);
        }
      }
    }
  } else {
#pragma omp parallel for schedule(auto) default(shared) private(i,index_outer,indexC,indexI)
    for (i=0; i<N_outer; i++) {
      _ArrayIndexnD_(ndims,i,bounds_outer,index_outer,0);
      _ArrayCopy1D_(index_outer,indexC,ndims);
      _ArrayCopy1D_(index_outer,indexI,ndims);
      for (indexI[dir] = 0; indexI[dir] < dim[dir]+1; indexI[dir]++) {
        int p; _ArrayIndex1D_(ndims,bounds_inter,indexI,0,p);
        int qm1,qp1,qp2,v;
        indexC[dir] = indexI[dir]-1; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qm1);
        indexC[dir] = indexI[dir]  ; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qp1);
        indexC[dir] = indexI[dir]+1; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qp2);
        for (v=0; v<nvars; v++)  {
          /* Defining stencil points */
          double m1, p1, p2;
          m1 = fC[qm1*nvars+v];
          p1 = fC[qp1*nvars+v];
          p2 = fC[qp2*nvars+v];

          double fdiff = p2 - p1;
          double bdiff = p1 - m1;
          double limit =  (3*fdiff*bdiff + muscl->eps) 
                        / (2*(fdiff-bdiff)*(fdiff-bdiff) + 3*fdiff*bdiff + muscl->eps);

          fI[p*nvars+v] = p1 -  limit * (one_third*fdiff + one_sixth*bdiff);
        }
      }
    }
  }
  
  return(0);
}

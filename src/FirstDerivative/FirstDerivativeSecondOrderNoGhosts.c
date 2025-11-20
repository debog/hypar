/*! @file FirstDerivativeSecondOrderNoGhosts.c
    @author Ping-Hsuan Tsai
    @brief Second order finite-difference approximation to first derivative
*/

#include <stdio.h>
#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <firstderivative.h>

#include <mpivars.h>
typedef MPIVariables  MPIContext;

#undef  _MINIMUM_GHOSTS_
/*! \def _MINIMUM_GHOSTS_
 * Minimum number of ghost points required.
*/
#define _MINIMUM_GHOSTS_ 1

#ifdef with_omp
#include <omp.h>
#endif

/*! Computes the second-order finite-difference approximation to the first derivative
    (\b Note: not divided by the grid spacing):
    \f{equation}{
      \left(\partial f\right)_i = \left\{ \begin{array}{ll} \frac{1}{2}\left(-3f_i+4f_{i+1}-f_{i+2}\right) & i = -g \\ \frac{1}{2}\left( f_{i+1} - f_{i-1} \right) & -g+1 \leq i \leq N+g-2 \\ \frac{1}{2}\left( f_{i-2} -4f_{i-1}+3f_i \right) & i = N+g-1 \end{array}\right.
    \f}
    where \f$i\f$ is the grid index along the spatial dimension of the derivative, \f$g\f$ is the number of ghost points, and \f$N\f$ is the number of grid points (not including the ghost points) in the spatial dimension of the derivative.
    \n\n
    Notes:
    + The first derivative is computed at the grid points or the cell centers.
    + \b Df and \b f are serialized 1D arrays containing the function and its computed
      derivatives on a multi-dimensional grid. The derivative along the specified dimension
      \b dir is computed by looping through all grid lines along \b dir.
    + \b Df and \b f must have the same number of ghost points.
    + Biased stencils are used at the physical boundaries and ghost point values of f are not used. So this
      function can be used when the ghost values at physical boundaries are not filled. The ghost values at
      internal (MPI) boundaries are still needed.
*/
int FirstDerivativeSecondOrderCentralNoGhosts(  double  *a_Df,    /*!< Array to hold the computed first derivative (with ghost points) */
                                                double  *a_f,     /*!< Array containing the grid point function values whose first
                                                                     derivative is to be computed (with ghost points) */
                                                int     a_dir,    /*!< The spatial dimension along which the derivative is computed */
                                                int     a_bias,   /*!< Forward or backward differencing for non-central
                                                                     finite-difference schemes (-1: backward, 1: forward)*/
                                                int     a_ndims,  /*!< Number of spatial/coordinate dimensions */
                                                int     *a_dim,   /*!< Local dimensions */
                                                int     a_ghosts, /*!< Number of ghost points */
                                                int     a_nvars,  /*!< Number of vector components at each grid points */
                                                void* a_m /*!< MPI object of type #MPIContext */ )
{
  MPIContext* mpi = (MPIContext*) a_m;
  int  i, j, v;

  if ((!a_Df) || (!a_f)) {
    fprintf(stderr, "Error in FirstDerivativeSecondOrder(): input arrays not allocated.\n");
    return(1);
  }
  if (a_ghosts < _MINIMUM_GHOSTS_) {
    fprintf(stderr, "Error in FirstDerivativeSecondOrderCentralNoGhosts(): insufficient number of a_ghosts.\n");
    return(1);
  }

  /* create index and bounds for the outer loop, i.e., to loop over all 1D lines along
     dimension "a_dir"                                                                    */
  int indexC[a_ndims], index_outer[a_ndims], bounds_outer[a_ndims];
  _ArrayCopy1D_(a_dim,bounds_outer,a_ndims); bounds_outer[a_dir] =  1;
  int N_outer;  _ArrayProduct1D_(bounds_outer,a_ndims,N_outer);

#pragma omp parallel for schedule(auto) default(shared) private(i,j,v,index_outer,indexC)
  for (j=0; j<N_outer; j++) {
    _ArrayIndexnD_(a_ndims,j,bounds_outer,index_outer,0);
    _ArrayCopy1D_(index_outer,indexC,a_ndims);
    for (i = 0; i < a_dim[a_dir]; i++) {
      int qC, qL, qR;
      indexC[a_dir] = i-1; _ArrayIndex1D_(a_ndims,a_dim,indexC,a_ghosts,qL);
      indexC[a_dir] = i  ; _ArrayIndex1D_(a_ndims,a_dim,indexC,a_ghosts,qC );
      indexC[a_dir] = i+1; _ArrayIndex1D_(a_ndims,a_dim,indexC,a_ghosts,qR);
      for (v=0; v<a_nvars; v++)  a_Df[qC*a_nvars+v] = 0.5 * (a_f[qR*a_nvars+v]-a_f[qL*a_nvars+v]);
    }
  }

  if (mpi->m_ip[a_dir] == 0) {
    /* left physical boundary: overwrite the leftmost value with biased finite-difference */
#pragma omp parallel for schedule(auto) default(shared) private(i,j,v,index_outer,indexC)
    for (j=0; j<N_outer; j++) {
      _ArrayIndexnD_(a_ndims,j,bounds_outer,index_outer,0);
      _ArrayCopy1D_(index_outer,indexC,a_ndims);
      i = 0;
      int qC, qR, qR2;
      indexC[a_dir] = i  ; _ArrayIndex1D_(a_ndims,a_dim,indexC,a_ghosts,qC );
      indexC[a_dir] = i+1; _ArrayIndex1D_(a_ndims,a_dim,indexC,a_ghosts,qR);
      indexC[a_dir] = i+2; _ArrayIndex1D_(a_ndims,a_dim,indexC,a_ghosts,qR2);
      for (v=0; v<a_nvars; v++)  a_Df[qC*a_nvars+v] = (-0.5*a_f[qR2*a_nvars+v]+2*a_f[qR*a_nvars+v]-1.5*a_f[qC*a_nvars+v]);
    }
  }

  if (mpi->m_ip[a_dir] == (mpi->m_iproc[a_dir]-1)) {
    /* right physical boundary: overwrite the rightmost value with biased finite-difference */
#pragma omp parallel for schedule(auto) default(shared) private(i,j,v,index_outer,indexC)
    for (j=0; j<N_outer; j++) {
      _ArrayIndexnD_(a_ndims,j,bounds_outer,index_outer,0);
      _ArrayCopy1D_(index_outer,indexC,a_ndims);
      i = a_dim[a_dir] - 1;
      int qC, qL, qL2;
      indexC[a_dir] = i-2; _ArrayIndex1D_(a_ndims,a_dim,indexC,a_ghosts,qL2);
      indexC[a_dir] = i-1; _ArrayIndex1D_(a_ndims,a_dim,indexC,a_ghosts,qL);
      indexC[a_dir] = i  ; _ArrayIndex1D_(a_ndims,a_dim,indexC,a_ghosts,qC );
      for (v=0; v<a_nvars; v++)  a_Df[qC*a_nvars+v] = (0.5*a_f[qL2*a_nvars+v]-2*a_f[qL*a_nvars+v]+1.5*a_f[qC*a_nvars+v]);
    }
  }

  return 0;
}

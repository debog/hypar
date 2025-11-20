/*! @file SecondDerivativeSecondOrderNoGhosts.c
    @brief 2nd order discretization of the second derivative.
    @author Ping-Hsuan Tsai, Debojyoti Ghosh
*/

#include <stdio.h>
#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <secondderivative.h>

#include <mpivars.h>
typedef MPIVariables  MPIContext;

#undef  _MINIMUM_GHOSTS_
/*! \def _MINIMUM_GHOSTS_
 * Minimum number of ghost points required.
*/
#define _MINIMUM_GHOSTS_ 1

/*! Computes the second-order finite-difference approximation to the second derivative (\b Note: not divided by the grid spacing):
    \f{equation}{
      \left(\partial^2 f\right)_i = f_{i-1} - 2f_i + f_{i+1}
    \f}
    where \f$i\f$ is the grid index along the spatial dimension of the derivative.

    \b Notes:
    + The second derivative is computed at the grid points or the cell centers.
    + Though the array D2f must include the same number of ghost points as the array f (including 0),
      the second derivative is \b not computed at these locations. Thus, array elements corresponding
      to the ghost points contain undefined values.
    + \a D2f and \a f are serialized 1D arrays containing the function and its computed derivatives on
      a multi-dimensional grid. The derivative along the specified dimension \b dir is computed by looping
      through all grid lines along \b dir.
    + Biased stencils are used at the physical boundaries and ghost point values of f are not used. So this
      function can be used when the ghost values at physical boundaries are not filled. The ghost values at
      internal (MPI) boundaries are still needed.
*/
int SecondDerivativeSecondOrderCentralNoGhosts( double*  a_D2f,   /*!< Array to hold the computed second derivative (with ghost points)
                                                                     (same size and layout as a_f) */
                                                double*  a_f,     /*!< Array containing the grid point function values whose second
                                                                     derivative is to be computed (with ghost points) */
                                                int     a_dir,    /*!< The spatial dimension along which the derivative is computed */
                                                int     a_ndims,  /*!< Number of spatial/coordinate dimensions */
                                                int*    a_dim,    /*!< Local dimensions */
                                                int     a_ghosts, /*!< Number of ghost points */
                                                int     a_nvars,  /*!< Number of vector components at each grid points */
                                                void*   a_m      /*!< MPI object of type #MPIContext */)
{
  MPIContext* mpi = (MPIContext*) a_m;
  int i, v;

  if ((!a_D2f) || (!a_f)) {
    fprintf(stderr, "Error in SecondDerivativeSecondOrder(): input arrays not allocated.\n");
    return(1);
  }
  if (a_ghosts < _MINIMUM_GHOSTS_) {
    fprintf(stderr, "Error in SecondDerivativeSecondOrderNoGhosts(): insufficient number of a_ghosts.\n");
    return(1);
  }

  /* create index and bounds for the outer loop, i.e., to loop over all 1D lines along
     dimension "a_dir"                                                                    */
  int indexC[a_ndims], index_outer[a_ndims], bounds_outer[a_ndims];
  _ArrayCopy1D_(a_dim,bounds_outer,a_ndims); bounds_outer[a_dir] =  1;

  int done = 0; _ArraySetValue_(index_outer,a_ndims,0);
  while (!done) {
    _ArrayCopy1D_(index_outer,indexC,a_ndims);
    for (i = 0; i < a_dim[a_dir]; i++) {
      int qL, qC, qR;
      indexC[a_dir] = i-1; _ArrayIndex1D_(a_ndims,a_dim,indexC,a_ghosts,qL);
      indexC[a_dir] = i  ; _ArrayIndex1D_(a_ndims,a_dim,indexC,a_ghosts,qC);
      indexC[a_dir] = i+1; _ArrayIndex1D_(a_ndims,a_dim,indexC,a_ghosts,qR);
      for (v=0; v<a_nvars; v++)  a_D2f[qC*a_nvars+v] = a_f[qL*a_nvars+v]-2*a_f[qC*a_nvars+v]+a_f[qR*a_nvars+v];
    }
    _ArrayIncrementIndex_(a_ndims,bounds_outer,index_outer,done);
  }

  if (mpi->m_ip[a_dir] == 0) {
    /* left physical boundary: overwrite the leftmost value with biased finite-difference */
    int done = 0; _ArraySetValue_(index_outer,a_ndims,0);
    while (!done) {
      _ArrayCopy1D_(index_outer,indexC,a_ndims);
      i = 0;
      int qC, qR, qR2, qR3;
      indexC[a_dir] = i  ; _ArrayIndex1D_(a_ndims,a_dim,indexC,a_ghosts,qC);
      indexC[a_dir] = i+1; _ArrayIndex1D_(a_ndims,a_dim,indexC,a_ghosts,qR);
      indexC[a_dir] = i+2; _ArrayIndex1D_(a_ndims,a_dim,indexC,a_ghosts,qR2);
      indexC[a_dir] = i+3; _ArrayIndex1D_(a_ndims,a_dim,indexC,a_ghosts,qR3);
      for (v=0; v<a_nvars; v++)  a_D2f[qC*a_nvars+v] = 2*a_f[qC*a_nvars+v]-5*a_f[qR*a_nvars+v]+4*a_f[qR2*a_nvars+v]-a_f[qR3*a_nvars+v];
      _ArrayIncrementIndex_(a_ndims,bounds_outer,index_outer,done);
    }
  }

  if (mpi->m_ip[a_dir] == (mpi->m_iproc[a_dir]-1)) {
    /* right physical boundary: overwrite the rightmost value with biased finite-difference */
    int done = 0; _ArraySetValue_(index_outer,a_ndims,0);
    while (!done) {
      _ArrayCopy1D_(index_outer,indexC,a_ndims);
      i = a_dim[a_dir]-1;
      int qL3, qL2, qL, qC;
      indexC[a_dir] = i-3; _ArrayIndex1D_(a_ndims,a_dim,indexC,a_ghosts,qL3);
      indexC[a_dir] = i-2; _ArrayIndex1D_(a_ndims,a_dim,indexC,a_ghosts,qL2);
      indexC[a_dir] = i-1; _ArrayIndex1D_(a_ndims,a_dim,indexC,a_ghosts,qL);
      indexC[a_dir] = i  ; _ArrayIndex1D_(a_ndims,a_dim,indexC,a_ghosts,qC);
      for (v=0; v<a_nvars; v++)  a_D2f[qC*a_nvars+v] = 2*a_f[qC*a_nvars+v]-5*a_f[qL*a_nvars+v]+4*a_f[qL2*a_nvars+v]-a_f[qL3*a_nvars+v];
      _ArrayIncrementIndex_(a_ndims,bounds_outer,index_outer,done);
    }
  }


  return 0;
}

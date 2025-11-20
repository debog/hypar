/*! @file ShallowWater1DModifiedSolution.c
    @author Debojyoti Ghosh
    @brief Contains the function to compute the modified solution for a balanced discretization scheme.
*/

#include <basic.h>
#include <arrayfunctions.h>
#include <physicalmodels/shallowwater1d.h>
#include <mpivars.h>
#include <hypar.h>

/*! Compute the modified solution for the upwinding step in a balanced conservative
    finite-difference algorithm for the 1D shallow water equations.
    \n\n
  Refer to:
  + Xing, Y., Shu, C.-W., "High order finite difference WENO schemes with the
    exact conservation property for the shallow water equations", Journal of
    Computational Physics, 208, 2005, pp. 206-227.
    http://dx.doi.org/10.1016/j.jcp.2005.02.006
*/
int ShallowWater1DModifiedSolution(
                            double  *a_uC, /*!< The modified solution (same array size and layout as a_u) */
                            double  *a_u,  /*!< The solution (conserved variables) */
                            int     a_d,   /*!< Spatial dimension (unused since this is a 1D system) */
                            void    *a_s,  /*!< Solver object of type #HyPar */
                            void    *a_m,  /*!< MPI object of type #MPIVariables */
                            double  waqt /*!< Current solution time */
                           )
{
  HyPar           *solver = (HyPar*) a_s;
  ShallowWater1D  *param  = (ShallowWater1D*) solver->m_physics;

  int ghosts = solver->m_ghosts;
  int     *dim    = solver->m_dim_local;
  int     ndims   = solver->m_ndims;
  int     index[ndims], bounds[ndims], offset[ndims];

  /* set bounds for array index to include ghost points */
  _ArrayCopy1D_(dim,bounds,ndims);
  int i; for (i=0; i<ndims; i++) bounds[i] += 2*ghosts;

  /* set offset such that index is compatible with ghost point arrangement */
  _ArraySetValue_(offset,ndims,-ghosts);

  int done = 0; _ArraySetValue_(index,ndims,0);
  while (!done) {
    int p; _ArrayIndex1DWO_(ndims,dim,index,offset,ghosts,p);
    double h,v;
    _ShallowWater1DGetFlowVar_((a_u+_MODEL_NVARS_*p),h,v);
    a_uC[_MODEL_NVARS_*p+0] = h + param->m_b[p];
    a_uC[_MODEL_NVARS_*p+1] = h * v;
    _ArrayIncrementIndex_(ndims,bounds,index,done);
  }
  return(0);
}

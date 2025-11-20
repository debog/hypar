/*! @file ShallowWater1DSourceUpwind.c
    @author Debojyoti Ghosh
    @brief Contains functions to compute the "upwind" source term at an interface (for a balanced finite-difference discretization of the 1D shallow water equations).
*/

#include <stdlib.h>
#include <math.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mathfunctions.h>
#include <physicalmodels/shallowwater1d.h>
#include <hypar.h>

/*! Compute the "upwind" source term in the balanced formulation introduced in the
    reference below. The "upwind" state is just the arithmetic average of the left
    and right states.
    + Xing, Y., Shu, C.-W., "High order finite difference WENO schemes with the
      exact conservation property for the shallow water equations", Journal of
      Computational Physics, 208, 2005, pp. 206-227.
      http://dx.doi.org/10.1016/j.jcp.2005.02.006
*/
int ShallowWater1DSourceUpwindLLF(
                                  double  *a_fI, /*!< Computed interface source term ("upwinded") */
                                  double  *a_fL, /*!< Left-biased interface source term */
                                  double  *a_fR, /*!< Right-biased interface source term */
                                  double  *a_u,  /*!< Solution (conserved variables) */
                                  int     a_dir, /*!< Spatial dimension (unused since this is a 1D case) */
                                  void    *a_s,  /*!< Solver object of type #HyPar */
                                  double  a_t   /*!< Current solution time */
                                 )
{
  HyPar     *solver = (HyPar*)    a_s;
  int       done,k;

  int ndims = solver->m_ndims;
  int *dim  = solver->m_dim_local;

  int index_outer[ndims], index_inter[ndims], bounds_outer[ndims], bounds_inter[ndims];
  _ArrayCopy1D_(dim,bounds_outer,ndims); bounds_outer[a_dir] =  1;
  _ArrayCopy1D_(dim,bounds_inter,ndims); bounds_inter[a_dir] += 1;

  done = 0; _ArraySetValue_(index_outer,ndims,0);
  while (!done) {
    _ArrayCopy1D_(index_outer,index_inter,ndims);
    for (index_inter[a_dir] = 0; index_inter[a_dir] < bounds_inter[a_dir]; index_inter[a_dir]++) {
      int p;  _ArrayIndex1D_(ndims,bounds_inter,index_inter,0,p);
      for (k = 0; k < _MODEL_NVARS_; k++)
        (a_fI+_MODEL_NVARS_*p)[k] = 0.5 * ((a_fL+_MODEL_NVARS_*p)[k] + (a_fR+_MODEL_NVARS_*p)[k]);
    }
    _ArrayIncrementIndex_(ndims,bounds_outer,index_outer,done);
  }

  return(0);
}

/*! Compute the "upwind" source term in the balanced formulation introduced in the
    reference below. The "upwind" state is just the arithmetic average of the left
    and right states.
    + Xing, Y., Shu, C.-W., "High order finite difference WENO schemes with the
      exact conservation property for the shallow water equations", Journal of
      Computational Physics, 208, 2005, pp. 206-227.
      http://dx.doi.org/10.1016/j.jcp.2005.02.006
*/
int ShallowWater1DSourceUpwindRoe(
                                  double  *a_fI, /*!< Computed interface source term ("upwinded") */
                                  double  *a_fL, /*!< Left-biased interface source term */
                                  double  *a_fR, /*!< Right-biased interface source term */
                                  double  *a_u,  /*!< Solution (conserved variables) */
                                  int     a_dir, /*!< Spatial dimension (unused since this is a 1D case) */
                                  void    *a_s,  /*!< Solver object of type #HyPar */
                                  double  a_t   /*!< Current solution time */
                                 )
{
  HyPar     *solver = (HyPar*)    a_s;
  int       done,k;

  int ndims = solver->m_ndims;
  int *dim  = solver->m_dim_local;

  int index_outer[ndims], index_inter[ndims], bounds_outer[ndims], bounds_inter[ndims];
  _ArrayCopy1D_(dim,bounds_outer,ndims); bounds_outer[a_dir] =  1;
  _ArrayCopy1D_(dim,bounds_inter,ndims); bounds_inter[a_dir] += 1;

  done = 0; _ArraySetValue_(index_outer,ndims,0);
  while (!done) {
    _ArrayCopy1D_(index_outer,index_inter,ndims);
    for (index_inter[a_dir] = 0; index_inter[a_dir] < bounds_inter[a_dir]; index_inter[a_dir]++) {
      int p;  _ArrayIndex1D_(ndims,bounds_inter,index_inter,0,p);
      for (k = 0; k < _MODEL_NVARS_; k++)
        (a_fI+_MODEL_NVARS_*p)[k] = 0.5 * ((a_fL+_MODEL_NVARS_*p)[k] + (a_fR+_MODEL_NVARS_*p)[k]);
    }
    _ArrayIncrementIndex_(ndims,bounds_outer,index_outer,done);
  }

  return(0);
}

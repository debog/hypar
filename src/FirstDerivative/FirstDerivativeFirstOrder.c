/*! @file FirstDerivativeFirstOrder.c
    @author Debojyoti Ghosh
    @brief First order approximation to the first derivative
*/

#include <stdio.h>
#include <stdlib.h>
#include <basic.h>
#include <mathfunctions.h>
#include <arrayfunctions.h>
#include <firstderivative.h>

#include <mpivars.h>
#include <hypar.h>
typedef MPIVariables  MPIContext;
typedef HyPar         SolverContext;

#ifdef with_omp
#include <omp.h>
#endif

/*! Computes the first-order finite-difference approximation to the first derivative
    (\b Note: not divided by the grid spacing):
    \f{equation}{
      \left(\partial f\right)_i = \left\{ \begin{array}{ll} f_{i+1} - f_i  & {\rm bias} = 1 \\ f_i - f_{i-1}  & {\rm bias} = -1 \end{array}\right.
    \f}
    where \f$i\f$ is the grid index along the spatial dimension of the derivative.
    \n\n
    Notes:
    + The first derivative is computed at the grid points or the cell centers.
    + The first derivative is computed at the ghost points too. Thus, biased schemes are used
      at and near the boundaries.
    + \b Df and \b f are 1D arrays containing the function and its computed derivatives on a multi-
      dimensional grid. The derivative along the specified dimension \b dir is computed by looping
      through all grid lines along \b dir.
*/
int FirstDerivativeFirstOrder(
                                double  *a_Df,  /*!< Array to hold the computed first derivative (with ghost points) */
                                double  *a_f,   /*!< Array containing the grid point function values whose first
                                                   derivative is to be computed (with ghost points) */
                                int     a_dir,  /*!< The spatial dimension along which the derivative is computed */
                                int     a_bias, /*!< Forward or backward differencing for non-central
                                                   finite-difference schemes (-1: backward, 1: forward)*/
                                void    *a_s,   /*!< Solver object of type #SolverContext */
                                void *a_m    /*!< MPI object of type #MPIContext */
                             )
{
  SolverContext *solver = (SolverContext*) a_s;
  int           i, j, v;

  int ghosts = solver->m_ghosts;
  int ndims  = solver->m_ndims;
  int nvars  = solver->m_nvars;
  int *dim   = solver->m_dim_local;


  if ((!a_Df) || (!a_f)) {
    fprintf(stderr, "Error in FirstDerivativeSecondOrder(): input arrays not allocated.\n");
    return(1);
  }

  /* create index and bounds for the outer loop, i.e., to loop over all 1D lines along
     dimension "a_dir"                                                                    */
  int indexC[ndims], index_outer[ndims], bounds_outer[ndims];
  _ArrayCopy1D_(dim,bounds_outer,ndims); bounds_outer[a_dir] =  1;
  int N_outer; _ArrayProduct1D_(bounds_outer,ndims,N_outer);

#pragma omp parallel for schedule(auto) default(shared) private(i,j,v,index_outer,indexC)
  for (j=0; j<N_outer; j++) {
    _ArrayIndexnD_(ndims,j,bounds_outer,index_outer,0);
    _ArrayCopy1D_(index_outer,indexC,ndims);
    /* left boundary */
    for (i = -ghosts; i < -ghosts+1; i++) {
      int qC, qR;
      indexC[a_dir] = i  ; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qC );
      indexC[a_dir] = i+1; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qR );
      for (v=0; v<nvars; v++)  a_Df[qC*nvars+v] = a_f[qR*nvars+v]-a_f[qC*nvars+v];
    }
    /* interior */
    for (i = -ghosts+1; i < dim[a_dir]+ghosts-1; i++) {
      int qC, qL, qR;
      indexC[a_dir] = i  ; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qC );
      indexC[a_dir] = i-1; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qL);
      indexC[a_dir] = i+1; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qR);
      for (v=0; v<nvars; v++)  a_Df[qC*nvars+v] = max(a_bias,0)*a_f[qR*nvars+v]-a_bias*a_f[qC*nvars+v]+min(a_bias,0)*a_f[qL*nvars+v];
    }
    /* right boundary */
    for (i = dim[a_dir]+ghosts-1; i < dim[a_dir]+ghosts; i++) {
      int qL, qC;
      indexC[a_dir] = i-1; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qL );
      indexC[a_dir] = i  ; _ArrayIndex1D_(ndims,dim,indexC,ghosts,qC );
      for (v=0; v<nvars; v++)  a_Df[qC*nvars+v] = a_f[qC*nvars+v]-a_f[qL*nvars+v];
    }
  }

  return(0);
}

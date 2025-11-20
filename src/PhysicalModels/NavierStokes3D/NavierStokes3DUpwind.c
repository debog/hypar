/*! @file NavierStokes3DUpwind.c
    @author Debojyoti Ghosh
    @brief Contains functions to compute the upwind flux at grid interfaces for the 3D Navier Stokes equations.
*/
#include <stdlib.h>
#include <math.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mathfunctions.h>
#include <matmult_native.h>
#include <physicalmodels/navierstokes3d.h>
#include <hypar.h>

#if defined(CPU_STAT)
#include <time.h>
#endif

static const int dummy = 1;

/*! Roe's upwinding scheme.
    \f{equation}{
      {\bf f}_{j+1/2} = \frac{1}{2}\left[ {\bf f}_{j+1/2}^L + {\bf f}_{j+1/2}^R
                         - \left| A\left({\bf u}_{j+1/2}^L,{\bf u}_{j+1/2}^R\right) \right|
                           \left( {\bf u}_{j+1/2}^R - {\bf u}_{j+1/2}^L  \right)\right]
    \f}
    + Roe, P. L., “Approximate Riemann solvers, parameter vectors, and difference schemes,” Journal of
    Computational Physics, Vol. 43, No. 2, 1981, pp. 357–372, http://dx.doi.org/10.1016/0021-9991(81)90128-5.

    This upwinding scheme is modified for the balanced discretization of the 3D Navier Stokes equations when
    there is a non-zero gravitational force. See the reference below. For flows without any gravitational forces,
    it reduces to its original form.
    + Ghosh, D., Constantinescu, E.M., Well-Balanced Formulation of Gravitational Source
      Terms for Conservative Finite-Difference Atmospheric Flow Solvers, AIAA Paper 2015-2889,
      7th AIAA Atmospheric and Space Environments Conference, June 22-26, 2015, Dallas, TX,
      http://dx.doi.org/10.2514/6.2015-2889
    + Ghosh, D., Constantinescu, E.M., A Well-Balanced, Conservative Finite-Difference Algorithm
      for Atmospheric Flows, 54 (4), 2016, pp. 1370-1385, http://dx.doi.org/10.2514/1.J054580

*/
int NavierStokes3DUpwindRoe(
                            double  *a_fI, /*!< Computed upwind interface flux */
                            double  *a_fL, /*!< Left-biased reconstructed interface flux */
                            double  *a_fR, /*!< Right-biased reconstructed interface flux */
                            double  *a_uL, /*!< Left-biased reconstructed interface solution */
                            double  *a_uR, /*!< Right-biased reconstructed interface solution */
                            double  *a_u,  /*!< Cell-centered solution */
                            int     a_dir, /*!< Spatial dimension (x, y, or z) */
                            void    *a_s,  /*!< Solver object of type #HyPar */
                            double  a_t   /*!< Current solution time */
                           )
{
  HyPar           *solver = (HyPar*)    a_s;
  NavierStokes3D  *param  = (NavierStokes3D*)  solver->m_physics;
  int             done;

  int *dim  = solver->m_dim_local;

  int bounds_outer[_MODEL_NDIMS_], bounds_inter[_MODEL_NDIMS_];
  _ArrayCopy1D3_(dim,bounds_outer,_MODEL_NDIMS_); bounds_outer[a_dir] =  1;
  _ArrayCopy1D3_(dim,bounds_inter,_MODEL_NDIMS_); bounds_inter[a_dir] += 1;
  static double R[_MODEL_NVARS_*_MODEL_NVARS_], D[_MODEL_NVARS_*_MODEL_NVARS_],
                L[_MODEL_NVARS_*_MODEL_NVARS_], DL[_MODEL_NVARS_*_MODEL_NVARS_],
                modA[_MODEL_NVARS_*_MODEL_NVARS_];

  done = 0; int index_outer[3] = {0,0,0}, index_inter[3];

#if defined(CPU_STAT)
  clock_t startEvent, stopEvent;
  startEvent = clock();
#endif

  while (!done) {
    _ArrayCopy1D3_(index_outer,index_inter,_MODEL_NDIMS_);
    for (index_inter[a_dir] = 0; index_inter[a_dir] < bounds_inter[a_dir]; index_inter[a_dir]++) {
      int p; _ArrayIndex1D3_(_MODEL_NDIMS_,bounds_inter,index_inter,0,p);
      int indexL[_MODEL_NDIMS_]; _ArrayCopy1D_(index_inter,indexL,_MODEL_NDIMS_); indexL[a_dir]--;
      int indexR[_MODEL_NDIMS_]; _ArrayCopy1D_(index_inter,indexR,_MODEL_NDIMS_);
      int pL; _ArrayIndex1D_(_MODEL_NDIMS_,dim,indexL,solver->m_ghosts,pL);
      int pR; _ArrayIndex1D_(_MODEL_NDIMS_,dim,indexR,solver->m_ghosts,pR);
      double udiff[_MODEL_NVARS_], uavg[_MODEL_NVARS_],udiss[_MODEL_NVARS_];

      /* Roe'a_s upwinding scheme */

      udiff[0] = 0.5 * (a_uR[_MODEL_NVARS_*p+0] - a_uL[_MODEL_NVARS_*p+0]);
      udiff[1] = 0.5 * (a_uR[_MODEL_NVARS_*p+1] - a_uL[_MODEL_NVARS_*p+1]);
      udiff[2] = 0.5 * (a_uR[_MODEL_NVARS_*p+2] - a_uL[_MODEL_NVARS_*p+2]);
      udiff[3] = 0.5 * (a_uR[_MODEL_NVARS_*p+3] - a_uL[_MODEL_NVARS_*p+3]);
      udiff[4] = 0.5 * (a_uR[_MODEL_NVARS_*p+4] - a_uL[_MODEL_NVARS_*p+4]);

      _NavierStokes3DRoeAverage_        (uavg,_NavierStokes3D_stride_,(a_u+_MODEL_NVARS_*pL),(a_u+_MODEL_NVARS_*pR),param->m_gamma);
      _NavierStokes3DEigenvalues_       (uavg,dummy,D,param->m_gamma,a_dir);
      _NavierStokes3DLeftEigenvectors_  (uavg,dummy,L,param->m_gamma,a_dir);
      _NavierStokes3DRightEigenvectors_ (uavg,dummy,R,param->m_gamma,a_dir);

      /* Harten'a_s Entropy Fix - Page 362 of Leveque */
      int k;
      double delta = 0.000001, delta2 = delta*delta;
      k=0;  D[k] = (absolute(D[k]) < delta ? (D[k]*D[k]+delta2)/(2*delta) : absolute(D[k]) );
      k=6;  D[k] = (absolute(D[k]) < delta ? (D[k]*D[k]+delta2)/(2*delta) : absolute(D[k]) );
      k=12; D[k] = (absolute(D[k]) < delta ? (D[k]*D[k]+delta2)/(2*delta) : absolute(D[k]) );
      k=18; D[k] = (absolute(D[k]) < delta ? (D[k]*D[k]+delta2)/(2*delta) : absolute(D[k]) );
      k=24; D[k] = (absolute(D[k]) < delta ? (D[k]*D[k]+delta2)/(2*delta) : absolute(D[k]) );

      MatMult5(_MODEL_NVARS_,DL,D,L);
      MatMult5(_MODEL_NVARS_,modA,R,DL);
      MatVecMult5(_MODEL_NVARS_,udiss,modA,udiff);

      a_fI[_MODEL_NVARS_*p+0] = 0.5 * (a_fL[_MODEL_NVARS_*p+0]+a_fR[_MODEL_NVARS_*p+0]) - udiss[0];
      a_fI[_MODEL_NVARS_*p+1] = 0.5 * (a_fL[_MODEL_NVARS_*p+1]+a_fR[_MODEL_NVARS_*p+1]) - udiss[1];
      a_fI[_MODEL_NVARS_*p+2] = 0.5 * (a_fL[_MODEL_NVARS_*p+2]+a_fR[_MODEL_NVARS_*p+2]) - udiss[2];
      a_fI[_MODEL_NVARS_*p+3] = 0.5 * (a_fL[_MODEL_NVARS_*p+3]+a_fR[_MODEL_NVARS_*p+3]) - udiss[3];
      a_fI[_MODEL_NVARS_*p+4] = 0.5 * (a_fL[_MODEL_NVARS_*p+4]+a_fR[_MODEL_NVARS_*p+4]) - udiss[4];
    }
    _ArrayIncrementIndex_(_MODEL_NDIMS_,bounds_outer,index_outer,done);
  }

#if defined(CPU_STAT)
  stopEvent = clock();
  printf("%-50s CPU time (secs) = %.6f a_dir = %d\n",
          "NavierStokes3DUpwindRoe", (double)(stopEvent-startEvent)/CLOCKS_PER_SEC, a_dir);
#endif

  return(0);
}

/*! Characteristic-based Roe-fixed upwinding scheme.
    \f{align}{
      \alpha_{j+1/2}^{k,L} &= \sum_{k=1}^3 {\bf l}_{j+1/2}^k \cdot {\bf f}_{j+1/2}^{k,L}, \\
      \alpha_{j+1/2}^{k,R} &= \sum_{k=1}^3 {\bf l}_{j+1/2}^k \cdot {\bf f}_{j+1/2}^{k,R}, \\
      v_{j+1/2}^{k,L} &= \sum_{k=1}^3 {\bf l}_{j+1/2}^k \cdot {\bf u}_{j+1/2}^{k,L}, \\
      v_{j+1/2}^{k,R} &= \sum_{k=1}^3 {\bf l}_{j+1/2}^k \cdot {\bf u}_{j+1/2}^{k,R}, \\
      \alpha_{j+1/2}^k &= \left\{ \begin{array}{cc} \alpha_{j+1/2}^{k,L} & {\rm if}\ \lambda_{j,j+1/2,j+1} > 0 \\ \alpha_{j+1/2}^{k,R} & {\rm if}\ \lambda_{j,j+1/2,j+1} < 0 \\ \frac{1}{2}\left[ \alpha_{j+1/2}^{k,L} + \alpha_{j+1/2}^{k,R} - \left(\max_{\left[j,j+1\right]} \lambda\right) \left( v_{j+1/2}^{k,R} - v_{j+1/2}^{k,L} \right) \right] & {\rm otherwise} \end{array}\right., \\
      {\bf f}_{j+1/2} &= \sum_{k=1}^3 \alpha_{j+1/2}^k {\bf r}_{j+1/2}^k
    \f}
    where \f${\bf l}\f$, \f${\bf r}\f$, and \f$\lambda\f$ are the left-eigenvectors, right-eigenvectors and eigenvalues. The subscripts denote the grid locations.
    + C.-W. Shu, and S. Osher, "Efficient implementation of essentially non-oscillatory schemes, II", J. Comput. Phys., 83 (1989), pp. 32–78, http://dx.doi.org/10.1016/0021-9991(89)90222-2.

    Note that this upwinding scheme cannot be used for solving flows with non-zero gravitational forces.
*/
int NavierStokes3DUpwindRF(
                            double  *a_fI, /*!< Computed upwind interface flux */
                            double  *a_fL, /*!< Left-biased reconstructed interface flux */
                            double  *a_fR, /*!< Right-biased reconstructed interface flux */
                            double  *a_uL, /*!< Left-biased reconstructed interface solution */
                            double  *a_uR, /*!< Right-biased reconstructed interface solution */
                            double  *a_u,  /*!< Cell-centered solution */
                            int     a_dir, /*!< Spatial dimension (x, y, or z) */
                            void    *a_s,  /*!< Solver object of type #HyPar */
                            double  a_t   /*!< Current solution time */
                          )
{
  HyPar           *solver = (HyPar*)    a_s;
  NavierStokes3D  *param  = (NavierStokes3D*)  solver->m_physics;
  int             done,k;

  int *dim  = solver->m_dim_local;

  int bounds_outer[_MODEL_NDIMS_], bounds_inter[_MODEL_NDIMS_];
  _ArrayCopy1D3_(dim,bounds_outer,_MODEL_NDIMS_); bounds_outer[a_dir] =  1;
  _ArrayCopy1D3_(dim,bounds_inter,_MODEL_NDIMS_); bounds_inter[a_dir] += 1;
  static double R[_MODEL_NVARS_*_MODEL_NVARS_], D[_MODEL_NVARS_*_MODEL_NVARS_], L[_MODEL_NVARS_*_MODEL_NVARS_];

  done = 0; int index_outer[3] = {0,0,0}, index_inter[3];
  while (!done) {
    _ArrayCopy1D3_(index_outer,index_inter,_MODEL_NDIMS_);
    for (index_inter[a_dir] = 0; index_inter[a_dir] < bounds_inter[a_dir]; index_inter[a_dir]++) {
      int p; _ArrayIndex1D3_(_MODEL_NDIMS_,bounds_inter,index_inter,0,p);
      double uavg[_MODEL_NVARS_], fcL[_MODEL_NVARS_], fcR[_MODEL_NVARS_],
             ucL[_MODEL_NVARS_], ucR[_MODEL_NVARS_], fc[_MODEL_NVARS_];

      /* Roe-Fixed upwinding scheme */

      _NavierStokes3DRoeAverage_(uavg,_NavierStokes3D_stride_,(a_uL+_MODEL_NVARS_*p),(a_uR+_MODEL_NVARS_*p),param->m_gamma);

      _NavierStokes3DLeftEigenvectors_(uavg,dummy,L,param->m_gamma,a_dir);
      _NavierStokes3DRightEigenvectors_(uavg,dummy,R,param->m_gamma,a_dir);

      /* calculate characteristic fluxes and variables */
      MatVecMult5(_MODEL_NVARS_,ucL,L,(a_uL+_MODEL_NVARS_*p));
      MatVecMult5(_MODEL_NVARS_,ucR,L,(a_uR+_MODEL_NVARS_*p));
      MatVecMult5(_MODEL_NVARS_,fcL,L,(a_fL+_MODEL_NVARS_*p));
      MatVecMult5(_MODEL_NVARS_,fcR,L,(a_fR+_MODEL_NVARS_*p));

      double eigL[_MODEL_NVARS_],eigC[_MODEL_NVARS_],eigR[_MODEL_NVARS_];
      _NavierStokes3DEigenvalues_((a_uL+_MODEL_NVARS_*p),_NavierStokes3D_stride_,D,param->m_gamma,a_dir);
      eigL[0] = D[0];
      eigL[1] = D[6];
      eigL[2] = D[12];
      eigL[3] = D[18];
      eigL[4] = D[24];
      _NavierStokes3DEigenvalues_((a_uR+_MODEL_NVARS_*p),_NavierStokes3D_stride_,D,param->m_gamma,a_dir);
      eigR[0] = D[0];
      eigR[1] = D[6];
      eigR[2] = D[12];
      eigR[3] = D[18];
      eigR[4] = D[24];
      _NavierStokes3DEigenvalues_(uavg,dummy,D,param->m_gamma,a_dir);
      eigC[0] = D[0];
      eigC[1] = D[6];
      eigC[2] = D[12];
      eigC[3] = D[18];
      eigC[4] = D[24];

      for (k = 0; k < _MODEL_NVARS_; k++) {
        if ((eigL[k] > 0) && (eigC[k] > 0) && (eigR[k] > 0))      fc[k] = fcL[k];
        else if ((eigL[k] < 0) && (eigC[k] < 0) && (eigR[k] < 0)) fc[k] = fcR[k];
        else {
          double alpha = max3(absolute(eigL[k]),absolute(eigC[k]),absolute(eigR[k]));
          fc[k] = 0.5 * (fcL[k] + fcR[k] + alpha * (ucL[k]-ucR[k]));
        }
      }

      /* calculate the interface flux from the characteristic flux */
      MatVecMult5(_MODEL_NVARS_,(a_fI+_MODEL_NVARS_*p),R,fc);
    }
    _ArrayIncrementIndex_(_MODEL_NDIMS_,bounds_outer,index_outer,done);
  }

  return(0);
}

/*! Characteristic-based local Lax-Friedrich upwinding scheme.
    \f{align}{
      \alpha_{j+1/2}^{k,L} &= \sum_{k=1}^3 {\bf l}_{j+1/2}^k \cdot {\bf f}_{j+1/2}^{k,L}, \\
      \alpha_{j+1/2}^{k,R} &= \sum_{k=1}^3 {\bf l}_{j+1/2}^k \cdot {\bf f}_{j+1/2}^{k,R}, \\
      v_{j+1/2}^{k,L} &= \sum_{k=1}^3 {\bf l}_{j+1/2}^k \cdot {\bf u}_{j+1/2}^{k,L}, \\
      v_{j+1/2}^{k,R} &= \sum_{k=1}^3 {\bf l}_{j+1/2}^k \cdot {\bf u}_{j+1/2}^{k,R}, \\
      \alpha_{j+1/2}^k &= \frac{1}{2}\left[ \alpha_{j+1/2}^{k,L} + \alpha_{j+1/2}^{k,R} - \left(\max_{\left[j,j+1\right]} \lambda\right) \left( v_{j+1/2}^{k,R} - v_{j+1/2}^{k,L} \right) \right], \\
      {\bf f}_{j+1/2} &= \sum_{k=1}^3 \alpha_{j+1/2}^k {\bf r}_{j+1/2}^k
    \f}
    where \f${\bf l}\f$, \f${\bf r}\f$, and \f$\lambda\f$ are the left-eigenvectors, right-eigenvectors and eigenvalues. The subscripts denote the grid locations.
    + C.-W. Shu, and S. Osher, "Efficient implementation of essentially non-oscillatory schemes, II", J. Comput. Phys., 83 (1989), pp. 32–78, http://dx.doi.org/10.1016/0021-9991(89)90222-2.

    This upwinding scheme is modified for the balanced discretization of the 3D Navier Stokes equations when
    there is a non-zero gravitational force. See the reference below. For flows without any gravitational forces,
    it reduces to its original form.
    + Ghosh, D., Constantinescu, E.M., Well-Balanced Formulation of Gravitational Source
      Terms for Conservative Finite-Difference Atmospheric Flow Solvers, AIAA Paper 2015-2889,
      7th AIAA Atmospheric and Space Environments Conference, June 22-26, 2015, Dallas, TX,
      http://dx.doi.org/10.2514/6.2015-2889
    + Ghosh, D., Constantinescu, E.M., A Well-Balanced, Conservative Finite-Difference Algorithm
      for Atmospheric Flows, 54 (4), 2016, pp. 1370-1385, http://dx.doi.org/10.2514/1.J054580
*/
int NavierStokes3DUpwindLLF(
                            double  *a_fI, /*!< Computed upwind interface flux */
                            double  *a_fL, /*!< Left-biased reconstructed interface flux */
                            double  *a_fR, /*!< Right-biased reconstructed interface flux */
                            double  *a_uL, /*!< Left-biased reconstructed interface solution */
                            double  *a_uR, /*!< Right-biased reconstructed interface solution */
                            double  *a_u,  /*!< Cell-centered solution */
                            int     a_dir, /*!< Spatial dimension (x, y, or z) */
                            void    *a_s,  /*!< Solver object of type #HyPar */
                            double  a_t   /*!< Current solution time */
                           )
{
  HyPar           *solver = (HyPar*)    a_s;
  NavierStokes3D  *param  = (NavierStokes3D*)  solver->m_physics;
  int             done;

  int *dim  = solver->m_dim_local;

  int bounds_outer[_MODEL_NDIMS_], bounds_inter[_MODEL_NDIMS_];
  _ArrayCopy1D3_(dim,bounds_outer,_MODEL_NDIMS_); bounds_outer[a_dir] =  1;
  _ArrayCopy1D3_(dim,bounds_inter,_MODEL_NDIMS_); bounds_inter[a_dir] += 1;
  static double R[_MODEL_NVARS_*_MODEL_NVARS_], D[_MODEL_NVARS_*_MODEL_NVARS_], L[_MODEL_NVARS_*_MODEL_NVARS_];

  done = 0; int index_outer[3] = {0,0,0}, index_inter[3];
  while (!done) {
    _ArrayCopy1D3_(index_outer,index_inter,_MODEL_NDIMS_);
    for (index_inter[a_dir] = 0; index_inter[a_dir] < bounds_inter[a_dir]; index_inter[a_dir]++) {
      int p; _ArrayIndex1D3_(_MODEL_NDIMS_,bounds_inter,index_inter,0,p);
      double uavg[_MODEL_NVARS_], fcL[_MODEL_NVARS_], fcR[_MODEL_NVARS_],
             ucL[_MODEL_NVARS_], ucR[_MODEL_NVARS_], fc[_MODEL_NVARS_];

      /* Local Lax-Friedrich upwinding scheme */

      _NavierStokes3DRoeAverage_(uavg,_NavierStokes3D_stride_,(a_uL+_MODEL_NVARS_*p),(a_uR+_MODEL_NVARS_*p),param->m_gamma);

      _NavierStokes3DLeftEigenvectors_(uavg,dummy,L,param->m_gamma,a_dir);
      _NavierStokes3DRightEigenvectors_(uavg,dummy,R,param->m_gamma,a_dir);

      /* calculate characteristic fluxes and variables */
      MatVecMult5(_MODEL_NVARS_,ucL,L,(a_uL+_MODEL_NVARS_*p));
      MatVecMult5(_MODEL_NVARS_,ucR,L,(a_uR+_MODEL_NVARS_*p));
      MatVecMult5(_MODEL_NVARS_,fcL,L,(a_fL+_MODEL_NVARS_*p));
      MatVecMult5(_MODEL_NVARS_,fcR,L,(a_fR+_MODEL_NVARS_*p));

      double eigL[_MODEL_NVARS_],eigC[_MODEL_NVARS_],eigR[_MODEL_NVARS_];
      _NavierStokes3DEigenvalues_((a_uL+_MODEL_NVARS_*p),_NavierStokes3D_stride_,D,param->m_gamma,a_dir);
      eigL[0] = D[0];
      eigL[1] = D[6];
      eigL[2] = D[12];
      eigL[3] = D[18];
      eigL[4] = D[24];
      _NavierStokes3DEigenvalues_((a_uR+_MODEL_NVARS_*p),_NavierStokes3D_stride_,D,param->m_gamma,a_dir);
      eigR[0] = D[0];
      eigR[1] = D[6];
      eigR[2] = D[12];
      eigR[3] = D[18];
      eigR[4] = D[24];
      _NavierStokes3DEigenvalues_(uavg,dummy,D,param->m_gamma,a_dir);
      eigC[0] = D[0];
      eigC[1] = D[6];
      eigC[2] = D[12];
      eigC[3] = D[18];
      eigC[4] = D[24];

      double alpha;
      alpha = max3(absolute(eigL[0]),absolute(eigC[0]),absolute(eigR[0]));
      fc[0] = 0.5 * (fcL[0] + fcR[0] + alpha * (ucL[0]-ucR[0]));
      alpha = max3(absolute(eigL[1]),absolute(eigC[1]),absolute(eigR[1]));
      fc[1] = 0.5 * (fcL[1] + fcR[1] + alpha * (ucL[1]-ucR[1]));
      alpha = max3(absolute(eigL[2]),absolute(eigC[2]),absolute(eigR[2]));
      fc[2] = 0.5 * (fcL[2] + fcR[2] + alpha * (ucL[2]-ucR[2]));
      alpha = max3(absolute(eigL[3]),absolute(eigC[3]),absolute(eigR[3]));
      fc[3] = 0.5 * (fcL[3] + fcR[3] + alpha * (ucL[3]-ucR[3]));
      alpha = max3(absolute(eigL[4]),absolute(eigC[4]),absolute(eigR[4]));
      fc[4] = 0.5 * (fcL[4] + fcR[4] + alpha * (ucL[4]-ucR[4]));

      /* calculate the interface flux from the characteristic flux */
      MatVecMult5(_MODEL_NVARS_,(a_fI+_MODEL_NVARS_*p),R,fc);
    }
    _ArrayIncrementIndex_(_MODEL_NDIMS_,bounds_outer,index_outer,done);
  }

  return(0);
}

/*! Rusanov's upwinding scheme.
    \f{equation}{
      {\bf f}_{j+1/2} = \frac{1}{2}\left[ {\bf f}_{j+1/2}^L + {\bf f}_{j+1/2}^R
                         - \max_{j,j+1} \nu_j \left( {\bf u}_{j+1/2}^R - {\bf u}_{j+1/2}^L  \right)\right]
    \f}
    where \f$\nu = c + \left|u\right|\f$.
    + Rusanov, V. V., "The calculation of the interaction of non-stationary shock waves and obstacles," USSR
    Computational Mathematics and Mathematical Physics, Vol. 1, No. 2, 1962, pp. 304–320

    This upwinding scheme is modified for the balanced discretization of the 3D Navier Stokes equations when
    there is a non-zero gravitational force. See the reference below. For flows without any gravitational forces,
    it reduces to its original form.
    + Ghosh, D., Constantinescu, E.M., Well-Balanced Formulation of Gravitational Source
      Terms for Conservative Finite-Difference Atmospheric Flow Solvers, AIAA Paper 2015-2889,
      7th AIAA Atmospheric and Space Environments Conference, June 22-26, 2015, Dallas, TX,
      http://dx.doi.org/10.2514/6.2015-2889
    + Ghosh, D., Constantinescu, E.M., A Well-Balanced, Conservative Finite-Difference Algorithm
      for Atmospheric Flows, 54 (4), 2016, pp. 1370-1385, http://dx.doi.org/10.2514/1.J054580

*/
int NavierStokes3DUpwindRusanov(
                                double  *a_fI, /*!< Computed upwind interface flux */
                                double  *a_fL, /*!< Left-biased reconstructed interface flux */
                                double  *a_fR, /*!< Right-biased reconstructed interface flux */
                                double  *a_uL, /*!< Left-biased reconstructed interface solution */
                                double  *a_uR, /*!< Right-biased reconstructed interface solution */
                                double  *a_u,  /*!< Cell-centered solution */
                                int     a_dir, /*!< Spatial dimension (x,y, or z) */
                                void    *a_s,  /*!< Solver object of type #HyPar */
                                double  a_t   /*!< Current solution time */
                               )
{
  HyPar           *solver = (HyPar*)          a_s;
  NavierStokes3D  *param  = (NavierStokes3D*) solver->m_physics;
  int             *dim    = solver->m_dim_local, done;

  int bounds_outer[_MODEL_NDIMS_], bounds_inter[_MODEL_NDIMS_];
  bounds_outer[0] = dim[0]; bounds_outer[1] = dim[1]; bounds_outer[2] = dim[2]; bounds_outer[a_dir] = 1;
  bounds_inter[0] = dim[0]; bounds_inter[1] = dim[1]; bounds_inter[2] = dim[2]; bounds_inter[a_dir]++;

  static double udiff[_MODEL_NVARS_],uavg[_MODEL_NVARS_];

  static int index_outer[_MODEL_NDIMS_], index_inter[_MODEL_NDIMS_],
             indexL[_MODEL_NDIMS_], indexR[_MODEL_NDIMS_];

  done = 0; _ArraySetValue_(index_outer,_MODEL_NDIMS_,0);
  while (!done) {
    _ArrayCopy1D_(index_outer,index_inter,_MODEL_NDIMS_);
    for (index_inter[a_dir] = 0; index_inter[a_dir] < bounds_inter[a_dir]; index_inter[a_dir]++) {

      int p; _ArrayIndex1D_(_MODEL_NDIMS_,bounds_inter,index_inter,0,p);
      _ArrayCopy1D_(index_inter,indexL,_MODEL_NDIMS_); indexL[a_dir]--;
      _ArrayCopy1D_(index_inter,indexR,_MODEL_NDIMS_);
      int pL; _ArrayIndex1D_(_MODEL_NDIMS_,dim,indexL,solver->m_ghosts,pL);
      int pR; _ArrayIndex1D_(_MODEL_NDIMS_,dim,indexR,solver->m_ghosts,pR);
      int q = p*_MODEL_NVARS_;

      /* Modified Rusanov'a_s upwinding scheme */

      udiff[0] = 0.5 * (a_uR[q+0] - a_uL[q+0]);
      udiff[1] = 0.5 * (a_uR[q+1] - a_uL[q+1]);
      udiff[2] = 0.5 * (a_uR[q+2] - a_uL[q+2]);
      udiff[3] = 0.5 * (a_uR[q+3] - a_uL[q+3]);
      udiff[4] = 0.5 * (a_uR[q+4] - a_uL[q+4]);

      _NavierStokes3DRoeAverage_(uavg,_NavierStokes3D_stride_,(a_u+_MODEL_NVARS_*pL),(a_u+_MODEL_NVARS_*pR),param->m_gamma);

      double c, vel[_MODEL_NDIMS_], rho,E,P;
      _NavierStokes3DGetFlowVar_((a_u+_MODEL_NVARS_*pL),_NavierStokes3D_stride_,rho,vel[0],vel[1],vel[2],E,P,param->m_gamma);
      c = sqrt(param->m_gamma*P/rho);
      double alphaL = c + absolute(vel[a_dir]);
      _NavierStokes3DGetFlowVar_((a_u+_MODEL_NVARS_*pR),_NavierStokes3D_stride_,rho,vel[0],vel[1],vel[2],E,P,param->m_gamma);
      c = sqrt(param->m_gamma*P/rho);
      double alphaR = c + absolute(vel[a_dir]);
      _NavierStokes3DGetFlowVar_(uavg,dummy,rho,vel[0],vel[1],vel[2],E,P,param->m_gamma);
      c = sqrt(param->m_gamma*P/rho);
      double alphaavg = c + absolute(vel[a_dir]);

      double kappa  = max(param->m_grav_field_g[pL],param->m_grav_field_g[pR]);
      double alpha  = kappa*max3(alphaL,alphaR,alphaavg);

      a_fI[q+0] = 0.5*(a_fL[q+0]+a_fR[q+0])-alpha*udiff[0];
      a_fI[q+1] = 0.5*(a_fL[q+1]+a_fR[q+1])-alpha*udiff[1];
      a_fI[q+2] = 0.5*(a_fL[q+2]+a_fR[q+2])-alpha*udiff[2];
      a_fI[q+3] = 0.5*(a_fL[q+3]+a_fR[q+3])-alpha*udiff[3];
      a_fI[q+4] = 0.5*(a_fL[q+4]+a_fR[q+4])-alpha*udiff[4];
    }
    _ArrayIncrementIndex_(_MODEL_NDIMS_,bounds_outer,index_outer,done);
  }

  return(0);
}

/*! The Roe upwinding scheme (#NavierStokes3DUpwindRoe) for the partitioned hyperbolic flux that comprises
    of the acoustic waves only (see #NavierStokes3DStiffFlux, #_NavierStokes3DSetStiffFlux_). Thus, only the
    characteristic fields / eigen-modes corresponding to \f$ u\pm a\f$ are used.
    Reference:
    + Ghosh, D., Constantinescu, E. M., Semi-Implicit Time Integration of Atmospheric Flows
      with Characteristic-Based Flux Partitioning, SIAM Journal on Scientific Computing,
      38 (3), 2016, A1848-A1875, http://dx.doi.org/10.1137/15M1044369.
*/
int NavierStokes3DUpwinddFRoe(
                              double  *a_fI, /*!< Computed upwind interface flux */
                              double  *a_fL, /*!< Left-biased reconstructed interface flux */
                              double  *a_fR, /*!< Right-biased reconstructed interface flux */
                              double  *a_uL, /*!< Left-biased reconstructed interface solution */
                              double  *a_uR, /*!< Right-biased reconstructed interface solution */
                              double  *a_u,  /*!< Cell-centered solution */
                              int     a_dir, /*!< Spatial dimension (x, y, or z) */
                              void    *a_s,  /*!< Solver object of type #HyPar */
                              double  a_t   /*!< Current solution time */
                             )
{
  HyPar           *solver = (HyPar*)    a_s;
  NavierStokes3D  *param  = (NavierStokes3D*)  solver->m_physics;
  int             done;

  int     *dim  = solver->m_dim_local;
  double  *uref = param->solution;

  int bounds_outer[_MODEL_NDIMS_], bounds_inter[_MODEL_NDIMS_];
  _ArrayCopy1D3_(dim,bounds_outer,_MODEL_NDIMS_); bounds_outer[a_dir] =  1;
  _ArrayCopy1D3_(dim,bounds_inter,_MODEL_NDIMS_); bounds_inter[a_dir] += 1;
  static double R[_MODEL_NVARS_*_MODEL_NVARS_], D[_MODEL_NVARS_*_MODEL_NVARS_],
                L[_MODEL_NVARS_*_MODEL_NVARS_], DL[_MODEL_NVARS_*_MODEL_NVARS_],
                modA[_MODEL_NVARS_*_MODEL_NVARS_];

  done = 0; int index_outer[3] = {0,0,0}, index_inter[3];
  while (!done) {
    _ArrayCopy1D3_(index_outer,index_inter,_MODEL_NDIMS_);
    for (index_inter[a_dir] = 0; index_inter[a_dir] < bounds_inter[a_dir]; index_inter[a_dir]++) {
      int p; _ArrayIndex1D3_(_MODEL_NDIMS_,bounds_inter,index_inter,0,p);
      int indexL[_MODEL_NDIMS_]; _ArrayCopy1D_(index_inter,indexL,_MODEL_NDIMS_); indexL[a_dir]--;
      int indexR[_MODEL_NDIMS_]; _ArrayCopy1D_(index_inter,indexR,_MODEL_NDIMS_);
      int pL; _ArrayIndex1D_(_MODEL_NDIMS_,dim,indexL,solver->m_ghosts,pL);
      int pR; _ArrayIndex1D_(_MODEL_NDIMS_,dim,indexR,solver->m_ghosts,pR);
      double udiff[_MODEL_NVARS_], uavg[_MODEL_NVARS_],udiss[_MODEL_NVARS_];

      /* Roe'a_s upwinding scheme */

      udiff[0] = 0.5 * (a_uR[_MODEL_NVARS_*p+0] - a_uL[_MODEL_NVARS_*p+0]);
      udiff[1] = 0.5 * (a_uR[_MODEL_NVARS_*p+1] - a_uL[_MODEL_NVARS_*p+1]);
      udiff[2] = 0.5 * (a_uR[_MODEL_NVARS_*p+2] - a_uL[_MODEL_NVARS_*p+2]);
      udiff[3] = 0.5 * (a_uR[_MODEL_NVARS_*p+3] - a_uL[_MODEL_NVARS_*p+3]);
      udiff[4] = 0.5 * (a_uR[_MODEL_NVARS_*p+4] - a_uL[_MODEL_NVARS_*p+4]);

      _NavierStokes3DRoeAverage_        (uavg,_NavierStokes3D_stride_,(uref+_MODEL_NVARS_*pL),(uref+_MODEL_NVARS_*pR),param->m_gamma);
      _NavierStokes3DEigenvalues_       (uavg,dummy,D,param->m_gamma,a_dir);
      _NavierStokes3DLeftEigenvectors_  (uavg,dummy,L,param->m_gamma,a_dir);
      _NavierStokes3DRightEigenvectors_ (uavg,dummy,R,param->m_gamma,a_dir);

      /* Harten'a_s Entropy Fix - Page 362 of Leveque */
      int k;
      double delta = 0.000001, delta2 = delta*delta;
      if (a_dir == _XDIR_) {
        k=0;  D[k] = 0.0;
        k=6;  D[k] = (absolute(D[k]) < delta ? (D[k]*D[k]+delta2)/(2*delta) : absolute(D[k]) );
        k=12; D[k] = 0.0;
        k=18; D[k] = 0.0;
        k=24; D[k] = (absolute(D[k]) < delta ? (D[k]*D[k]+delta2)/(2*delta) : absolute(D[k]) );
      } else if (a_dir == _YDIR_) {
        k=0;  D[k] = 0.0;
        k=6;  D[k] = 0.0;
        k=12; D[k] = (absolute(D[k]) < delta ? (D[k]*D[k]+delta2)/(2*delta) : absolute(D[k]) );
        k=18; D[k] = 0.0;
        k=24; D[k] = (absolute(D[k]) < delta ? (D[k]*D[k]+delta2)/(2*delta) : absolute(D[k]) );
      } else if (a_dir == _ZDIR_) {
        k=0;  D[k] = 0.0;
        k=6;  D[k] = 0.0;
        k=12; D[k] = 0.0;
        k=18; D[k] = (absolute(D[k]) < delta ? (D[k]*D[k]+delta2)/(2*delta) : absolute(D[k]) );
        k=24; D[k] = (absolute(D[k]) < delta ? (D[k]*D[k]+delta2)/(2*delta) : absolute(D[k]) );
      }

      MatMult5(_MODEL_NVARS_,DL,D,L);
      MatMult5(_MODEL_NVARS_,modA,R,DL);
      MatVecMult5(_MODEL_NVARS_,udiss,modA,udiff);

      a_fI[_MODEL_NVARS_*p+0] = 0.5 * (a_fL[_MODEL_NVARS_*p+0]+a_fR[_MODEL_NVARS_*p+0]) - udiss[0];
      a_fI[_MODEL_NVARS_*p+1] = 0.5 * (a_fL[_MODEL_NVARS_*p+1]+a_fR[_MODEL_NVARS_*p+1]) - udiss[1];
      a_fI[_MODEL_NVARS_*p+2] = 0.5 * (a_fL[_MODEL_NVARS_*p+2]+a_fR[_MODEL_NVARS_*p+2]) - udiss[2];
      a_fI[_MODEL_NVARS_*p+3] = 0.5 * (a_fL[_MODEL_NVARS_*p+3]+a_fR[_MODEL_NVARS_*p+3]) - udiss[3];
      a_fI[_MODEL_NVARS_*p+4] = 0.5 * (a_fL[_MODEL_NVARS_*p+4]+a_fR[_MODEL_NVARS_*p+4]) - udiss[4];
    }
    _ArrayIncrementIndex_(_MODEL_NDIMS_,bounds_outer,index_outer,done);
  }

  return(0);
}

/*! The Roe upwinding scheme (#NavierStokes3DUpwindRoe) for the partitioned hyperbolic flux that comprises
    of the entropy waves only (see #NavierStokes3DNonStiffFlux, #_NavierStokes3DSetStiffFlux_). Thus, only the
    characteristic fields / eigen-modes corresponding to \f$u\f$ are used.
    Reference:
    + Ghosh, D., Constantinescu, E. M., Semi-Implicit Time Integration of Atmospheric Flows
      with Characteristic-Based Flux Partitioning, SIAM Journal on Scientific Computing,
      38 (3), 2016, A1848-A1875, http://dx.doi.org/10.1137/15M1044369.
*/
int NavierStokes3DUpwindFdFRoe(
                                double  *a_fI, /*!< Computed upwind interface flux */
                                double  *a_fL, /*!< Left-biased reconstructed interface flux */
                                double  *a_fR, /*!< Right-biased reconstructed interface flux */
                                double  *a_uL, /*!< Left-biased reconstructed interface solution */
                                double  *a_uR, /*!< Right-biased reconstructed interface solution */
                                double  *a_u,  /*!< Cell-centered solution */
                                int     a_dir, /*!< Spatial dimension (x, y, or z) */
                                void    *a_s,  /*!< Solver object of type #HyPar */
                                double  a_t   /*!< Current solution time */
                              )
{
  HyPar           *solver = (HyPar*)    a_s;
  NavierStokes3D  *param  = (NavierStokes3D*)  solver->m_physics;
  int             done;

  int     *dim  = solver->m_dim_local;
  double  *uref = param->solution;

  int bounds_outer[_MODEL_NDIMS_], bounds_inter[_MODEL_NDIMS_];
  _ArrayCopy1D3_(dim,bounds_outer,_MODEL_NDIMS_); bounds_outer[a_dir] =  1;
  _ArrayCopy1D3_(dim,bounds_inter,_MODEL_NDIMS_); bounds_inter[a_dir] += 1;
  static double R[_MODEL_NVARS_*_MODEL_NVARS_], D[_MODEL_NVARS_*_MODEL_NVARS_],
                L[_MODEL_NVARS_*_MODEL_NVARS_], DL[_MODEL_NVARS_*_MODEL_NVARS_],
                modA[_MODEL_NVARS_*_MODEL_NVARS_];

  done = 0; int index_outer[3] = {0,0,0}, index_inter[3];
  while (!done) {
    _ArrayCopy1D3_(index_outer,index_inter,_MODEL_NDIMS_);
    for (index_inter[a_dir] = 0; index_inter[a_dir] < bounds_inter[a_dir]; index_inter[a_dir]++) {
      int p; _ArrayIndex1D3_(_MODEL_NDIMS_,bounds_inter,index_inter,0,p);
      int indexL[_MODEL_NDIMS_]; _ArrayCopy1D_(index_inter,indexL,_MODEL_NDIMS_); indexL[a_dir]--;
      int indexR[_MODEL_NDIMS_]; _ArrayCopy1D_(index_inter,indexR,_MODEL_NDIMS_);
      int pL; _ArrayIndex1D_(_MODEL_NDIMS_,dim,indexL,solver->m_ghosts,pL);
      int pR; _ArrayIndex1D_(_MODEL_NDIMS_,dim,indexR,solver->m_ghosts,pR);
      int k;
      double udiff[_MODEL_NVARS_],uavg[_MODEL_NVARS_],udiss[_MODEL_NVARS_],
             udiss_total[_MODEL_NVARS_],udiss_stiff[_MODEL_NVARS_];
      double delta = 0.000001, delta2 = delta*delta;

      /* Roe'a_s upwinding scheme */

      udiff[0] = 0.5 * (a_uR[_MODEL_NVARS_*p+0] - a_uL[_MODEL_NVARS_*p+0]);
      udiff[1] = 0.5 * (a_uR[_MODEL_NVARS_*p+1] - a_uL[_MODEL_NVARS_*p+1]);
      udiff[2] = 0.5 * (a_uR[_MODEL_NVARS_*p+2] - a_uL[_MODEL_NVARS_*p+2]);
      udiff[3] = 0.5 * (a_uR[_MODEL_NVARS_*p+3] - a_uL[_MODEL_NVARS_*p+3]);
      udiff[4] = 0.5 * (a_uR[_MODEL_NVARS_*p+4] - a_uL[_MODEL_NVARS_*p+4]);

      /* Compute total dissipation */
      _NavierStokes3DRoeAverage_        (uavg,_NavierStokes3D_stride_,(a_u+_MODEL_NVARS_*pL),(a_u+_MODEL_NVARS_*pR),param->m_gamma);
      _NavierStokes3DEigenvalues_       (uavg,dummy,D,param->m_gamma,a_dir);
      _NavierStokes3DLeftEigenvectors_  (uavg,dummy,L,param->m_gamma,a_dir);
      _NavierStokes3DRightEigenvectors_ (uavg,dummy,R,param->m_gamma,a_dir);
      k=0;  D[k] = (absolute(D[k]) < delta ? (D[k]*D[k]+delta2)/(2*delta) : absolute(D[k]) );
      k=6;  D[k] = (absolute(D[k]) < delta ? (D[k]*D[k]+delta2)/(2*delta) : absolute(D[k]) );
      k=12; D[k] = (absolute(D[k]) < delta ? (D[k]*D[k]+delta2)/(2*delta) : absolute(D[k]) );
      k=18; D[k] = (absolute(D[k]) < delta ? (D[k]*D[k]+delta2)/(2*delta) : absolute(D[k]) );
      k=24; D[k] = (absolute(D[k]) < delta ? (D[k]*D[k]+delta2)/(2*delta) : absolute(D[k]) );
      MatMult5(_MODEL_NVARS_,DL,D,L);
      MatMult5(_MODEL_NVARS_,modA,R,DL);
      MatVecMult5(_MODEL_NVARS_,udiss_total,modA,udiff);

      /* Compute dissipation corresponding to acoustic modes */
      _NavierStokes3DRoeAverage_        (uavg,_NavierStokes3D_stride_,(uref+_MODEL_NVARS_*pL),(uref+_MODEL_NVARS_*pR),param->m_gamma);
      _NavierStokes3DEigenvalues_       (uavg,dummy,D,param->m_gamma,a_dir);
      _NavierStokes3DLeftEigenvectors_  (uavg,dummy,L,param->m_gamma,a_dir);
      _NavierStokes3DRightEigenvectors_ (uavg,dummy,R,param->m_gamma,a_dir);
      if (a_dir == _XDIR_) {
        k=0;  D[k] = 0.0;
        k=6;  D[k] = (absolute(D[k]) < delta ? (D[k]*D[k]+delta2)/(2*delta) : absolute(D[k]) );
        k=12; D[k] = 0.0;
        k=18; D[k] = 0.0;
        k=24; D[k] = (absolute(D[k]) < delta ? (D[k]*D[k]+delta2)/(2*delta) : absolute(D[k]) );
      } else if (a_dir == _YDIR_) {
        k=0;  D[k] = 0.0;
        k=6;  D[k] = 0.0;
        k=12; D[k] = (absolute(D[k]) < delta ? (D[k]*D[k]+delta2)/(2*delta) : absolute(D[k]) );
        k=18; D[k] = 0.0;
        k=24; D[k] = (absolute(D[k]) < delta ? (D[k]*D[k]+delta2)/(2*delta) : absolute(D[k]) );
      } else if (a_dir == _ZDIR_) {
        k=0;  D[k] = 0.0;
        k=6;  D[k] = 0.0;
        k=12; D[k] = 0.0;
        k=18; D[k] = (absolute(D[k]) < delta ? (D[k]*D[k]+delta2)/(2*delta) : absolute(D[k]) );
        k=24; D[k] = (absolute(D[k]) < delta ? (D[k]*D[k]+delta2)/(2*delta) : absolute(D[k]) );
      }
      MatMult5(_MODEL_NVARS_,DL,D,L);
      MatMult5(_MODEL_NVARS_,modA,R,DL);
      MatVecMult5(_MODEL_NVARS_,udiss_stiff,modA,udiff);

     /* Compute the dissipation term for the entropy modes */
      _ArraySubtract1D_(udiss,udiss_total,udiss_stiff,_MODEL_NVARS_);

      a_fI[_MODEL_NVARS_*p+0] = 0.5 * (a_fL[_MODEL_NVARS_*p+0]+a_fR[_MODEL_NVARS_*p+0]) - udiss[0];
      a_fI[_MODEL_NVARS_*p+1] = 0.5 * (a_fL[_MODEL_NVARS_*p+1]+a_fR[_MODEL_NVARS_*p+1]) - udiss[1];
      a_fI[_MODEL_NVARS_*p+2] = 0.5 * (a_fL[_MODEL_NVARS_*p+2]+a_fR[_MODEL_NVARS_*p+2]) - udiss[2];
      a_fI[_MODEL_NVARS_*p+3] = 0.5 * (a_fL[_MODEL_NVARS_*p+3]+a_fR[_MODEL_NVARS_*p+3]) - udiss[3];
      a_fI[_MODEL_NVARS_*p+4] = 0.5 * (a_fL[_MODEL_NVARS_*p+4]+a_fR[_MODEL_NVARS_*p+4]) - udiss[4];
    }
    _ArrayIncrementIndex_(_MODEL_NDIMS_,bounds_outer,index_outer,done);
  }

  return(0);
}

/*! Modified Rusanov's upwinding scheme: NavierStokes3DUpwindRusanov() modified as described in the
    following paper (for consistent characteristic-based splitting):
    + Ghosh, D., Constantinescu, E. M., "Semi-Implicit Time Integration of Atmospheric Flows with
      Characteristic-Based Flux Partitioning", Submitted (http://arxiv.org/abs/1510.05751).

*/
int NavierStokes3DUpwindRusanovModified(
                                        double  *a_fI, /*!< Computed upwind interface flux */
                                        double  *a_fL, /*!< Left-biased reconstructed interface flux */
                                        double  *a_fR, /*!< Right-biased reconstructed interface flux */
                                        double  *a_uL, /*!< Left-biased reconstructed interface solution */
                                        double  *a_uR, /*!< Right-biased reconstructed interface solution */
                                        double  *a_u,  /*!< Cell-centered solution */
                                        int     a_dir, /*!< Spatial dimension (x,y, or z) */
                                        void    *a_s,  /*!< Solver object of type #HyPar */
                                        double  a_t   /*!< Current solution time */
                                       )
{
  HyPar           *solver = (HyPar*)          a_s;
  NavierStokes3D  *param  = (NavierStokes3D*) solver->m_physics;
  int             *dim    = solver->m_dim_local, done;

  static int bounds_outer[_MODEL_NDIMS_], bounds_inter[_MODEL_NDIMS_];
  bounds_outer[0] = dim[0]; bounds_outer[1] = dim[1]; bounds_outer[2] = dim[2]; bounds_outer[a_dir] = 1;
  bounds_inter[0] = dim[0]; bounds_inter[1] = dim[1]; bounds_inter[2] = dim[2]; bounds_inter[a_dir]++;

  static double R[_MODEL_NVARS_*_MODEL_NVARS_], D[_MODEL_NVARS_*_MODEL_NVARS_],
                L[_MODEL_NVARS_*_MODEL_NVARS_], DL[_MODEL_NVARS_*_MODEL_NVARS_],
                modA[_MODEL_NVARS_*_MODEL_NVARS_];

  static int indexL[_MODEL_NDIMS_], indexR[_MODEL_NDIMS_],
             index_outer[_MODEL_NDIMS_], index_inter[_MODEL_NDIMS_];

  static double udiff[_MODEL_NVARS_],uavg[_MODEL_NVARS_],udiss[_MODEL_NVARS_];

  done = 0; _ArraySetValue_(index_outer,_MODEL_NDIMS_,0);
  while (!done) {
    _ArrayCopy1D_(index_outer,index_inter,_MODEL_NDIMS_);
    for (index_inter[a_dir] = 0; index_inter[a_dir] < bounds_inter[a_dir]; index_inter[a_dir]++) {

      int p; _ArrayIndex1D_(_MODEL_NDIMS_,bounds_inter,index_inter,0,p);
      _ArrayCopy1D_(index_inter,indexL,_MODEL_NDIMS_); indexL[a_dir]--;
      _ArrayCopy1D_(index_inter,indexR,_MODEL_NDIMS_);
      int pL; _ArrayIndex1D_(_MODEL_NDIMS_,dim,indexL,solver->m_ghosts,pL);
      int pR; _ArrayIndex1D_(_MODEL_NDIMS_,dim,indexR,solver->m_ghosts,pR);
      int q = p*_MODEL_NVARS_;

      /* Modified Rusanov'a_s upwinding scheme */

      udiff[0] = 0.5 * (a_uR[q+0] - a_uL[q+0]);
      udiff[1] = 0.5 * (a_uR[q+1] - a_uL[q+1]);
      udiff[2] = 0.5 * (a_uR[q+2] - a_uL[q+2]);
      udiff[3] = 0.5 * (a_uR[q+3] - a_uL[q+3]);
      udiff[4] = 0.5 * (a_uR[q+4] - a_uL[q+4]);

      _NavierStokes3DRoeAverage_        (uavg,_NavierStokes3D_stride_,(a_u+_MODEL_NVARS_*pL),(a_u+_MODEL_NVARS_*pR),param->m_gamma);
      _NavierStokes3DLeftEigenvectors_  (uavg,dummy,L,param->m_gamma,a_dir);
      _NavierStokes3DRightEigenvectors_ (uavg,dummy,R,param->m_gamma,a_dir);

      double c, vel[_MODEL_NDIMS_], rho,E,P;

      _NavierStokes3DGetFlowVar_((a_u+_MODEL_NVARS_*pL),_NavierStokes3D_stride_,rho,vel[0],vel[1],vel[2],E,P,param->m_gamma);
      c = sqrt(param->m_gamma*P/rho);
      double alphaL = c + absolute(vel[a_dir]);
      double betaL = absolute(vel[a_dir]);

      _NavierStokes3DGetFlowVar_((a_u+_MODEL_NVARS_*pR),_NavierStokes3D_stride_,rho,vel[0],vel[1],vel[2],E,P,param->m_gamma);
      c = sqrt(param->m_gamma*P/rho);
      double alphaR = c + absolute(vel[a_dir]);
      double betaR = absolute(vel[a_dir]);

      _NavierStokes3DGetFlowVar_(uavg,dummy,rho,vel[0],vel[1],vel[2],E,P,param->m_gamma);
      c = sqrt(param->m_gamma*P/rho);
      double alphaavg = c + absolute(vel[a_dir]);
      double betaavg = absolute(vel[a_dir]);

      double kappa  = max(param->m_grav_field_g[pL],param->m_grav_field_g[pR]);
      double alpha  = kappa*max3(alphaL,alphaR,alphaavg);
      double beta   = kappa*max3(betaL,betaR,betaavg);

      _ArraySetValue_(D,_MODEL_NVARS_*_MODEL_NVARS_,0.0);
      if (a_dir == _XDIR_) {
        D[0]  = beta;
        D[6]  = alpha;
        D[12] = beta;
        D[18] = beta;
        D[24] = alpha;
      } else if (a_dir == _YDIR_) {
        D[0]  = beta;
        D[6]  = beta;
        D[12] = alpha;
        D[18] = beta;
        D[24] = alpha;
      } else if (a_dir == _ZDIR_) {
        D[0]  = beta;
        D[6]  = beta;
        D[12] = beta;
        D[18] = alpha;
        D[24] = alpha;
      }
      MatMult5    (_MODEL_NVARS_,DL,D,L);
      MatMult5    (_MODEL_NVARS_,modA,R,DL);
      MatVecMult5 (_MODEL_NVARS_,udiss,modA,udiff);

      a_fI[q+0] = 0.5*(a_fL[q+0]+a_fR[q+0])-udiss[0];
      a_fI[q+1] = 0.5*(a_fL[q+1]+a_fR[q+1])-udiss[1];
      a_fI[q+2] = 0.5*(a_fL[q+2]+a_fR[q+2])-udiss[2];
      a_fI[q+3] = 0.5*(a_fL[q+3]+a_fR[q+3])-udiss[3];
      a_fI[q+4] = 0.5*(a_fL[q+4]+a_fR[q+4])-udiss[4];
    }
    _ArrayIncrementIndex_(_MODEL_NDIMS_,bounds_outer,index_outer,done);
  }

  return(0);
}

/*! The modified Rusanov upwinding scheme (NavierStokes3DUpwindRusanovModified()) for the partitioned hyperbolic flux that comprises
    of the acoustic waves only (see #NavierStokes3DStiffFlux, #_NavierStokes3DSetStiffFlux_). Thus, only the
    characteristic fields / eigen-modes corresponding to \f$ u\pm a\f$ are used.
    Reference:
    + Ghosh, D., Constantinescu, E. M., Semi-Implicit Time Integration of Atmospheric Flows
      with Characteristic-Based Flux Partitioning, SIAM Journal on Scientific Computing,
      38 (3), 2016, A1848-A1875, http://dx.doi.org/10.1137/15M1044369.
*/
int NavierStokes3DUpwinddFRusanovModified(
                                          double  *a_fI, /*!< Computed upwind interface flux */
                                          double  *a_fL, /*!< Left-biased reconstructed interface flux */
                                          double  *a_fR, /*!< Right-biased reconstructed interface flux */
                                          double  *a_uL, /*!< Left-biased reconstructed interface solution */
                                          double  *a_uR, /*!< Right-biased reconstructed interface solution */
                                          double  *a_u,  /*!< Cell-centered solution */
                                          int     a_dir, /*!< Spatial dimension (x,y, or z) */
                                          void    *a_s,  /*!< Solver object of type #HyPar */
                                          double  a_t   /*!< Current solution time */
                                         )
{
  HyPar           *solver = (HyPar*)          a_s;
  NavierStokes3D  *param  = (NavierStokes3D*) solver->m_physics;
  int             *dim    = solver->m_dim_local, done;
  double          *uref   = param->solution;

  static int bounds_outer[_MODEL_NDIMS_], bounds_inter[_MODEL_NDIMS_];
  bounds_outer[0] = dim[0]; bounds_outer[1] = dim[1]; bounds_outer[2] = dim[2]; bounds_outer[a_dir] = 1;
  bounds_inter[0] = dim[0]; bounds_inter[1] = dim[1]; bounds_inter[2] = dim[2]; bounds_inter[a_dir]++;

  static double R[_MODEL_NVARS_*_MODEL_NVARS_], D[_MODEL_NVARS_*_MODEL_NVARS_],
                L[_MODEL_NVARS_*_MODEL_NVARS_], DL[_MODEL_NVARS_*_MODEL_NVARS_],
                modA[_MODEL_NVARS_*_MODEL_NVARS_];

  static int indexL[_MODEL_NDIMS_], indexR[_MODEL_NDIMS_],
             index_outer[_MODEL_NDIMS_], index_inter[_MODEL_NDIMS_];

  static double udiff[_MODEL_NVARS_],uavg[_MODEL_NVARS_],udiss[_MODEL_NVARS_];

  done = 0; _ArraySetValue_(index_outer,_MODEL_NDIMS_,0);
  while (!done) {
    _ArrayCopy1D_(index_outer,index_inter,_MODEL_NDIMS_);
    for (index_inter[a_dir] = 0; index_inter[a_dir] < bounds_inter[a_dir]; index_inter[a_dir]++) {

      int p; _ArrayIndex1D_(_MODEL_NDIMS_,bounds_inter,index_inter,0,p);
      _ArrayCopy1D_(index_inter,indexL,_MODEL_NDIMS_); indexL[a_dir]--;
      _ArrayCopy1D_(index_inter,indexR,_MODEL_NDIMS_);
      int pL; _ArrayIndex1D_(_MODEL_NDIMS_,dim,indexL,solver->m_ghosts,pL);
      int pR; _ArrayIndex1D_(_MODEL_NDIMS_,dim,indexR,solver->m_ghosts,pR);
      int q = p*_MODEL_NVARS_;

      /* Modified Rusanov'a_s upwinding scheme */

      udiff[0] = 0.5 * (a_uR[q+0] - a_uL[q+0]);
      udiff[1] = 0.5 * (a_uR[q+1] - a_uL[q+1]);
      udiff[2] = 0.5 * (a_uR[q+2] - a_uL[q+2]);
      udiff[3] = 0.5 * (a_uR[q+3] - a_uL[q+3]);
      udiff[4] = 0.5 * (a_uR[q+4] - a_uL[q+4]);

      _NavierStokes3DRoeAverage_        (uavg,_NavierStokes3D_stride_,(uref+_MODEL_NVARS_*pL),(uref+_MODEL_NVARS_*pR),param->m_gamma);
      _NavierStokes3DLeftEigenvectors_  (uavg,dummy,L,param->m_gamma,a_dir);
      _NavierStokes3DRightEigenvectors_ (uavg,dummy,R,param->m_gamma,a_dir);

      double c, vel[_MODEL_NDIMS_], rho,E,P;

      _NavierStokes3DGetFlowVar_((uref+_MODEL_NVARS_*pL),_NavierStokes3D_stride_,rho,vel[0],vel[1],vel[2],E,P,param->m_gamma);
      c = sqrt(param->m_gamma*P/rho);
      double alphaL = c + absolute(vel[a_dir]);

      _NavierStokes3DGetFlowVar_((uref+_MODEL_NVARS_*pR),_NavierStokes3D_stride_,rho,vel[0],vel[1],vel[2],E,P,param->m_gamma);
      c = sqrt(param->m_gamma*P/rho);
      double alphaR = c + absolute(vel[a_dir]);

      _NavierStokes3DGetFlowVar_(uavg,dummy,rho,vel[0],vel[1],vel[2],E,P,param->m_gamma);
      c = sqrt(param->m_gamma*P/rho);
      double alphaavg = c + absolute(vel[a_dir]);

      double kappa  = max(param->m_grav_field_g[pL],param->m_grav_field_g[pR]);
      double alpha  = kappa*max3(alphaL,alphaR,alphaavg);

      _ArraySetValue_(D,_MODEL_NVARS_*_MODEL_NVARS_,0.0);
      if (a_dir == _XDIR_) {
        D[6]  = alpha;
        D[24] = alpha;
      } else if (a_dir == _YDIR_) {
        D[12] = alpha;
        D[24] = alpha;
      } else if (a_dir == _ZDIR_) {
        D[18] = alpha;
        D[24] = alpha;
      }
      MatMult5    (_MODEL_NVARS_,DL,D,L);
      MatMult5    (_MODEL_NVARS_,modA,R,DL);
      MatVecMult5 (_MODEL_NVARS_,udiss,modA,udiff);

      a_fI[q+0] = 0.5*(a_fL[q+0]+a_fR[q+0])-udiss[0];
      a_fI[q+1] = 0.5*(a_fL[q+1]+a_fR[q+1])-udiss[1];
      a_fI[q+2] = 0.5*(a_fL[q+2]+a_fR[q+2])-udiss[2];
      a_fI[q+3] = 0.5*(a_fL[q+3]+a_fR[q+3])-udiss[3];
      a_fI[q+4] = 0.5*(a_fL[q+4]+a_fR[q+4])-udiss[4];
    }
    _ArrayIncrementIndex_(_MODEL_NDIMS_,bounds_outer,index_outer,done);
  }

  return(0);
}

/*! The modified Rusanov upwinding scheme (NavierStokes3DUpwindRusanovModified()) for the partitioned hyperbolic flux that comprises
    of the entropy waves only (see #NavierStokes3DNonStiffFlux, #_NavierStokes3DSetNonStiffFlux_). Thus, only the
    characteristic fields / eigen-modes corresponding to \f$u\f$ are used.
    Reference:
    + Ghosh, D., Constantinescu, E. M., Semi-Implicit Time Integration of Atmospheric Flows
      with Characteristic-Based Flux Partitioning, SIAM Journal on Scientific Computing,
      38 (3), 2016, A1848-A1875, http://dx.doi.org/10.1137/15M1044369.
*/
int NavierStokes3DUpwindFdFRusanovModified(
                                            double  *a_fI, /*!< Computed upwind interface flux */
                                            double  *a_fL, /*!< Left-biased reconstructed interface flux */
                                            double  *a_fR, /*!< Right-biased reconstructed interface flux */
                                            double  *a_uL, /*!< Left-biased reconstructed interface solution */
                                            double  *a_uR, /*!< Right-biased reconstructed interface solution */
                                            double  *a_u,  /*!< Cell-centered solution */
                                            int     a_dir, /*!< Spatial dimension (x,y, or z) */
                                            void    *a_s,  /*!< Solver object of type #HyPar */
                                            double  a_t   /*!< Current solution time */
                                          )
{
  HyPar           *solver = (HyPar*)          a_s;
  NavierStokes3D  *param  = (NavierStokes3D*) solver->m_physics;
  int             *dim    = solver->m_dim_local, done;
  double          *uref   = param->solution;

  static int bounds_outer[_MODEL_NDIMS_], bounds_inter[_MODEL_NDIMS_];
  bounds_outer[0] = dim[0]; bounds_outer[1] = dim[1]; bounds_outer[2] = dim[2]; bounds_outer[a_dir] = 1;
  bounds_inter[0] = dim[0]; bounds_inter[1] = dim[1]; bounds_inter[2] = dim[2]; bounds_inter[a_dir]++;

  static double R[_MODEL_NVARS_*_MODEL_NVARS_], D[_MODEL_NVARS_*_MODEL_NVARS_],
                L[_MODEL_NVARS_*_MODEL_NVARS_], DL[_MODEL_NVARS_*_MODEL_NVARS_],
                modA[_MODEL_NVARS_*_MODEL_NVARS_];

  static int indexL[_MODEL_NDIMS_], indexR[_MODEL_NDIMS_],
             index_outer[_MODEL_NDIMS_], index_inter[_MODEL_NDIMS_];

  static double udiff[_MODEL_NVARS_],uavg[_MODEL_NVARS_],udiss[_MODEL_NVARS_];

  done = 0; _ArraySetValue_(index_outer,_MODEL_NDIMS_,0);
  while (!done) {
    _ArrayCopy1D_(index_outer,index_inter,_MODEL_NDIMS_);
    for (index_inter[a_dir] = 0; index_inter[a_dir] < bounds_inter[a_dir]; index_inter[a_dir]++) {

      int p; _ArrayIndex1D_(_MODEL_NDIMS_,bounds_inter,index_inter,0,p);
      _ArrayCopy1D_(index_inter,indexL,_MODEL_NDIMS_); indexL[a_dir]--;
      _ArrayCopy1D_(index_inter,indexR,_MODEL_NDIMS_);
      int pL; _ArrayIndex1D_(_MODEL_NDIMS_,dim,indexL,solver->m_ghosts,pL);
      int pR; _ArrayIndex1D_(_MODEL_NDIMS_,dim,indexR,solver->m_ghosts,pR);
      int q = p*_MODEL_NVARS_;
      double udiff[_MODEL_NVARS_], udiss[_MODEL_NVARS_], uavg[_MODEL_NVARS_],
             udiss_total[_MODEL_NVARS_],udiss_acoustic[_MODEL_NVARS_];

      /* Modified Rusanov'a_s upwinding scheme */
      /* Modified Rusanov'a_s upwinding scheme */
      double c, vel[_MODEL_NDIMS_], rho,E,P, alphaL, alphaR, alphaavg,
             betaL, betaR, betaavg, alpha, beta,
             kappa  = max(param->m_grav_field_g[pL],param->m_grav_field_g[pR]);

      udiff[0] = 0.5 * (a_uR[q+0] - a_uL[q+0]);
      udiff[1] = 0.5 * (a_uR[q+1] - a_uL[q+1]);
      udiff[2] = 0.5 * (a_uR[q+2] - a_uL[q+2]);
      udiff[3] = 0.5 * (a_uR[q+3] - a_uL[q+3]);
      udiff[4] = 0.5 * (a_uR[q+4] - a_uL[q+4]);

      /* Compute total dissipation */
      _NavierStokes3DRoeAverage_        (uavg,_NavierStokes3D_stride_,(a_u+_MODEL_NVARS_*pL),(a_u+_MODEL_NVARS_*pR),param->m_gamma);
      _NavierStokes3DLeftEigenvectors_  (uavg,dummy,L,param->m_gamma,a_dir);
      _NavierStokes3DRightEigenvectors_ (uavg,dummy,R,param->m_gamma,a_dir);

      _NavierStokes3DGetFlowVar_((a_u+_MODEL_NVARS_*pL),_NavierStokes3D_stride_,rho,vel[0],vel[1],vel[2],E,P,param->m_gamma);
      c = sqrt(param->m_gamma*P/rho);
      alphaL = c + absolute(vel[a_dir]);
      betaL = absolute(vel[a_dir]);

      _NavierStokes3DGetFlowVar_((a_u+_MODEL_NVARS_*pR),_NavierStokes3D_stride_,rho,vel[0],vel[1],vel[2],E,P,param->m_gamma);
      c = sqrt(param->m_gamma*P/rho);
      alphaR = c + absolute(vel[a_dir]);
      betaR = absolute(vel[a_dir]);

      _NavierStokes3DGetFlowVar_(uavg,dummy,rho,vel[0],vel[1],vel[2],E,P,param->m_gamma);
      c = sqrt(param->m_gamma*P/rho);
      alphaavg = c + absolute(vel[a_dir]);
      betaavg = absolute(vel[a_dir]);

      kappa  = max(param->m_grav_field_g[pL],param->m_grav_field_g[pR]);
      alpha  = kappa*max3(alphaL,alphaR,alphaavg);
      beta   = kappa*max3(betaL,betaR,betaavg);

      _ArraySetValue_(D,_MODEL_NVARS_*_MODEL_NVARS_,0.0);
      if (a_dir == _XDIR_) {
        D[0]  = beta;
        D[6]  = alpha;
        D[12] = beta;
        D[18] = beta;
        D[24] = alpha;
      } else if (a_dir == _YDIR_) {
        D[0]  = beta;
        D[6]  = beta;
        D[12] = alpha;
        D[18] = beta;
        D[24] = alpha;
      } else if (a_dir == _ZDIR_) {
        D[0]  = beta;
        D[6]  = beta;
        D[12] = beta;
        D[18] = alpha;
        D[24] = alpha;
      }
      MatMult5    (_MODEL_NVARS_,DL,D,L);
      MatMult5    (_MODEL_NVARS_,modA,R,DL);
      MatVecMult5 (_MODEL_NVARS_,udiss_total,modA,udiff);

      /* Compute dissipation for the linearized acoustic modes */
      _NavierStokes3DRoeAverage_        (uavg,_NavierStokes3D_stride_,(uref+_MODEL_NVARS_*pL),(uref+_MODEL_NVARS_*pR),param->m_gamma);
      _NavierStokes3DLeftEigenvectors_  (uavg,dummy,L,param->m_gamma,a_dir);
      _NavierStokes3DRightEigenvectors_ (uavg,dummy,R,param->m_gamma,a_dir);

      _NavierStokes3DGetFlowVar_((uref+_MODEL_NVARS_*pL),_NavierStokes3D_stride_,rho,vel[0],vel[1],vel[2],E,P,param->m_gamma);
      c = sqrt(param->m_gamma*P/rho);
      alphaL = c + absolute(vel[a_dir]);

      _NavierStokes3DGetFlowVar_((uref+_MODEL_NVARS_*pR),_NavierStokes3D_stride_,rho,vel[0],vel[1],vel[2],E,P,param->m_gamma);
      c = sqrt(param->m_gamma*P/rho);
      alphaR = c + absolute(vel[a_dir]);

      _NavierStokes3DGetFlowVar_(uavg,dummy,rho,vel[0],vel[1],vel[2],E,P,param->m_gamma);
      c = sqrt(param->m_gamma*P/rho);
      alphaavg = c + absolute(vel[a_dir]);

      kappa  = max(param->m_grav_field_g[pL],param->m_grav_field_g[pR]);
      alpha  = kappa*max3(alphaL,alphaR,alphaavg);

      _ArraySetValue_(D,_MODEL_NVARS_*_MODEL_NVARS_,0.0);
      if (a_dir == _XDIR_) {
        D[6]  = alpha;
        D[24] = alpha;
      } else if (a_dir == _YDIR_) {
        D[12] = alpha;
        D[24] = alpha;
      } else if (a_dir == _ZDIR_) {
        D[18] = alpha;
        D[24] = alpha;
      }
      MatMult5    (_MODEL_NVARS_,DL,D,L);
      MatMult5    (_MODEL_NVARS_,modA,R,DL);
      MatVecMult5 (_MODEL_NVARS_,udiss_acoustic,modA,udiff);

      /* Compute dissipation for the entropy modes */
      _ArraySubtract1D_(udiss,udiss_total,udiss_acoustic,_MODEL_NVARS_);

      a_fI[q+0] = 0.5*(a_fL[q+0]+a_fR[q+0])-udiss[0];
      a_fI[q+1] = 0.5*(a_fL[q+1]+a_fR[q+1])-udiss[1];
      a_fI[q+2] = 0.5*(a_fL[q+2]+a_fR[q+2])-udiss[2];
      a_fI[q+3] = 0.5*(a_fL[q+3]+a_fR[q+3])-udiss[3];
      a_fI[q+4] = 0.5*(a_fL[q+4]+a_fR[q+4])-udiss[4];
    }
    _ArrayIncrementIndex_(_MODEL_NDIMS_,bounds_outer,index_outer,done);
  }

  return(0);
}

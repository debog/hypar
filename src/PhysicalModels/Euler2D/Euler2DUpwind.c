#include <stdlib.h>
#include <math.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <physicalmodels/euler2d.h>
#include <mathfunctions.h>
#include <matmult_native.h>
#include <hypar.h>

int Euler2DUpwindRoe(double *a_fI,double *a_fL,double *a_fR,double *a_uL,double *a_uR,double *a_u,int a_dir,void *a_s,double a_t)
{
  HyPar    *solver = (HyPar*)    a_s;
  Euler2D  *param  = (Euler2D*)  solver->m_physics;
  int      done;

  int *dim  = solver->m_dim_local;

  int bounds_outer[2], bounds_inter[2];
  bounds_outer[0] = dim[0]; bounds_outer[1] = dim[1]; bounds_outer[a_dir] = 1;
  bounds_inter[0] = dim[0]; bounds_inter[1] = dim[1]; bounds_inter[a_dir]++;
  static double R[_MODEL_NVARS_*_MODEL_NVARS_], D[_MODEL_NVARS_*_MODEL_NVARS_],
                L[_MODEL_NVARS_*_MODEL_NVARS_], DL[_MODEL_NVARS_*_MODEL_NVARS_],
                modA[_MODEL_NVARS_*_MODEL_NVARS_];

  done = 0; int index_outer[2] = {0,0}; int index_inter[2];
  while (!done) {
    index_inter[0] = index_outer[0]; index_inter[1] = index_outer[1];
    for (index_inter[a_dir] = 0; index_inter[a_dir] < bounds_inter[a_dir]; index_inter[a_dir]++) {
      int p; _ArrayIndex1D2_(_MODEL_NDIMS_,bounds_inter,index_inter,0,p);
      double udiff[_MODEL_NVARS_], uavg[_MODEL_NVARS_],udiss[_MODEL_NVARS_];

      /* Roe'a_s upwinding scheme */

      udiff[0] = 0.5 * (a_uR[_MODEL_NVARS_*p+0] - a_uL[_MODEL_NVARS_*p+0]);
      udiff[1] = 0.5 * (a_uR[_MODEL_NVARS_*p+1] - a_uL[_MODEL_NVARS_*p+1]);
      udiff[2] = 0.5 * (a_uR[_MODEL_NVARS_*p+2] - a_uL[_MODEL_NVARS_*p+2]);
      udiff[3] = 0.5 * (a_uR[_MODEL_NVARS_*p+3] - a_uL[_MODEL_NVARS_*p+3]);

      _Euler2DRoeAverage_(uavg,(a_uL+_MODEL_NVARS_*p),(a_uR+_MODEL_NVARS_*p),param);

      _Euler2DEigenvalues_(uavg,D,param,a_dir);
      _Euler2DLeftEigenvectors_(uavg,L,param,a_dir);
      _Euler2DRightEigenvectors_(uavg,R,param,a_dir);

       /* Harten'a_s Entropy Fix - Page 362 of Leveque */
      int k;
      double delta = 0.000001, delta2 = delta*delta;
      k=0;  D[k] = (absolute(D[k]) < delta ? (D[k]*D[k]+delta2)/(2*delta) : absolute(D[k]) );
      k=5;  D[k] = (absolute(D[k]) < delta ? (D[k]*D[k]+delta2)/(2*delta) : absolute(D[k]) );
      k=10; D[k] = (absolute(D[k]) < delta ? (D[k]*D[k]+delta2)/(2*delta) : absolute(D[k]) );
      k=15; D[k] = (absolute(D[k]) < delta ? (D[k]*D[k]+delta2)/(2*delta) : absolute(D[k]) );

      MatMult4(_MODEL_NVARS_,DL,D,L);
      MatMult4(_MODEL_NVARS_,modA,R,DL);
      MatVecMult4(_MODEL_NVARS_,udiss,modA,udiff);

      a_fI[_MODEL_NVARS_*p+0] = 0.5 * (a_fL[_MODEL_NVARS_*p+0]+a_fR[_MODEL_NVARS_*p+0]) - udiss[0];
      a_fI[_MODEL_NVARS_*p+1] = 0.5 * (a_fL[_MODEL_NVARS_*p+1]+a_fR[_MODEL_NVARS_*p+1]) - udiss[1];
      a_fI[_MODEL_NVARS_*p+2] = 0.5 * (a_fL[_MODEL_NVARS_*p+2]+a_fR[_MODEL_NVARS_*p+2]) - udiss[2];
      a_fI[_MODEL_NVARS_*p+3] = 0.5 * (a_fL[_MODEL_NVARS_*p+3]+a_fR[_MODEL_NVARS_*p+3]) - udiss[3];
    }
    _ArrayIncrementIndex_(_MODEL_NDIMS_,bounds_outer,index_outer,done);
  }

  return(0);
}

int Euler2DUpwindRF(double *a_fI,double *a_fL,double *a_fR,double *a_uL,double *a_uR,double *a_u,int a_dir,void *a_s,double a_t)
{
  HyPar    *solver = (HyPar*)    a_s;
  Euler2D  *param  = (Euler2D*)  solver->m_physics;
  int      done,k;

  int *dim  = solver->m_dim_local;

  int bounds_outer[2], bounds_inter[2];
  bounds_outer[0] = dim[0]; bounds_outer[1] = dim[1]; bounds_outer[a_dir] = 1;
  bounds_inter[0] = dim[0]; bounds_inter[1] = dim[1]; bounds_inter[a_dir]++;
  static double R[_MODEL_NVARS_*_MODEL_NVARS_], D[_MODEL_NVARS_*_MODEL_NVARS_],
                L[_MODEL_NVARS_*_MODEL_NVARS_];

  done = 0; int index_outer[2] = {0,0}, index_inter[2];
  while (!done) {
    index_inter[0] = index_outer[0]; index_inter[1] = index_outer[1];
    for (index_inter[a_dir] = 0; index_inter[a_dir] < bounds_inter[a_dir]; index_inter[a_dir]++) {
      int p; _ArrayIndex1D2_(_MODEL_NDIMS_,bounds_inter,index_inter,0,p);
      double uavg[_MODEL_NVARS_], fcL[_MODEL_NVARS_], fcR[_MODEL_NVARS_],
             ucL[_MODEL_NVARS_], ucR[_MODEL_NVARS_], fc[_MODEL_NVARS_];

      /* Roe-Fixed upwinding scheme */

      _Euler2DRoeAverage_(uavg,(a_uL+_MODEL_NVARS_*p),(a_uR+_MODEL_NVARS_*p),param);

      _Euler2DEigenvalues_(uavg,D,param,a_dir);
      _Euler2DLeftEigenvectors_ (uavg,L,param,a_dir);
      _Euler2DRightEigenvectors_(uavg,R,param,a_dir);

      /* calculate characteristic fluxes and variables */
      MatVecMult4(_MODEL_NVARS_,ucL,L,(a_uL+_MODEL_NVARS_*p));
      MatVecMult4(_MODEL_NVARS_,ucR,L,(a_uR+_MODEL_NVARS_*p));
      MatVecMult4(_MODEL_NVARS_,fcL,L,(a_fL+_MODEL_NVARS_*p));
      MatVecMult4(_MODEL_NVARS_,fcR,L,(a_fR+_MODEL_NVARS_*p));

      double eigL[4],eigC[4],eigR[4];
      _Euler2DEigenvalues_((a_uL+_MODEL_NVARS_*p),D,param,a_dir);
      eigL[0] = D[0];
      eigL[1] = D[5];
      eigL[2] = D[10];
      eigL[3] = D[15];
      _Euler2DEigenvalues_((a_uR+_MODEL_NVARS_*p),D,param,a_dir);
      eigR[0] = D[0];
      eigR[1] = D[5];
      eigR[2] = D[10];
      eigR[3] = D[15];
      _Euler2DEigenvalues_(uavg,D,param,a_dir);
      eigC[0] = D[0];
      eigC[1] = D[5];
      eigC[2] = D[10];
      eigC[3] = D[15];

      for (k = 0; k < _MODEL_NVARS_; k++) {
        if ((eigL[k] > 0) && (eigC[k] > 0) && (eigR[k] > 0))      fc[k] = fcL[k];
        else if ((eigL[k] < 0) && (eigC[k] < 0) && (eigR[k] < 0)) fc[k] = fcR[k];
        else {
          double alpha = max3(absolute(eigL[k]),absolute(eigC[k]),absolute(eigR[k]));
          fc[k] = 0.5 * (fcL[k] + fcR[k] + alpha * (ucL[k]-ucR[k]));
        }
      }

      /* calculate the interface flux from the characteristic flux */
      MatVecMult4(_MODEL_NVARS_,(a_fI+_MODEL_NVARS_*p),R,fc);
    }
    _ArrayIncrementIndex_(_MODEL_NDIMS_,bounds_outer,index_outer,done);
  }

  return(0);
}

int Euler2DUpwindLLF(double *a_fI,double *a_fL,double *a_fR,double *a_uL,double *a_uR,double *a_u,int a_dir,void *a_s,double a_t)
{
  HyPar    *solver = (HyPar*)    a_s;
  Euler2D  *param  = (Euler2D*)  solver->m_physics;
  int      done;

  int *dim  = solver->m_dim_local;

  int bounds_outer[2], bounds_inter[2];
  bounds_outer[0] = dim[0]; bounds_outer[1] = dim[1]; bounds_outer[a_dir] = 1;
  bounds_inter[0] = dim[0]; bounds_inter[1] = dim[1]; bounds_inter[a_dir]++;
  static double R[_MODEL_NVARS_*_MODEL_NVARS_], D[_MODEL_NVARS_*_MODEL_NVARS_],
                L[_MODEL_NVARS_*_MODEL_NVARS_];

  done = 0; int index_outer[2] = {0,0}, index_inter[2];
  while (!done) {
    index_inter[0] = index_outer[0]; index_inter[1] = index_outer[1];
    for (index_inter[a_dir] = 0; index_inter[a_dir] < bounds_inter[a_dir]; index_inter[a_dir]++) {
      int p; _ArrayIndex1D2_(_MODEL_NDIMS_,bounds_inter,index_inter,0,p);
      double uavg[_MODEL_NVARS_], fcL[_MODEL_NVARS_], fcR[_MODEL_NVARS_],
             ucL[_MODEL_NVARS_], ucR[_MODEL_NVARS_], fc[_MODEL_NVARS_];

      /* Local Lax-Friedrich upwinding scheme */

      _Euler2DRoeAverage_(uavg,(a_uL+_MODEL_NVARS_*p),(a_uR+_MODEL_NVARS_*p),param);

      _Euler2DEigenvalues_(uavg,D,param,a_dir);
      _Euler2DLeftEigenvectors_ (uavg,L,param,a_dir);
      _Euler2DRightEigenvectors_(uavg,R,param,a_dir);

      /* calculate characteristic fluxes and variables */
      MatVecMult4(_MODEL_NVARS_,ucL,L,(a_uL+_MODEL_NVARS_*p));
      MatVecMult4(_MODEL_NVARS_,ucR,L,(a_uR+_MODEL_NVARS_*p));
      MatVecMult4(_MODEL_NVARS_,fcL,L,(a_fL+_MODEL_NVARS_*p));
      MatVecMult4(_MODEL_NVARS_,fcR,L,(a_fR+_MODEL_NVARS_*p));

      double eigL[4],eigC[4],eigR[4];
      _Euler2DEigenvalues_((a_uL+_MODEL_NVARS_*p),D,param,a_dir);
      eigL[0] = D[0];
      eigL[1] = D[5];
      eigL[2] = D[10];
      eigL[3] = D[15];
      _Euler2DEigenvalues_((a_uR+_MODEL_NVARS_*p),D,param,a_dir);
      eigR[0] = D[0];
      eigR[1] = D[5];
      eigR[2] = D[10];
      eigR[3] = D[15];
      _Euler2DEigenvalues_(uavg,D,param,a_dir);
      eigC[0] = D[0];
      eigC[1] = D[5];
      eigC[2] = D[10];
      eigC[3] = D[15];

      double alpha;
      alpha = max3(absolute(eigL[0]),absolute(eigC[0]),absolute(eigR[0]));
      fc[0] = 0.5 * (fcL[0] + fcR[0] + alpha * (ucL[0]-ucR[0]));
      alpha = max3(absolute(eigL[1]),absolute(eigC[1]),absolute(eigR[1]));
      fc[1] = 0.5 * (fcL[1] + fcR[1] + alpha * (ucL[1]-ucR[1]));
      alpha = max3(absolute(eigL[2]),absolute(eigC[2]),absolute(eigR[2]));
      fc[2] = 0.5 * (fcL[2] + fcR[2] + alpha * (ucL[2]-ucR[2]));
      alpha = max3(absolute(eigL[3]),absolute(eigC[3]),absolute(eigR[3]));
      fc[3] = 0.5 * (fcL[3] + fcR[3] + alpha * (ucL[3]-ucR[3]));

      /* calculate the interface flux from the characteristic flux */
      MatVecMult4(_MODEL_NVARS_,(a_fI+_MODEL_NVARS_*p),R,fc);
    }
    _ArrayIncrementIndex_(_MODEL_NDIMS_,bounds_outer,index_outer,done);
  }

  return(0);
}

int Euler2DUpwindSWFS(double *a_fI,double *a_fL,double *a_fR,double *a_uL,double *a_uR,double *a_u,int a_dir,void *a_s,double a_t)
{
  HyPar     *solver = (HyPar*)    a_s;
  Euler2D   *param  = (Euler2D*)  solver->m_physics;
  int       done,k;
  _DECLARE_IERR_;

  int ndims = solver->m_ndims;
  int *dim  = solver->m_dim_local;

  int index_outer[ndims], index_inter[ndims], bounds_outer[ndims], bounds_inter[ndims];
  _ArrayCopy1D_(dim,bounds_outer,ndims); bounds_outer[a_dir] =  1;
  _ArrayCopy1D_(dim,bounds_inter,ndims); bounds_inter[a_dir] += 1;
  static double fp[_MODEL_NVARS_], fm[_MODEL_NVARS_],uavg[_MODEL_NVARS_];

  done = 0; _ArraySetValue_(index_outer,ndims,0);
  while (!done) {
    _ArrayCopy1D_(index_outer,index_inter,ndims);
    for (index_inter[a_dir] = 0; index_inter[a_dir] < bounds_inter[a_dir]; index_inter[a_dir]++) {
      int p; _ArrayIndex1D_(ndims,bounds_inter,index_inter,0,p);
      double rho,vx,vy,e,P,c,gamma=param->m_gamma,term,Mach,lp[_MODEL_NVARS_],lm[_MODEL_NVARS_];

      /* Steger Warming flux splitting */
      _Euler2DRoeAverage_(uavg,(a_uL+_MODEL_NVARS_*p),(a_uR+_MODEL_NVARS_*p),param);
      _Euler2DGetFlowVar_(uavg,rho,vx,vy,e,P,param);
      Mach = (a_dir==_XDIR_ ? vx : vy) / sqrt(gamma*P/rho);

      if (Mach < -1.0) {

        _ArrayCopy1D_((a_fR+_MODEL_NVARS_*p),(a_fI+_MODEL_NVARS_*p),_MODEL_NVARS_);

      } else if (Mach < 1.0) {

        double kx = 0, ky = 0;
        kx = (a_dir==_XDIR_ ? 1.0 : 0.0);
        ky = (a_dir==_YDIR_ ? 1.0 : 0.0);

        _Euler2DGetFlowVar_((a_uL+_MODEL_NVARS_*p),rho,vx,vy,e,P,param);
        c = sqrt(gamma*P/rho);
        term = rho/(2.0*gamma);
        lp[0] = lp[1] = kx*vx + ky*vy;
        lp[2] = lp[0] + c;
        lp[3] = lp[0] - c;
        for (k=0; k<_MODEL_NVARS_; k++) if (lp[k] < 0.0) lp[k] = 0.0;

        fp[0] = term * (2.0*(gamma-1.0)*lp[0] + lp[2] + lp[3]);
        fp[1] = term * (2.0*(gamma-1.0)*lp[0]*vx + lp[2]*(vx+c*kx) + lp[3]*(vx-c*kx));
        fp[2] = term * (2.0*(gamma-1.0)*lp[0]*vy + lp[2]*(vy+c*ky) + lp[3]*(vy-c*ky));
        fp[3] = term * ((gamma-1.0)*lp[0]*(vx*vx+vy*vy) + 0.5*lp[2]*((vx+c*kx)*(vx+c*kx) + (vy+c*ky)*(vy+c*ky))
                        + 0.5*lp[3]*((vx-c*kx)*(vx-c*kx) + (vy-c*ky)*(vy-c*ky))
                        + ((3.0-gamma)*(lp[2]+lp[3])*c*c)/(2.0*(gamma-1.0)) );

        _Euler2DGetFlowVar_((a_uR+_MODEL_NVARS_*p),rho,vx,vy,e,P,param);
        c = sqrt(gamma*P/rho);
        term = rho/(2.0*gamma);
        lm[0] = lm[1] = kx*vx + ky*vy;
        lm[2] = lm[0] + c;
        lm[3] = lm[0] - c;
        for (k=0; k<_MODEL_NVARS_; k++) if (lm[k] > 0.0) lm[k] = 0.0;

        fm[0] = term * (2.0*(gamma-1.0)*lm[0] + lm[2] + lm[3]);
        fm[1] = term * (2.0*(gamma-1.0)*lm[0]*vx + lm[2]*(vx+c*kx) + lm[3]*(vx-c*kx));
        fm[2] = term * (2.0*(gamma-1.0)*lm[0]*vy + lm[2]*(vy+c*ky) + lm[3]*(vy-c*ky));
        fm[3] = term * ((gamma-1.0)*lm[0]*(vx*vx+vy*vy) + 0.5*lm[2]*((vx+c*kx)*(vx+c*kx) + (vy+c*ky)*(vy+c*ky))
                        + 0.5*lm[3]*((vx-c*kx)*(vx-c*kx) + (vy-c*ky)*(vy-c*ky))
                        + ((3.0-gamma)*(lm[2]+lm[3])*c*c)/(2.0*(gamma-1.0)) );

        _ArrayAdd1D_((a_fI+_MODEL_NVARS_*p),fp,fm,_MODEL_NVARS_);

      } else {

        _ArrayCopy1D_((a_fL+_MODEL_NVARS_*p),(a_fI+_MODEL_NVARS_*p),_MODEL_NVARS_);

      }

    }
    _ArrayIncrementIndex_(ndims,bounds_outer,index_outer,done);
  }

  return(0);
}

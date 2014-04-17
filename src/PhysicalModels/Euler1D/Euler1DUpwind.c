#include <stdlib.h>
#include <math.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mathfunctions.h>
#include <physicalmodels/euler1d.h>
#include <hypar.h>

int Euler1DUpwindRoe(double *fI,double *fL,double *fR,double *uL,double *uR,double *u,int dir,void *s,double t)
{
  HyPar     *solver = (HyPar*)    s;
  Euler1D   *param  = (Euler1D*)  solver->physics;
  int       done,k;
  _DECLARE_IERR_;

  int ndims = solver->ndims;
  int *dim  = solver->dim_local;

  int index_outer[ndims], index_inter[ndims], bounds_outer[ndims], bounds_inter[ndims];
  _ArrayCopy1D_(dim,bounds_outer,ndims); bounds_outer[dir] =  1;
  _ArrayCopy1D_(dim,bounds_inter,ndims); bounds_inter[dir] += 1;
  static double R[_MODEL_NVARS_*_MODEL_NVARS_], D[_MODEL_NVARS_*_MODEL_NVARS_], L[_MODEL_NVARS_*_MODEL_NVARS_], 
                DL[_MODEL_NVARS_*_MODEL_NVARS_], modA[_MODEL_NVARS_*_MODEL_NVARS_];

  done = 0; _ArraySetValue_(index_outer,ndims,0);
  while (!done) {
    _ArrayCopy1D_(index_outer,index_inter,ndims);
    for (index_inter[dir] = 0; index_inter[dir] < bounds_inter[dir]; index_inter[dir]++) {
      int p; _ArrayIndex1D_(ndims,bounds_inter,index_inter,0,p);
      double udiff[_MODEL_NVARS_], uavg[_MODEL_NVARS_],udiss[_MODEL_NVARS_];

      /* Roe's upwinding scheme */

      udiff[0] = 0.5 * (uR[_MODEL_NVARS_*p+0] - uL[_MODEL_NVARS_*p+0]);
      udiff[1] = 0.5 * (uR[_MODEL_NVARS_*p+1] - uL[_MODEL_NVARS_*p+1]);
      udiff[2] = 0.5 * (uR[_MODEL_NVARS_*p+2] - uL[_MODEL_NVARS_*p+2]);

      _Euler1DRoeAverage_(uavg,(uL+_MODEL_NVARS_*p),(uR+_MODEL_NVARS_*p),param); 

      _Euler1DEigenvalues_(uavg,D,param,0);
      _Euler1DLeftEigenvectors_(uavg,L,param,0);
      _Euler1DRightEigenvectors_ (uavg,R,param,0);

      for (k=0; k<_MODEL_NVARS_; k++) D[k*_MODEL_NVARS_+k] = absolute(D[k*_MODEL_NVARS_+k]);
      MatMult3(3,DL,D,L);
      MatMult3(3,modA,R,DL);

      udiss[0] = modA[0*_MODEL_NVARS_+0]*udiff[0] + modA[0*_MODEL_NVARS_+1]*udiff[1] + modA[0*_MODEL_NVARS_+2]*udiff[2];
      udiss[1] = modA[1*_MODEL_NVARS_+0]*udiff[0] + modA[1*_MODEL_NVARS_+1]*udiff[1] + modA[1*_MODEL_NVARS_+2]*udiff[2];
      udiss[2] = modA[2*_MODEL_NVARS_+0]*udiff[0] + modA[2*_MODEL_NVARS_+1]*udiff[1] + modA[2*_MODEL_NVARS_+2]*udiff[2];

      fI[_MODEL_NVARS_*p+0] = 0.5 * (fL[_MODEL_NVARS_*p+0]+fR[_MODEL_NVARS_*p+0]) - udiss[0];
      fI[_MODEL_NVARS_*p+1] = 0.5 * (fL[_MODEL_NVARS_*p+1]+fR[_MODEL_NVARS_*p+1]) - udiss[1];
      fI[_MODEL_NVARS_*p+2] = 0.5 * (fL[_MODEL_NVARS_*p+2]+fR[_MODEL_NVARS_*p+2]) - udiss[2];
    }
    _ArrayIncrementIndex_(ndims,bounds_outer,index_outer,done);
  }

  return(0);
}

int Euler1DUpwindRF(double *fI,double *fL,double *fR,double *uL,double *uR,double *u,int dir,void *s,double t)
{
  HyPar     *solver = (HyPar*)    s;
  Euler1D   *param  = (Euler1D*)  solver->physics;
  int       done,k;
  _DECLARE_IERR_;

  int ndims = solver->ndims;
  int *dim  = solver->dim_local;

  int index_outer[ndims], index_inter[ndims], bounds_outer[ndims], bounds_inter[ndims];
  _ArrayCopy1D_(dim,bounds_outer,ndims); bounds_outer[dir] =  1;
  _ArrayCopy1D_(dim,bounds_inter,ndims); bounds_inter[dir] += 1;
  static double R[_MODEL_NVARS_*_MODEL_NVARS_], D[_MODEL_NVARS_*_MODEL_NVARS_], L[_MODEL_NVARS_*_MODEL_NVARS_];

  done = 0; _ArraySetValue_(index_outer,ndims,0);
  while (!done) {
    _ArrayCopy1D_(index_outer,index_inter,ndims);
    for (index_inter[dir] = 0; index_inter[dir] < bounds_inter[dir]; index_inter[dir]++) {
      int p; _ArrayIndex1D_(ndims,bounds_inter,index_inter,0,p);
      double uavg[_MODEL_NVARS_], fcL[_MODEL_NVARS_], fcR[_MODEL_NVARS_], 
             ucL[_MODEL_NVARS_], ucR[_MODEL_NVARS_], fc[_MODEL_NVARS_];

      /* Roe-Fixed upwinding scheme */

      _Euler1DRoeAverage_(uavg,(uL+_MODEL_NVARS_*p),(uR+_MODEL_NVARS_*p),param);

      _Euler1DEigenvalues_(uavg,D,param,0);
      _Euler1DLeftEigenvectors_(uavg,L,param,0);
      _Euler1DRightEigenvectors_(uavg,R,param,0);

      /* calculate characteristic fluxes and variables */
      MatVecMult3(_MODEL_NVARS_,ucL,L,(uL+_MODEL_NVARS_*p));
      MatVecMult3(_MODEL_NVARS_,ucR,L,(uR+_MODEL_NVARS_*p));
      MatVecMult3(_MODEL_NVARS_,fcL,L,(fL+_MODEL_NVARS_*p));
      MatVecMult3(_MODEL_NVARS_,fcR,L,(fR+_MODEL_NVARS_*p));

      for (k = 0; k < _MODEL_NVARS_; k++) {
        double eigL,eigC,eigR;
        _Euler1DEigenvalues_((uL+_MODEL_NVARS_*p),D,param,0);
        eigL = D[k*_MODEL_NVARS_+k];
        _Euler1DEigenvalues_((uR+_MODEL_NVARS_*p),D,param,0);
        eigR = D[k*_MODEL_NVARS_+k];
        _Euler1DEigenvalues_(uavg,D,param,0);
        eigC = D[k*_MODEL_NVARS_+k];

        if ((eigL > 0) && (eigC > 0) && (eigR > 0)) {
          fc[k] = fcL[k];
        } else if ((eigL < 0) && (eigC < 0) && (eigR < 0)) {
          fc[k] = fcR[k];
        } else {
          double alpha = max3(absolute(eigL),absolute(eigC),absolute(eigR));
          fc[k] = 0.5 * (fcL[k] + fcR[k] + alpha * (ucL[k]-ucR[k]));
        }
      }

      /* calculate the interface flux from the characteristic flux */
      MatVecMult3(_MODEL_NVARS_,(fI+_MODEL_NVARS_*p),R,fc);
    }
    _ArrayIncrementIndex_(ndims,bounds_outer,index_outer,done);
  }

  return(0);
}

int Euler1DUpwindLLF(double *fI,double *fL,double *fR,double *uL,double *uR,double *u,int dir,void *s,double t)
{
  HyPar     *solver = (HyPar*)    s;
  Euler1D   *param  = (Euler1D*)  solver->physics;
  int       done,k;
  _DECLARE_IERR_;

  int ndims = solver->ndims;
  int *dim  = solver->dim_local;

  int index_outer[ndims], index_inter[ndims], bounds_outer[ndims], bounds_inter[ndims];
  _ArrayCopy1D_(dim,bounds_outer,ndims); bounds_outer[dir] =  1;
  _ArrayCopy1D_(dim,bounds_inter,ndims); bounds_inter[dir] += 1;
  static double R[_MODEL_NVARS_*_MODEL_NVARS_], D[_MODEL_NVARS_*_MODEL_NVARS_], L[_MODEL_NVARS_*_MODEL_NVARS_];

  done = 0; _ArraySetValue_(index_outer,ndims,0);
  while (!done) {
    _ArrayCopy1D_(index_outer,index_inter,ndims);
    for (index_inter[dir] = 0; index_inter[dir] < bounds_inter[dir]; index_inter[dir]++) {
      int p; _ArrayIndex1D_(ndims,bounds_inter,index_inter,0,p);
      double uavg[_MODEL_NVARS_], fcL[_MODEL_NVARS_], fcR[_MODEL_NVARS_], 
             ucL[_MODEL_NVARS_], ucR[_MODEL_NVARS_], fc[_MODEL_NVARS_];

      /* Local Lax-Friedrich upwinding scheme */

      _Euler1DRoeAverage_(uavg,(uL+_MODEL_NVARS_*p),(uR+_MODEL_NVARS_*p),param);

      _Euler1DEigenvalues_(uavg,D,param,0);
      _Euler1DLeftEigenvectors_(uavg,L,param,0);
      _Euler1DRightEigenvectors_(uavg,R,param,0);

      /* calculate characteristic fluxes and variables */
      MatVecMult3(_MODEL_NVARS_,ucL,L,(uL+_MODEL_NVARS_*p));
      MatVecMult3(_MODEL_NVARS_,ucR,L,(uR+_MODEL_NVARS_*p));
      MatVecMult3(_MODEL_NVARS_,fcL,L,(fL+_MODEL_NVARS_*p));
      MatVecMult3(_MODEL_NVARS_,fcR,L,(fR+_MODEL_NVARS_*p));

      for (k = 0; k < _MODEL_NVARS_; k++) {
        double eigL,eigC,eigR;
        _Euler1DEigenvalues_((uL+_MODEL_NVARS_*p),D,param,0);
        eigL = absolute(D[k*_MODEL_NVARS_+k]);
        _Euler1DEigenvalues_((uR+_MODEL_NVARS_*p),D,param,0);
        eigR = absolute(D[k*_MODEL_NVARS_+k]);
        _Euler1DEigenvalues_(uavg,D,param,0);
        eigC = absolute(D[k*_MODEL_NVARS_+k]);

        double alpha = max3(eigL,eigC,eigR);
        fc[k] = 0.5 * (fcL[k] + fcR[k] + alpha * (ucL[k]-ucR[k]));
      }

      /* calculate the interface flux from the characteristic flux */
      MatVecMult3(_MODEL_NVARS_,(fI+_MODEL_NVARS_*p),R,fc);
    }
    _ArrayIncrementIndex_(ndims,bounds_outer,index_outer,done);
  }

  return(0);
}

int Euler1DUpwindSWFS(double *fI,double *fL,double *fR,double *uL,double *uR,double *u,int dir,void *s,double t)
{
  HyPar     *solver = (HyPar*)    s;
  Euler1D   *param  = (Euler1D*)  solver->physics;
  int       done,k;
  _DECLARE_IERR_;

  int ndims = solver->ndims;
  int *dim  = solver->dim_local;

  int index_outer[ndims], index_inter[ndims], bounds_outer[ndims], bounds_inter[ndims];
  _ArrayCopy1D_(dim,bounds_outer,ndims); bounds_outer[dir] =  1;
  _ArrayCopy1D_(dim,bounds_inter,ndims); bounds_inter[dir] += 1;
  static double fp[_MODEL_NVARS_], fm[_MODEL_NVARS_],uavg[_MODEL_NVARS_];

  done = 0; _ArraySetValue_(index_outer,ndims,0);
  while (!done) {
    _ArrayCopy1D_(index_outer,index_inter,ndims);
    for (index_inter[dir] = 0; index_inter[dir] < bounds_inter[dir]; index_inter[dir]++) {
      int p; _ArrayIndex1D_(ndims,bounds_inter,index_inter,0,p);
      double rho,v,e,P,c,gamma=param->gamma,term,Mach;

      /* Steger Warming flux splitting */
      _Euler1DRoeAverage_(uavg,(uL+_MODEL_NVARS_*p),(uR+_MODEL_NVARS_*p),param); 
      _Euler1DGetFlowVar_(uavg,rho,v,e,P,param);
      Mach = v/sqrt(gamma*P/rho);

      if (Mach < -1.0) {

        _ArrayCopy1D3_((fR+_MODEL_NVARS_*p),(fI+_MODEL_NVARS_*p),_MODEL_NVARS_);

      } else if (Mach < 1.0) {

        _Euler1DGetFlowVar_((uL+_MODEL_NVARS_*p),rho,v,e,P,param);
        c = sqrt(gamma*P/rho);
        term = rho/(2.0*gamma);

        fp[0] = term * (2*gamma*v + c - v);
        fp[1] = term * (2*(gamma-1.0)*v*v + (v+c)*(v+c));
        fp[2] = term * ((gamma-1.0)*v*v*v + 0.5*(v+c)*(v+c)*(v+c) + ((3.0-gamma)*(v+c)*c*c)/(2.0*(gamma-1.0)));

        _Euler1DGetFlowVar_((uR+_MODEL_NVARS_*p),rho,v,e,P,param);
        c = sqrt(gamma*P/rho);
        term = rho/(2.0*gamma);

        fm[0] = term * (v - c);
        fm[1] = term * (v-c) * (v-c);
        fm[2] = term * (0.5*(v-c)*(v-c)*(v-c) + ((3.0-gamma)*(v-c)*c*c)/(2.0*(gamma-1.0)));

        _ArrayAdd1D_((fI+_MODEL_NVARS_*p),fp,fm,_MODEL_NVARS_);

      } else {

        _ArrayCopy1D3_((fL+_MODEL_NVARS_*p),(fI+_MODEL_NVARS_*p),_MODEL_NVARS_);

      }

    }
    _ArrayIncrementIndex_(ndims,bounds_outer,index_outer,done);
  }

  return(0);
}

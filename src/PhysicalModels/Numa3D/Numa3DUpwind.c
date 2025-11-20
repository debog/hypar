#include <stdlib.h>
#include <math.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mathfunctions.h>
#include <physicalmodels/numa3d.h>
#include <hypar.h>

int Numa3DRusanovFlux(double *a_fI,double *a_fL,double *a_fR,double *a_uL,double *a_uR,double *a_u,int a_dir,void *a_s,double a_t)
{
  HyPar   *solver = (HyPar*)  a_s;
  Numa3D  *param  = (Numa3D*) solver->m_physics;
  int      done;

  int *dim   = solver->m_dim_local;
  int ghosts = solver->m_ghosts;

  int bounds_outer[_MODEL_NDIMS_], bounds_inter[_MODEL_NDIMS_];
  _ArrayCopy1D3_(dim,bounds_outer,_MODEL_NDIMS_); bounds_outer[a_dir] =  1;
  _ArrayCopy1D3_(dim,bounds_inter,_MODEL_NDIMS_); bounds_inter[a_dir] += 1;

  done = 0; int index_outer[3] = {0,0,0}, index_inter[3];
  while (!done) {
    _ArrayCopy1D3_(index_outer,index_inter,_MODEL_NDIMS_);
    for (index_inter[a_dir] = 0; index_inter[a_dir] < bounds_inter[a_dir]; index_inter[a_dir]++) {
      int p; _ArrayIndex1D3_(_MODEL_NDIMS_,bounds_inter,index_inter,0,p);
      double udiff[_MODEL_NVARS_],drho,vel[3],dT,rho0,P0,T0,EP,c;
      double zcoordL, zcoordR;

      if (a_dir == _ZDIR_) {
        _GetCoordinate_(_ZDIR_,(index_inter[_ZDIR_]-1),dim,ghosts,solver->m_x,zcoordL);
        _GetCoordinate_(_ZDIR_,(index_inter[_ZDIR_]  ),dim,ghosts,solver->m_x,zcoordR);
      } else {
        _GetCoordinate_(_ZDIR_,(index_inter[_ZDIR_]  ),dim,ghosts,solver->m_x,zcoordL);
        zcoordR = zcoordL;
      }

      /* Rusanov'a_s upwinding scheme */

      udiff[0] = 0.5 * (a_uR[_MODEL_NVARS_*p+0] - a_uL[_MODEL_NVARS_*p+0]);
      udiff[1] = 0.5 * (a_uR[_MODEL_NVARS_*p+1] - a_uL[_MODEL_NVARS_*p+1]);
      udiff[2] = 0.5 * (a_uR[_MODEL_NVARS_*p+2] - a_uL[_MODEL_NVARS_*p+2]);
      udiff[3] = 0.5 * (a_uR[_MODEL_NVARS_*p+3] - a_uL[_MODEL_NVARS_*p+3]);
      udiff[4] = 0.5 * (a_uR[_MODEL_NVARS_*p+4] - a_uL[_MODEL_NVARS_*p+4]);

      /* left of the interface */
      param->StandardAtmosphere(param,zcoordL,&EP,&P0,&rho0,&T0);
      _Numa3DGetFlowVars_         ((a_uL+_MODEL_NVARS_*p),drho,vel[0],vel[1],vel[2],dT,rho0);
      _Numa3DComputeSpeedofSound_ (param->m_gamma,param->m_R,T0,dT,rho0,drho,EP,c);
      double alphaL = c + absolute(vel[a_dir]);

      /* right of the interface */
      param->StandardAtmosphere(param,zcoordR,&EP,&P0,&rho0,&T0);
      _Numa3DGetFlowVars_         ((a_uR+_MODEL_NVARS_*p),drho,vel[0],vel[1],vel[2],dT,rho0);
      _Numa3DComputeSpeedofSound_ (param->m_gamma,param->m_R,T0,dT,rho0,drho,EP,c);
      double alphaR = c + absolute(vel[a_dir]);

      double alpha = max(alphaL,alphaR);
      a_fI[_MODEL_NVARS_*p+0] = 0.5*(a_fL[_MODEL_NVARS_*p+0]+a_fR[_MODEL_NVARS_*p+0])-alpha*udiff[0];
      a_fI[_MODEL_NVARS_*p+1] = 0.5*(a_fL[_MODEL_NVARS_*p+1]+a_fR[_MODEL_NVARS_*p+1])-alpha*udiff[1];
      a_fI[_MODEL_NVARS_*p+2] = 0.5*(a_fL[_MODEL_NVARS_*p+2]+a_fR[_MODEL_NVARS_*p+2])-alpha*udiff[2];
      a_fI[_MODEL_NVARS_*p+3] = 0.5*(a_fL[_MODEL_NVARS_*p+3]+a_fR[_MODEL_NVARS_*p+3])-alpha*udiff[3];
      a_fI[_MODEL_NVARS_*p+4] = 0.5*(a_fL[_MODEL_NVARS_*p+4]+a_fR[_MODEL_NVARS_*p+4])-alpha*udiff[4];
    }
    _ArrayIncrementIndex_(_MODEL_NDIMS_,bounds_outer,index_outer,done);
  }

  return(0);
}

int Numa3DRusanovLinearFlux(double *a_fI,double *a_fL,double *a_fR,double *a_uL,double *a_uR,double *a_u,int a_dir,void *a_s,double a_t)
{
  HyPar   *solver = (HyPar*)  a_s;
  Numa3D  *param  = (Numa3D*) solver->m_physics;
  int      done;

  int *dim   = solver->m_dim_local;
  int ghosts = solver->m_ghosts;

  int bounds_outer[_MODEL_NDIMS_], bounds_inter[_MODEL_NDIMS_];
  _ArrayCopy1D3_(dim,bounds_outer,_MODEL_NDIMS_); bounds_outer[a_dir] =  1;
  _ArrayCopy1D3_(dim,bounds_inter,_MODEL_NDIMS_); bounds_inter[a_dir] += 1;

  done = 0; int index_outer[3] = {0,0,0}, index_inter[3];
  while (!done) {
    _ArrayCopy1D3_(index_outer,index_inter,_MODEL_NDIMS_);
    for (index_inter[a_dir] = 0; index_inter[a_dir] < bounds_inter[a_dir]; index_inter[a_dir]++) {
      int p; _ArrayIndex1D3_(_MODEL_NDIMS_,bounds_inter,index_inter,0,p);
      double udiff[_MODEL_NVARS_],drho,vel[3],dT,rho0,P0,T0,EP,c;
      double zcoordL, zcoordR;

      if (a_dir == _ZDIR_) {
        _GetCoordinate_(_ZDIR_,(index_inter[_ZDIR_]-1),dim,ghosts,solver->m_x,zcoordL);
        _GetCoordinate_(_ZDIR_,(index_inter[_ZDIR_]  ),dim,ghosts,solver->m_x,zcoordR);
      } else {
        _GetCoordinate_(_ZDIR_,(index_inter[_ZDIR_]  ),dim,ghosts,solver->m_x,zcoordL);
        zcoordR = zcoordL;
      }

      /* Rusanov'a_s upwinding scheme */

      udiff[0] = 0.5 * (a_uR[_MODEL_NVARS_*p+0] - a_uL[_MODEL_NVARS_*p+0]);
      udiff[1] = 0.5 * (a_uR[_MODEL_NVARS_*p+1] - a_uL[_MODEL_NVARS_*p+1]);
      udiff[2] = 0.5 * (a_uR[_MODEL_NVARS_*p+2] - a_uL[_MODEL_NVARS_*p+2]);
      udiff[3] = 0.5 * (a_uR[_MODEL_NVARS_*p+3] - a_uL[_MODEL_NVARS_*p+3]);
      udiff[4] = 0.5 * (a_uR[_MODEL_NVARS_*p+4] - a_uL[_MODEL_NVARS_*p+4]);

      /* left of the interface */
      param->StandardAtmosphere   (param,zcoordL,&EP,&P0,&rho0,&T0);
      _Numa3DGetFlowVars_         ((a_uL+_MODEL_NVARS_*p),drho,vel[0],vel[1],vel[2],dT,rho0);
      _Numa3DComputeLinearizedSpeedofSound_ (param->m_gamma,param->m_R,T0,rho0,EP,c);
      double alphaL = c + absolute(vel[a_dir]);

      /* right of the interface */
      param->StandardAtmosphere(param,zcoordR,&EP,&P0,&rho0,&T0);
      _Numa3DGetFlowVars_         ((a_uR+_MODEL_NVARS_*p),drho,vel[0],vel[1],vel[2],dT,rho0);
      _Numa3DComputeLinearizedSpeedofSound_ (param->m_gamma,param->m_R,T0,rho0,EP,c);
      double alphaR = c + absolute(vel[a_dir]);

      double alpha = max(alphaL,alphaR);
      a_fI[_MODEL_NVARS_*p+0] = 0.5*(a_fL[_MODEL_NVARS_*p+0]+a_fR[_MODEL_NVARS_*p+0])-alpha*udiff[0];
      a_fI[_MODEL_NVARS_*p+1] = 0.5*(a_fL[_MODEL_NVARS_*p+1]+a_fR[_MODEL_NVARS_*p+1])-alpha*udiff[1];
      a_fI[_MODEL_NVARS_*p+2] = 0.5*(a_fL[_MODEL_NVARS_*p+2]+a_fR[_MODEL_NVARS_*p+2])-alpha*udiff[2];
      a_fI[_MODEL_NVARS_*p+3] = 0.5*(a_fL[_MODEL_NVARS_*p+3]+a_fR[_MODEL_NVARS_*p+3])-alpha*udiff[3];
      a_fI[_MODEL_NVARS_*p+4] = 0.5*(a_fL[_MODEL_NVARS_*p+4]+a_fR[_MODEL_NVARS_*p+4])-alpha*udiff[4];

      vel[0] = dT; /* useless statement to avoid compiler warning */
    }
    _ArrayIncrementIndex_(_MODEL_NDIMS_,bounds_outer,index_outer,done);
  }

  return(0);
}

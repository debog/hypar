#include <stdlib.h>
#include <math.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mathfunctions.h>
#include <physicalmodels/numa2d.h>
#include <hypar.h>

int Numa2DRusanovFlux(double *a_fI,double *a_fL,double *a_fR,double *a_uL,double *a_uR,double *a_u,int a_dir,void *a_s,double a_t)
{
  HyPar   *solver = (HyPar*)  a_s;
  Numa2D  *param  = (Numa2D*) solver->m_physics;
  int      done;

  int *dim   = solver->m_dim_local;
  int ghosts = solver->m_ghosts;

  int bounds_outer[_MODEL_NDIMS_], bounds_inter[_MODEL_NDIMS_];
  _ArrayCopy1D2_(dim,bounds_outer,_MODEL_NDIMS_); bounds_outer[a_dir] =  1;
  _ArrayCopy1D2_(dim,bounds_inter,_MODEL_NDIMS_); bounds_inter[a_dir] += 1;

  done = 0; int index_outer[2] = {0,0}, index_inter[2];
  while (!done) {
    _ArrayCopy1D2_(index_outer,index_inter,_MODEL_NDIMS_);
    for (index_inter[a_dir] = 0; index_inter[a_dir] < bounds_inter[a_dir]; index_inter[a_dir]++) {
      int p; _ArrayIndex1D2_(_MODEL_NDIMS_,bounds_inter,index_inter,0,p);
      double udiff[_MODEL_NVARS_],drho,vel[2],dT,rho0,P0,T0,EP,c;
      double ycoordL, ycoordR;

      if (a_dir == _YDIR_) {
        _GetCoordinate_(_YDIR_,(index_inter[_YDIR_]-1),dim,ghosts,solver->m_x,ycoordL);
        _GetCoordinate_(_YDIR_,(index_inter[_YDIR_]  ),dim,ghosts,solver->m_x,ycoordR);
      } else {
        _GetCoordinate_(_YDIR_,(index_inter[_YDIR_]  ),dim,ghosts,solver->m_x,ycoordL);
        ycoordR = ycoordL;
      }

      /* Rusanov'a_s upwinding scheme */

      udiff[0] = 0.5 * (a_uR[_MODEL_NVARS_*p+0] - a_uL[_MODEL_NVARS_*p+0]);
      udiff[1] = 0.5 * (a_uR[_MODEL_NVARS_*p+1] - a_uL[_MODEL_NVARS_*p+1]);
      udiff[2] = 0.5 * (a_uR[_MODEL_NVARS_*p+2] - a_uL[_MODEL_NVARS_*p+2]);
      udiff[3] = 0.5 * (a_uR[_MODEL_NVARS_*p+3] - a_uL[_MODEL_NVARS_*p+3]);

      /* left of the interface */
      param->StandardAtmosphere   (param,ycoordL,&EP,&P0,&rho0,&T0);
      _Numa2DGetFlowVars_         ((a_uL+_MODEL_NVARS_*p),drho,vel[0],vel[1],dT,rho0);
      _Numa2DComputeSpeedofSound_ (param->m_gamma,param->m_R,T0,dT,rho0,drho,EP,c);
      double alphaL = c + absolute(vel[a_dir]);

      /* right of the interface */
      param->StandardAtmosphere   (param,ycoordR,&EP,&P0,&rho0,&T0);
      _Numa2DGetFlowVars_         ((a_uR+_MODEL_NVARS_*p),drho,vel[0],vel[1],dT,rho0);
      _Numa2DComputeSpeedofSound_ (param->m_gamma,param->m_R,T0,dT,rho0,drho,EP,c);
      double alphaR = c + absolute(vel[a_dir]);

      double alpha = max(alphaL,alphaR);
      a_fI[_MODEL_NVARS_*p+0] = 0.5*(a_fL[_MODEL_NVARS_*p+0]+a_fR[_MODEL_NVARS_*p+0])-alpha*udiff[0];
      a_fI[_MODEL_NVARS_*p+1] = 0.5*(a_fL[_MODEL_NVARS_*p+1]+a_fR[_MODEL_NVARS_*p+1])-alpha*udiff[1];
      a_fI[_MODEL_NVARS_*p+2] = 0.5*(a_fL[_MODEL_NVARS_*p+2]+a_fR[_MODEL_NVARS_*p+2])-alpha*udiff[2];
      a_fI[_MODEL_NVARS_*p+3] = 0.5*(a_fL[_MODEL_NVARS_*p+3]+a_fR[_MODEL_NVARS_*p+3])-alpha*udiff[3];
    }
    _ArrayIncrementIndex_(_MODEL_NDIMS_,bounds_outer,index_outer,done);
  }

  return(0);
}

int Numa2DRusanovLinearFlux(double *a_fI,double *a_fL,double *a_fR,double *a_uL,double *a_uR,double *a_u,int a_dir,void *a_s,double a_t)
{
  HyPar   *solver = (HyPar*)  a_s;
  Numa2D  *param  = (Numa2D*) solver->m_physics;
  int      done;

  int *dim   = solver->m_dim_local;
  int ghosts = solver->m_ghosts;

  int bounds_outer[_MODEL_NDIMS_], bounds_inter[_MODEL_NDIMS_];
  _ArrayCopy1D2_(dim,bounds_outer,_MODEL_NDIMS_); bounds_outer[a_dir] =  1;
  _ArrayCopy1D2_(dim,bounds_inter,_MODEL_NDIMS_); bounds_inter[a_dir] += 1;

  done = 0; int index_outer[2] = {0,0}, index_inter[2];
  while (!done) {
    _ArrayCopy1D2_(index_outer,index_inter,_MODEL_NDIMS_);
    for (index_inter[a_dir] = 0; index_inter[a_dir] < bounds_inter[a_dir]; index_inter[a_dir]++) {
      int p; _ArrayIndex1D2_(_MODEL_NDIMS_,bounds_inter,index_inter,0,p);
      double udiff[_MODEL_NVARS_],drho,vel[2],dT,rho0,P0,T0,EP,c;
      double ycoordL, ycoordR;

      if (a_dir == _YDIR_) {
        _GetCoordinate_(_YDIR_,(index_inter[_YDIR_]-1),dim,ghosts,solver->m_x,ycoordL);
        _GetCoordinate_(_YDIR_,(index_inter[_YDIR_]  ),dim,ghosts,solver->m_x,ycoordR);
      } else {
        _GetCoordinate_(_YDIR_,(index_inter[_YDIR_]  ),dim,ghosts,solver->m_x,ycoordL);
        ycoordR = ycoordL;
      }

      /* Rusanov'a_s upwinding scheme */

      udiff[0] = 0.5 * (a_uR[_MODEL_NVARS_*p+0] - a_uL[_MODEL_NVARS_*p+0]);
      udiff[1] = 0.5 * (a_uR[_MODEL_NVARS_*p+1] - a_uL[_MODEL_NVARS_*p+1]);
      udiff[2] = 0.5 * (a_uR[_MODEL_NVARS_*p+2] - a_uL[_MODEL_NVARS_*p+2]);
      udiff[3] = 0.5 * (a_uR[_MODEL_NVARS_*p+3] - a_uL[_MODEL_NVARS_*p+3]);

      /* left of the interface */
      param->StandardAtmosphere             (param,ycoordL,&EP,&P0,&rho0,&T0);
      _Numa2DGetFlowVars_                   ((a_uL+_MODEL_NVARS_*p),drho,vel[0],vel[1],dT,rho0);
      _Numa2DComputeLinearizedSpeedofSound_ (param->m_gamma,param->m_R,T0,rho0,EP,c);
      double alphaL = c;

      /* right of the interface */
      param->StandardAtmosphere             (param,ycoordR,&EP,&P0,&rho0,&T0);
      _Numa2DGetFlowVars_                   ((a_uR+_MODEL_NVARS_*p),drho,vel[0],vel[1],dT,rho0);
      _Numa2DComputeLinearizedSpeedofSound_ (param->m_gamma,param->m_R,T0,rho0,EP,c);
      double alphaR = c;

      double alpha = max(alphaL,alphaR);
      a_fI[_MODEL_NVARS_*p+0] = 0.5*(a_fL[_MODEL_NVARS_*p+0]+a_fR[_MODEL_NVARS_*p+0])-alpha*udiff[0];
      a_fI[_MODEL_NVARS_*p+1] = 0.5*(a_fL[_MODEL_NVARS_*p+1]+a_fR[_MODEL_NVARS_*p+1])-alpha*udiff[1];
      a_fI[_MODEL_NVARS_*p+2] = 0.5*(a_fL[_MODEL_NVARS_*p+2]+a_fR[_MODEL_NVARS_*p+2])-alpha*udiff[2];
      a_fI[_MODEL_NVARS_*p+3] = 0.5*(a_fL[_MODEL_NVARS_*p+3]+a_fR[_MODEL_NVARS_*p+3])-alpha*udiff[3];

      /* some harmless statements to avoid compiler warnings */
      vel[0] = vel[1];
      dT = vel[0];
      vel[1] = dT;
    }
    _ArrayIncrementIndex_(_MODEL_NDIMS_,bounds_outer,index_outer,done);
  }

  return(0);
}


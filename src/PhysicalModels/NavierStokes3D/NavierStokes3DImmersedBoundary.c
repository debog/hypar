/*! @file NavierStokes3DImmersedBoundary.c
    @brief Immersed boundary treatment for 3D Navier-Stokes equations
    @author Debojyoti Ghosh
*/

#include <math.h>
#include <string.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <immersedboundaries.h>
#include <physicalmodels/navierstokes3d.h>
#include <mpivars.h>
#include <hypar.h>

/*! Apply no-slip adiabatic wall boundary conditions on the immersed boundary
    points (grid points within the immersed body that are within
    stencil-width distance of interior points, i.e., points in the
    interior of the computational domain). */
int NavierStokes3DIBAdiabatic(void    *s, /*!< Solver object of type #HyPar */
                              void    *m, /*!< Solver object of type #HyPar */
                              double  *u, /*!< Array with the solution vector */
                              double  t   /*!< Current simulation time */
                             )
{
  HyPar             *solver   = (HyPar*)   s;
  MPIVariables      *mpi      = (MPIVariables*) m;
  ImmersedBoundary  *IB       = (ImmersedBoundary*) solver->ib;
  IBNode            *boundary = IB->boundary;
  NavierStokes3D    *param    = (NavierStokes3D*) solver->physics;
  static double     v[_MODEL_NVARS_];
  int               n, j, k, nb = IB->n_boundary_nodes;

  if (!solver->flag_ib) return(0);

  /* Ideally, this shouldn't be here - But this function is called everywhere
     (through ApplyIBConditions()) *before* MPIExchangeBoundariesnD is called! */
  MPIExchangeBoundariesnD(_MODEL_NDIMS_,_MODEL_NVARS_,solver->dim_local,solver->ghosts,mpi,u);

  double inv_gamma_m1 = 1.0 / (param->gamma - 1.0);

  double ramp_fac = 1.0;
  if (param->t_ib_ramp > 0) {
    double x = t/param->t_ib_ramp;
    if (!strcmp(param->ib_ramp_type,_IB_RAMP_LINEAR_)) {
      ramp_fac = x;
    } else if (!strcmp(param->ib_ramp_type,_IB_RAMP_SMOOTHEDSLAB_)) {
      double a = 0.0;
      double b = 1.0;
      double c = 0.5;
      double r = param->t_ib_ramp/param->t_ib_width;
      ramp_fac = (a*exp(c*r)+b*exp(r*x))/(exp(c*r)+exp(r*x));
    } else if (!strcmp(param->ib_ramp_type,_IB_RAMP_DISABLE_)) {
      ramp_fac = 0.0;
    } else {
      fprintf(stderr,"Error in NavierStokes3DImmersedBoundary():\n");
      fprintf(stderr,"  Ramp type %s not recognized.\n", param->ib_ramp_type);
      return 1;
    }
  }

  for (n=0; n<nb; n++) {

    int     node_index = boundary[n].p;
    double  *alpha = &(boundary[n].interp_coeffs[0]);
    int     *nodes = &(boundary[n].interp_nodes[0]);
    double  factor = boundary[n].surface_distance / boundary[n].interp_node_distance;

    _ArraySetValue_(v,_MODEL_NVARS_,0.0);
    for (j=0; j<_IB_NNODES_; j++) {
      for (k=0; k<_MODEL_NVARS_; k++) {
        v[k] += ( alpha[j] * u[_MODEL_NVARS_*nodes[j]+k] );
      }
    }

    double rho, uvel, vvel, wvel, energy, pressure;
    _NavierStokes3DGetFlowVar_(v,_NavierStokes3D_stride_,rho,uvel,vvel,wvel,energy,pressure,param->gamma);

    double rho_gpt, uvel_gpt, vvel_gpt, wvel_gpt, energy_gpt, pressure_gpt;
    _NavierStokes3DGetFlowVar_( (u+_MODEL_NVARS_*node_index),
                                _NavierStokes3D_stride_,
                                rho_gpt,
                                uvel_gpt,
                                vvel_gpt,
                                wvel_gpt,
                                energy_gpt,
                                pressure_gpt,
                                param->gamma );

    double rho_ib_target, uvel_ib_target, vvel_ib_target, wvel_ib_target, pressure_ib_target;
    rho_ib_target = rho;
    pressure_ib_target = pressure;
    uvel_ib_target = -uvel * factor;
    vvel_ib_target = -vvel * factor;
    wvel_ib_target = -wvel * factor;

    double rho_ib, uvel_ib, vvel_ib, wvel_ib, energy_ib, pressure_ib;
    rho_ib      = ramp_fac * rho_ib_target      + (1.0-ramp_fac) * rho_gpt;
    pressure_ib = ramp_fac * pressure_ib_target + (1.0-ramp_fac) * pressure_gpt;
    uvel_ib     = ramp_fac * uvel_ib_target     + (1.0-ramp_fac) * uvel_gpt;
    vvel_ib     = ramp_fac * vvel_ib_target     + (1.0-ramp_fac) * vvel_gpt;
    wvel_ib     = ramp_fac * wvel_ib_target     + (1.0-ramp_fac) * wvel_gpt;
    energy_ib   = inv_gamma_m1*pressure_ib
                  + 0.5*rho_ib*(uvel_ib*uvel_ib+vvel_ib*vvel_ib+wvel_ib*wvel_ib);

    u[_MODEL_NVARS_*node_index+0] = rho_ib;
    u[_MODEL_NVARS_*node_index+1] = rho_ib * uvel_ib;
    u[_MODEL_NVARS_*node_index+2] = rho_ib * vvel_ib;
    u[_MODEL_NVARS_*node_index+3] = rho_ib * wvel_ib;
    u[_MODEL_NVARS_*node_index+4] = energy_ib;
  }

  return 0;
}

/*! Apply no-slip isothermal wall boundary conditions on the immersed boundary
    points (grid points within the immersed body that are within
    stencil-width distance of interior points, i.e., points in the
    interior of the computational domain). */
int NavierStokes3DIBIsothermal( void    *s, /*!< Solver object of type #HyPar */
                                void    *m, /*!< Solver object of type #HyPar */
                                double  *u, /*!< Array with the solution vector */
                                double  t   /*!< Current simulation time */
                              )
{
  HyPar             *solver   = (HyPar*)   s;
  MPIVariables      *mpi      = (MPIVariables*) m;
  ImmersedBoundary  *IB       = (ImmersedBoundary*) solver->ib;
  IBNode            *boundary = IB->boundary;
  NavierStokes3D    *param    = (NavierStokes3D*) solver->physics;
  static double     v[_MODEL_NVARS_];
  int               n, j, k, nb = IB->n_boundary_nodes;

  if (!solver->flag_ib) return(0);

  /* Ideally, this shouldn't be here - But this function is called everywhere
     (through ApplyIBConditions()) *before* MPIExchangeBoundariesnD is called! */
  MPIExchangeBoundariesnD(_MODEL_NDIMS_,_MODEL_NVARS_,solver->dim_local,solver->ghosts,mpi,u);

  double inv_gamma_m1 = 1.0 / (param->gamma - 1.0);

  double ramp_fac = 1.0;
  if (param->t_ib_ramp > 0) {
    double x = t/param->t_ib_ramp;
    if (!strcmp(param->ib_ramp_type,_IB_RAMP_LINEAR_)) {
      ramp_fac = x;
    } else if (!strcmp(param->ib_ramp_type,_IB_RAMP_SMOOTHEDSLAB_)) {
      double a = 0.0;
      double b = 1.0;
      double c = 0.5;
      double r = param->t_ib_ramp/param->t_ib_width;
      ramp_fac = (a*exp(c*r)+b*exp(r*x))/(exp(c*r)+exp(r*x));
    } else if (!strcmp(param->ib_ramp_type,_IB_RAMP_DISABLE_)) {
      ramp_fac = 0.0;
    } else {
      fprintf(stderr,"Error in NavierStokes3DImmersedBoundary():\n");
      fprintf(stderr,"  Ramp type %s not recognized.\n", param->ib_ramp_type);
      return 1;
    }
  }

  for (n=0; n<nb; n++) {

    int     node_index = boundary[n].p;
    double  *alpha = &(boundary[n].interp_coeffs[0]);
    int     *nodes = &(boundary[n].interp_nodes[0]);
    double  factor = boundary[n].surface_distance / boundary[n].interp_node_distance;

    _ArraySetValue_(v,_MODEL_NVARS_,0.0);
    for (j=0; j<_IB_NNODES_; j++) {
      for (k=0; k<_MODEL_NVARS_; k++) {
        v[k] += ( alpha[j] * u[_MODEL_NVARS_*nodes[j]+k] );
      }
    }

    double rho, uvel, vvel, wvel, energy, pressure, temperature;
    _NavierStokes3DGetFlowVar_(v,_NavierStokes3D_stride_,rho,uvel,vvel,wvel,energy,pressure,param->gamma);
    temperature = pressure / rho;

    double rho_gpt, uvel_gpt, vvel_gpt, wvel_gpt, energy_gpt, pressure_gpt, temperature_gpt;
    _NavierStokes3DGetFlowVar_( (u+_MODEL_NVARS_*node_index),
                                _NavierStokes3D_stride_,
                                rho_gpt,
                                uvel_gpt,
                                vvel_gpt,
                                wvel_gpt,
                                energy_gpt,
                                pressure_gpt,
                                param->gamma );
    temperature_gpt = pressure_gpt / rho_gpt;

    double  rho_ib_target,
            uvel_ib_target, vvel_ib_target, wvel_ib_target,
            pressure_ib_target,
            temperature_ib_target;
    temperature_ib_target = (1.0+factor)*param->T_ib_wall - factor * temperature;
    if (    (temperature_ib_target < param->T_ib_wall/param->ib_T_tol)
         || (temperature_ib_target > param->T_ib_wall*param->ib_T_tol) ) {
      temperature_ib_target = param->T_ib_wall;
    }
    pressure_ib_target = pressure;
    rho_ib_target = pressure_ib_target / temperature_ib_target;
    uvel_ib_target = - factor * uvel;
    vvel_ib_target = - factor * vvel;
    wvel_ib_target = - factor * wvel;

    double rho_ib, uvel_ib, vvel_ib, wvel_ib, energy_ib, pressure_ib;
    rho_ib      = ramp_fac * rho_ib_target      + (1.0-ramp_fac) * rho_gpt;
    pressure_ib = ramp_fac * pressure_ib_target + (1.0-ramp_fac) * pressure_gpt;
    uvel_ib     = ramp_fac * uvel_ib_target     + (1.0-ramp_fac) * uvel_gpt;
    vvel_ib     = ramp_fac * vvel_ib_target     + (1.0-ramp_fac) * vvel_gpt;
    wvel_ib     = ramp_fac * wvel_ib_target     + (1.0-ramp_fac) * wvel_gpt;
    energy_ib   = inv_gamma_m1*pressure_ib
                  + 0.5*rho_ib*(uvel_ib*uvel_ib+vvel_ib*vvel_ib+wvel_ib*wvel_ib);

    u[_MODEL_NVARS_*node_index+0] = rho_ib;
    u[_MODEL_NVARS_*node_index+1] = rho_ib * uvel_ib;
    u[_MODEL_NVARS_*node_index+2] = rho_ib * vvel_ib;
    u[_MODEL_NVARS_*node_index+3] = rho_ib * wvel_ib;
    u[_MODEL_NVARS_*node_index+4] = energy_ib;
  }

  return(0);
}

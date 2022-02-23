/*! @file NavierStokes2DParabolicFunction_GPU.cu
    @author Youngdae Kim
    @brief Compute the viscous terms for the 2D Navier Stokes equations
*/
#include <basic_gpu.h>
#include <arrayfunctions_gpu.h>
#include <mathfunctions.h>
#include <physicalmodels/navierstokes2d.h>
#include <mpivars.h>
#include <hypar.h>

#ifdef CUDA_VAR_ORDERDING_AOS

/*! Kernel for gpuNavierStokes2DParabolicFunction() */
__global__
void NavierStokes2DParabolicFunction_Q_kernel(
    int ngrid_points,
    int nvars,
    double gamma,
    const double *u,
    double *Q
)
{
    int p = threadIdx.x + (blockDim.x * blockIdx.x);
    if (p < ngrid_points) {
        double energy, pressure;
        p *= nvars;
        _NavierStokes2DGetFlowVar_( (u+p),Q[p+0],Q[p+1],Q[p+2],energy,
                                    pressure,gamma);
        Q[p+3] = gamma*pressure/Q[p+0];
    }

    return;
}

/*!
    Compute the viscous terms in the 2D Navier Stokes equations: this function computes
    the following:
    \f{equation}{
      \frac {\partial} {\partial x} \left[\begin{array}{c} 0 \\ \tau_{xx} \\ \tau_{yx} \\ u \tau_{xx} + v \tau_{yx} - q_x \end{array}\right]
      + \frac {\partial} {\partial y} \left[\begin{array}{c} 0 \\ \tau_{xy} \\ \tau_{yy} \\ u \tau_{xy} + v \tau_{yy} - q_y \end{array}\right]
    \f}
    where
    \f{align}{
      \tau_{xx} &= \frac{2}{3}\left(\frac{\mu}{Re}\right)\left(2\frac{\partial u}{\partial x} - \frac{\partial v}{\partial y}\right),\\
      \tau_{xy} &= \left(\frac{\mu}{Re}\right)\left(\frac{\partial u}{\partial y} + \frac{\partial v}{\partial x}\right),\\
      \tau_{yx} &= \left(\frac{\mu}{Re}\right)\left(\frac{\partial u}{\partial y} + \frac{\partial v}{\partial x}\right),\\
      \tau_{yy} &= \frac{2}{3}\left(\frac{\mu}{Re}\right)\left(-\frac{\partial u}{\partial x} +2\frac{\partial v}{\partial y}\right),\\
      q_x &= -\frac{mu}{\left(\gamma-1\right)Re Pr}\frac{\partial T}{\partial x}, \\
      q_y &= -\frac{mu}{\left(\gamma-1\right)Re Pr}\frac{\partial T}{\partial y}
    \f}
    and the temperature is \f$T = \gamma p/\rho\f$. \f$Re\f$ and \f$Pr\f$ are the Reynolds and Prandtl numbers, respectively. Note that this function
    computes the entire parabolic term, and thus bypasses HyPar's parabolic function calculation interfaces. NavierStokes2DInitialize() assigns this
    function to #HyPar::ParabolicFunction.
    \n\n
    Reference:
    + Tannehill, Anderson and Pletcher, Computational Fluid Mechanics and Heat Transfer,
      Chapter 5, Section 5.1.7 (However, the non-dimensional velocity and the Reynolds
      number is based on speed of sound, instead of the freestream velocity).
*/
extern "C" int gpuNavierStokes2DParabolicFunction(
    double  *par, /*!< Array to hold the computed viscous terms */
    double  *u,   /*!< Solution vector array */
    void    *s,   /*!< Solver object of type #HyPar */
    void    *m,   /*!< MPI object of type #MPIVariables */
    double  t     /*!< Current simulation time */
)
{
  HyPar           *solver   = (HyPar*) s;
  NavierStokes2D  *physics  = (NavierStokes2D*) solver->physics;
  _DECLARE_IERR_;

  int ghosts = solver->ghosts;
  int imax   = solver->dim_local[0];
  int jmax   = solver->dim_local[1];
  int nvars  = solver->nvars;
  int size   = (imax+2*ghosts)*(jmax+2*ghosts)*nvars;

  gpuArraySetValue(par, size, 0.0);
  if (physics->Re <= 0) return(0); /* inviscid flow */

  printf("gpuNavierStokes2DParabolicFunction() hasn't been fully implemented yet.\n");
  return 1;
}

#endif

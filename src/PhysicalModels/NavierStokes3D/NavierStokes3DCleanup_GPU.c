/*! @file NavierStokes3DCleanup_GPU.c
    @author Youngdae Kim
    @brief Clean up the 3D Navier Stokes module in the GPU memory
*/

#include <stdlib.h>
#include <physicalmodels/navierstokes3d.h>
#include <arrayfunctions_gpu.h>

/*! Function to clean up all allocations in the 3D Navier
    Stokes module.
*/
int gpuNavierStokes3DCleanup(void *s /*!< Object of type #NavierStokes3D*/)
{
  NavierStokes3D  *param  = (NavierStokes3D*) s;

  gpuFree(param->gpu_Q);
  gpuFree(param->gpu_QDerivX);
  gpuFree(param->gpu_QDerivY);
  gpuFree(param->gpu_QDerivZ);
  gpuFree(param->gpu_FViscous);
  gpuFree(param->gpu_FDeriv);
  gpuFree(param->gpu_grav_field_f);
  gpuFree(param->gpu_grav_field_g);
  gpuFree(param->gpu_fast_jac);
  gpuFree(param->gpu_solution);

  return(0);
}

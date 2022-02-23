/*! @file NavierStokes2DInitialize_GPU.c
    @author Youngdae Kim
    @brief Initialization of the physics-related variables and function pointers for the 2D Navier-Stokes system
*/
#include <stdio.h>
#include <arrayfunctions_gpu.h>
#include <physicalmodels/navierstokes2d.h>
#include <hypar.h>

/*! Initialize GPU-related arrays. */
int gpuNavierStokes2DInitialize(
    void *s,  /*!< Solver object of type #HyPar */
    void *m   /*!< MPI object of type #MPIVariables */
)
{
    HyPar          *solver     = (HyPar*)          s;
    NavierStokes2D *physics    = (NavierStokes2D*) solver->physics;

    int *dim = solver->dim_local;
    int ghosts = solver->ghosts;
    int d, size = 1; for (d = 0; d <_MODEL_NDIMS_; d++) size *= (dim[d] + 2*ghosts);

    gpuMalloc((void**)&physics->gpu_grav_field_f, size*sizeof(double));
    gpuMalloc((void**)&physics->gpu_grav_field_g, size*sizeof(double));
    gpuMalloc((void**)&physics->gpu_fast_jac, 2*size*_MODEL_NVARS_*_MODEL_NVARS_*sizeof(double));
    gpuMalloc((void**)&physics->gpu_solution, size*_MODEL_NVARS_*sizeof(double));
    gpuMemcpy(physics->gpu_grav_field_f, physics->grav_field_f, size*sizeof(double), gpuMemcpyHostToDevice);
    gpuMemcpy(physics->gpu_grav_field_g, physics->grav_field_g, size*sizeof(double), gpuMemcpyHostToDevice);
    gpuMemset(physics->gpu_fast_jac, 0, 2*size*_MODEL_NVARS_*_MODEL_NVARS_*sizeof(double));
    gpuMemset(physics->gpu_solution, 0, size*_MODEL_NVARS_*sizeof(double));

    return (0);
}

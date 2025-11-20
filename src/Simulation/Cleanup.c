/*! @file Cleanup.c
    @author Debojyoti Ghosh
    @brief Clean up and free memory after simulation is complete.
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <basic.h>
#include <bandedmatrix.h>
#include <tridiagLU.h>
#include <boundaryconditions.h>
#include <immersedboundaries.h>
#include <timeintegration.h>
#include <interpolation.h>
#include <mpivars.h>
#include <simulation_object.h>

/* include header files for each physical model */
#include <physicalmodels/linearadr.h>
#include <physicalmodels/fpdoublewell.h>
#include <physicalmodels/fppowersystem.h>
#include <physicalmodels/fppowersystem1bus.h>
#include <physicalmodels/fppowersystem3bus.h>
#include <physicalmodels/euler1d.h>
#include <physicalmodels/euler2d.h>
#include <physicalmodels/navierstokes2d.h>
#include <physicalmodels/navierstokes3d.h>
#include <physicalmodels/numa2d.h>
#include <physicalmodels/numa3d.h>
#include <physicalmodels/shallowwater1d.h>
#include <physicalmodels/shallowwater2d.h>
#include <physicalmodels/vlasov.h>

#if defined(HAVE_CUDA)
#include <arrayfunctions_gpu.h>
#endif

/*! Cleans up and frees the memory after the completion of the simulation. */
int Cleanup(  void  *s,   /*!< Array of simulation objects of type #SimulationObject */
              int   nsims /*!< number of simulation objects */
           )
{
  SimulationObject* sim = (SimulationObject*) s;
  int ns;
  _DECLARE_IERR_;

  if (nsims == 0) return 0;

  if (!sim[0].mpi.m_rank) {
    printf("Deallocating arrays.\n");
  }

  for (ns = 0; ns < nsims; ns++) {

    if (sim[ns].is_barebones == 1) {
      fprintf(stderr, "Error in Cleanup(): object is barebones type.\n");
      return 1;
    }

    HyPar* solver = &(sim[ns].solver);
    MPIVariables* mpi = &(sim[ns].mpi);
    DomainBoundary* boundary = (DomainBoundary*) solver->m_boundary;
    int i;

    /* Clean up boundary zones */
    for (i = 0; i < solver->m_n_boundary_zones; i++) {
#if defined(HAVE_CUDA)
      BCCleanup(&boundary[i], solver->m_use_gpu);
#else
      BCCleanup(&boundary[i], 0);
#endif
    }
    free(solver->m_boundary);

    /* Clean up immersed boundaries */
    if (solver->m_flag_ib) {
      IERR IBCleanup(solver->m_ib);
      free(solver->m_ib);
    }

    /* Clean up any allocations in physical model */
    if (!strcmp(solver->m_model,_LINEAR_ADVECTION_DIFFUSION_REACTION_)) {
      IERR LinearADRCleanup(solver->m_physics); CHECKERR(ierr);
    } else if (!strcmp(solver->m_model,_FP_DOUBLE_WELL_)) {
      IERR FPDoubleWellCleanup(solver->m_physics); CHECKERR(ierr);
    } else if (!strcmp(solver->m_model,_FP_POWER_SYSTEM_)) {
      IERR FPPowerSystemCleanup(solver->m_physics); CHECKERR(ierr);
    } else if (!strcmp(solver->m_model,_FP_POWER_SYSTEM_1BUS_)) {
      IERR FPPowerSystem1BusCleanup(solver->m_physics); CHECKERR(ierr);
    } else if (!strcmp(solver->m_model,_FP_POWER_SYSTEM_3BUS_)) {
      IERR FPPowerSystem3BusCleanup(solver->m_physics); CHECKERR(ierr);
    } else if (!strcmp(solver->m_model,_EULER_1D_)) {
      IERR Euler1DCleanup(solver->m_physics); CHECKERR(ierr);
    } else if (!strcmp(solver->m_model,_EULER_2D_)) {
      IERR Euler2DCleanup(solver->m_physics); CHECKERR(ierr);
    } else if (!strcmp(solver->m_model,_NAVIER_STOKES_2D_)) {
      IERR NavierStokes2DCleanup(solver->m_physics); CHECKERR(ierr);
    } else if (!strcmp(solver->m_model,_NAVIER_STOKES_3D_)) {
      IERR NavierStokes3DCleanup(solver->m_physics); CHECKERR(ierr);
#if defined(HAVE_CUDA)
      if (solver->m_use_gpu) {
        IERR gpuNavierStokes3DCleanup(solver->m_physics); CHECKERR(ierr);
      }
#endif
    } else if (!strcmp(solver->m_model,_NUMA2D_)) {
      IERR Numa2DCleanup(solver->m_physics); CHECKERR(ierr);
    } else if (!strcmp(solver->m_model,_NUMA3D_)) {
      IERR Numa3DCleanup(solver->m_physics); CHECKERR(ierr);
    } else if (!strcmp(solver->m_model,_SHALLOW_WATER_1D_)) {
      IERR ShallowWater1DCleanup(solver->m_physics); CHECKERR(ierr);
    } else if (!strcmp(solver->m_model,_SHALLOW_WATER_2D_)) {
      IERR ShallowWater2DCleanup(solver->m_physics); CHECKERR(ierr);
    } else if (!strcmp(solver->m_model,_VLASOV_)) {
      IERR VlasovCleanup(solver->m_physics); CHECKERR(ierr);
    }
    free(solver->m_physics);

    /* Clean up any allocations from time-integration */
#ifdef with_petsc
    if (!solver->m_use_petsc_ts) {
      if (!strcmp(solver->m_time_scheme,_RK_)) {
        IERR TimeExplicitRKCleanup(solver->m_msti); CHECKERR(ierr);
        free(solver->m_msti);
      } else if (!strcmp(solver->m_time_scheme,_GLM_GEE_)) {
        IERR TimeGLMGEECleanup(solver->m_msti); CHECKERR(ierr);
        free(solver->m_msti);
      }
    }
#else
    if (!strcmp(solver->m_time_scheme,_RK_)) {
      IERR TimeExplicitRKCleanup(solver->m_msti); CHECKERR(ierr);
      free(solver->m_msti);
    } else if (!strcmp(solver->m_time_scheme,_GLM_GEE_)) {
      IERR TimeGLMGEECleanup(solver->m_msti); CHECKERR(ierr);
      free(solver->m_msti);
    }
#endif

    /* Clean up any spatial reconstruction related allocations */
    if (   (!strcmp(solver->m_spatial_scheme_hyp,_FIFTH_ORDER_WENO_  ))
        || (!strcmp(solver->m_spatial_scheme_hyp,_FIFTH_ORDER_CRWENO_))
        || (!strcmp(solver->m_spatial_scheme_hyp,_FIFTH_ORDER_HCWENO_)) ) {
#if defined(HAVE_CUDA)
      IERR WENOCleanup(solver->m_interp, solver->m_use_gpu); CHECKERR(ierr);
#else
      IERR WENOCleanup(solver->m_interp, 0); CHECKERR(ierr);
#endif
    }
    if (solver->m_interp)   free(solver->m_interp);
    if (   (!strcmp(solver->m_spatial_scheme_hyp,_FIFTH_ORDER_COMPACT_UPWIND_ ))
        || (!strcmp(solver->m_spatial_scheme_hyp,_FIFTH_ORDER_CRWENO_         ))
        || (!strcmp(solver->m_spatial_scheme_hyp,_FIFTH_ORDER_HCWENO_         )) ) {
      IERR CompactSchemeCleanup(solver->m_compact); CHECKERR(ierr);
    }
    if (solver->m_compact)  free(solver->m_compact);
    if (solver->m_lusolver) free(solver->m_lusolver);

    /* Free the communicators created */
    IERR MPIFreeCommunicators(solver->m_ndims,mpi); CHECKERR(ierr);

    /* These variables are allocated in Initialize.c */
    free(solver->m_dim_global);
    free(solver->m_dim_global_ex);
    free(solver->m_dim_local);
    free(solver->m_index);
    free(solver->m_u);
#ifdef with_petsc
    if (solver->m_u0)     free(solver->m_u0);
    if (solver->m_uref)   free(solver->m_uref);
    if (solver->m_rhsref) free(solver->m_rhsref);
    if (solver->m_rhs)    free(solver->m_rhs);
#endif
#ifdef with_librom
    free(solver->m_u_rom_predicted);
#endif
    free(solver->m_iblank);
    free(solver->m_x);
    free(solver->m_dxinv);
    free(solver->m_is_periodic);
    free(mpi->m_iproc);
    free(mpi->m_ip);
    free(mpi->m_is);
    free(mpi->m_ie);
    free(mpi->m_bcperiodic);
    free(mpi->m_sendbuf);
    free(mpi->m_recvbuf);
    free(solver->m_volume_integral);
    free(solver->m_volume_integral_initial);
    free(solver->m_total_boundary_integral);
    free(solver->m_conservation_error);
    free(solver->m_stride_with_ghosts);
    free(solver->m_stride_without_ghosts);

#if defined(HAVE_CUDA)
    if (solver->m_use_gpu) {
      gpuFree(solver->m_hyp);
      gpuFree(solver->m_par);
      gpuFree(solver->m_source);
      gpuFree(solver->m_u_c);
      gpuFree(solver->m_flux_c);
      gpuFree(solver->m_deriv1);
      gpuFree(solver->m_deriv2);
      gpuFree(solver->m_flux_i);
      gpuFree(solver->m_u_l);
      gpuFree(solver->m_u_r);
      gpuFree(solver->m_f_l);
      gpuFree(solver->m_f_r);
      gpuFree(solver->m_stage_boundary_buffer);
      gpuFree(solver->m_stage_boundary_integral);
      gpuFree(solver->m_step_boundary_integral);

      gpuFree(solver->m_gpu_dim_local);
      gpuFree(solver->m_gpu_iblank);
      gpuFree(solver->m_gpu_x);
      gpuFree(solver->m_gpu_dxinv);
      gpuFree(solver->m_gpu_u);
    } else {
#endif
      free(solver->m_hyp);
      free(solver->m_par);
      free(solver->m_source);
      free(solver->m_u_c);
      free(solver->m_flux_c);
      free(solver->m_deriv1);
      free(solver->m_deriv2);
      free(solver->m_flux_i);
      free(solver->m_u_l);
      free(solver->m_u_r);
      free(solver->m_f_l);
      free(solver->m_f_r);
      free(solver->m_stage_boundary_integral);
      free(solver->m_step_boundary_integral);
#if defined(HAVE_CUDA)
    }
#endif

    if (solver->m_filename_index) free(solver->m_filename_index);

  }

  return(0);
}

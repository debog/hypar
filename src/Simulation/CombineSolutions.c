/*! @file CombineSolutions.c
    @brief Functions to combine solutions on multiple grids on a target grid
    @author Debojyoti Ghosh
*/

#include <arrayfunctions.h>
#include <mathfunctions.h>
#include <mpivars.h>
#include <simulation_object.h>

/*! This function combines solutions on multiple grids on a target grid
    using given coefficients.

    The source grids may have varying processor layouts, so they are
    all gathered on to rank 0, interpolated on to the target grid, and
    combined. The result is partitioned to the destination processor
    layout.
*/
void CombineSolutions(  SimulationObject*     a_sims_src, /*!< Array of simulation objects of type #SimulationObject */
                        double* const* const  a_u_src,    /*!< Array of source solutions */
                        const int             a_nsims,    /*!< Number of simulation objects */
                        SimulationObject*     a_sim_dst,  /*!< simulation object corresponding to destination */
                        double* const         a_u_dst,    /*!< Destination solution array */
                        const double* const   a_coeffs    /*!< Combination coefficients for the source solutions */
                     )
{
  HyPar* solver_dst = &(a_sim_dst->solver);
  MPIVariables* mpi_dst = &(a_sim_dst->mpi);

  int ndims = solver_dst->ndims;
  int nvars = solver_dst->nvars;
  int ghosts = solver_dst->ghosts;

  /* array size */
  int* dim_global_dst = solver_dst->dim_global;
  int dim_dst_wgpt[ndims];
  _ArrayCopy1D_(dim_global_dst, dim_dst_wgpt, ndims);
  long size_dst_wgpt = nvars;
  for (int n=0; n<ndims; n++) {
    dim_dst_wgpt[n] += 2*ghosts;
    size_dst_wgpt *= dim_dst_wgpt[n];
  }

  /* allocate global data array for combined solution */
  double *ug_dst_wgpt = NULL;
  if (!mpi_dst->rank) {
    ug_dst_wgpt = (double*) calloc (size_dst_wgpt, sizeof(double));
    _ArraySetValue_(ug_dst_wgpt, size_dst_wgpt, 0.0);
  }

  /* for each sparse grids, interpolate onto the dst grid dimension
   * and add to the solution */
  for (int n = 0; n < a_nsims; n++) {

    int* dim_global_src = a_sims_src[n].solver.dim_global;
    int dim_global_src_wgpt[ndims];
    _ArrayCopy1D_(dim_global_src, dim_global_src_wgpt, ndims);
    long size_src_wgpt = nvars;
    for (int n=0; n<ndims; n++) {
      dim_global_src_wgpt[n] += 2*ghosts;
      size_src_wgpt *= dim_global_src_wgpt[n];
    }

    double* ug_src_wgpt = NULL;
    if (!a_sims_src[n].mpi.rank) {
      ug_src_wgpt = (double*) calloc (size_src_wgpt, sizeof(double));
      _ArraySetValue_(ug_src_wgpt, size_src_wgpt, 0.0);
    }
    MPIGatherArraynDwGhosts( ndims,
                             (void*) &(a_sims_src[n].mpi),
                             ug_src_wgpt,
                             a_u_src[n],
                             dim_global_src,
                             a_sims_src[n].solver.dim_local,
                             ghosts,
                             nvars );

    if (!mpi_dst->rank) {
      double *u_src_interpolated = NULL;
      InterpolateGlobalnDVar( dim_global_dst,
                              &u_src_interpolated,
                              a_sims_src[n].solver.dim_global,
                              ug_src_wgpt,
                              nvars,
                              ghosts,
                              ndims,
                              a_sims_src[n].solver.isPeriodic );
      _ArrayAXPY_(u_src_interpolated, a_coeffs[n], ug_dst_wgpt, size_dst_wgpt);
      free(u_src_interpolated);
    }

  }

  /* now partition the global combined solution to its processor */
  MPIPartitionArraynDwGhosts( ndims,
                              mpi_dst,
                              (mpi_dst->rank ? NULL : ug_dst_wgpt),
                              a_u_dst,
                              solver_dst->dim_global,
                              solver_dst->dim_local,
                              ghosts,
                              nvars );

  /* free memory */
  if (!mpi_dst->rank) free(ug_dst_wgpt);

  /* done */
  return;
}


/*! @file SparseGridsSetSolverParameters.cpp
    @author Debojyoti Ghosh, John Loffeld, Lee Ricketson
    @brief Set the solver parameters of a simulation object
*/

#include <sparse_grids_simulation.h>

/*! Set the solver parameters (stuff that is usually read in from solver.inp) for
 *  a simulation object, where the global grid sizes and processor distribution are
 *  specified, and all other parameters are same as a source simulation object */
int SparseGridsSimulation::SetSolverParameters( SimulationObject&       a_dst_sim,    /*!< Simulation object */
                                                const GridDimensions&   a_dim_global, /*!< Specified global grid sizes */
                                                const ProcDistribution& a_iproc,      /*!< Specified processor distibution */
                                                const SimulationObject& a_src_sim,    /*!< Source simulation object */
                                                const int               a_idx,        /*!< Index or serial number */
                                                const int               a_nsims       /*!< Total number of simulations */
                                              )
{
  a_dst_sim.solver.m_my_idx = a_idx;
  a_dst_sim.solver.m_nsims = a_nsims;

  a_dst_sim.mpi.m_rank = m_rank;
  a_dst_sim.mpi.m_nproc = m_nproc;

  a_dst_sim.solver.m_ndims = a_src_sim.solver.m_ndims;
  a_dst_sim.solver.m_nvars = a_src_sim.solver.m_nvars;
  a_dst_sim.solver.m_ghosts = a_src_sim.solver.m_ghosts;

  a_dst_sim.solver.m_dim_global = (int*) calloc (a_dst_sim.solver.m_ndims,sizeof(int));
  a_dst_sim.mpi.m_iproc         = (int*) calloc (a_dst_sim.solver.m_ndims,sizeof(int));
  for (int d = 0; d < a_dst_sim.solver.m_ndims; d++) {
    a_dst_sim.solver.m_dim_global[d] = a_dim_global[d];
    a_dst_sim.mpi.m_iproc[d] = a_iproc[d];
  }

  a_dst_sim.solver.m_n_iter       = a_src_sim.solver.m_n_iter;
  a_dst_sim.solver.m_restart_iter = a_src_sim.solver.m_restart_iter;

  strcpy(a_dst_sim.solver.m_time_scheme, a_src_sim.solver.m_time_scheme);
  strcpy(a_dst_sim.solver.m_time_scheme_type, a_src_sim.solver.m_time_scheme_type);
  strcpy(a_dst_sim.solver.m_spatial_scheme_hyp, a_src_sim.solver.m_spatial_scheme_hyp);
  strcpy(a_dst_sim.solver.m_split_hyperbolic_flux, a_src_sim.solver.m_split_hyperbolic_flux);
  strcpy(a_dst_sim.solver.m_interp_type, a_src_sim.solver.m_interp_type);
  strcpy(a_dst_sim.solver.m_spatial_type_par, a_src_sim.solver.m_spatial_type_par);
  strcpy(a_dst_sim.solver.m_spatial_scheme_par, a_src_sim.solver.m_spatial_scheme_par);

  a_dst_sim.solver.m_dt = a_src_sim.solver.m_dt;

  strcpy(a_dst_sim.solver.m_conservation_check, a_src_sim.solver.m_conservation_check);

  a_dst_sim.solver.m_screen_op_iter = a_src_sim.solver.m_screen_op_iter;
  a_dst_sim.solver.m_file_op_iter = a_src_sim.solver.m_file_op_iter;

  strcpy(a_dst_sim.solver.m_op_file_format, a_src_sim.solver.m_op_file_format);
  strcpy(a_dst_sim.solver.m_ip_file_type, a_src_sim.solver.m_ip_file_type);

  strcpy(a_dst_sim.solver.m_input_mode, a_src_sim.solver.m_input_mode);
  strcpy(a_dst_sim.solver.m_output_mode, a_src_sim.solver.m_output_mode);
  a_dst_sim.mpi.m_N_IORanks = a_src_sim.mpi.m_N_IORanks;

  strcpy(a_dst_sim.solver.m_op_overwrite, a_src_sim.solver.m_op_overwrite);
  strcpy(a_dst_sim.solver.m_plot_solution, a_src_sim.solver.m_plot_solution);
  strcpy(a_dst_sim.solver.m_model, a_src_sim.solver.m_model);
  strcpy(a_dst_sim.solver.m_ib_filename, a_src_sim.solver.m_ib_filename);

  a_dst_sim.solver.m_flag_ib = a_src_sim.solver.m_flag_ib;

#ifdef with_petsc
  a_dst_sim.solver.m_use_petsc_ts = a_src_sim.solver.m_use_petsc_ts;
#endif

  return(0);
}

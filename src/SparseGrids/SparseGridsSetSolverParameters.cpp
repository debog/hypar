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
  a_dst_sim.solver.my_idx = a_idx;
  a_dst_sim.solver.nsims = a_nsims;

  a_dst_sim.mpi.rank = m_rank;
  a_dst_sim.mpi.nproc = m_nproc;

  a_dst_sim.solver.ndims = a_src_sim.solver.ndims;
  a_dst_sim.solver.nvars = a_src_sim.solver.nvars;
  a_dst_sim.solver.ghosts = a_src_sim.solver.ghosts;

  a_dst_sim.solver.dim_global = (int*) calloc (a_dst_sim.solver.ndims,sizeof(int));
  a_dst_sim.mpi.iproc         = (int*) calloc (a_dst_sim.solver.ndims,sizeof(int));
  for (int d = 0; d < a_dst_sim.solver.ndims; d++) {
    a_dst_sim.solver.dim_global[d] = a_dim_global[d];
    a_dst_sim.mpi.iproc[d] = a_iproc[d];
  }

  a_dst_sim.solver.n_iter       = a_src_sim.solver.n_iter;
  a_dst_sim.solver.restart_iter = a_src_sim.solver.restart_iter;

  strcpy(a_dst_sim.solver.time_scheme, a_src_sim.solver.time_scheme);
  strcpy(a_dst_sim.solver.time_scheme_type, a_src_sim.solver.time_scheme_type);
  strcpy(a_dst_sim.solver.spatial_scheme_hyp, a_src_sim.solver.spatial_scheme_hyp);
  strcpy(a_dst_sim.solver.SplitHyperbolicFlux, a_src_sim.solver.SplitHyperbolicFlux);
  strcpy(a_dst_sim.solver.interp_type, a_src_sim.solver.interp_type);
  strcpy(a_dst_sim.solver.spatial_type_par, a_src_sim.solver.spatial_type_par);
  strcpy(a_dst_sim.solver.spatial_scheme_par, a_src_sim.solver.spatial_scheme_par);

  a_dst_sim.solver.dt = a_src_sim.solver.dt;

  strcpy(a_dst_sim.solver.ConservationCheck, a_src_sim.solver.ConservationCheck);

  a_dst_sim.solver.screen_op_iter = a_src_sim.solver.screen_op_iter;
  a_dst_sim.solver.file_op_iter = a_src_sim.solver.file_op_iter;

  strcpy(a_dst_sim.solver.op_file_format, a_src_sim.solver.op_file_format);
  strcpy(a_dst_sim.solver.ip_file_type, a_src_sim.solver.ip_file_type);

  strcpy(a_dst_sim.solver.input_mode, a_src_sim.solver.input_mode);
  strcpy(a_dst_sim.solver.output_mode, a_src_sim.solver.output_mode);
  a_dst_sim.mpi.N_IORanks = a_src_sim.mpi.N_IORanks;

  strcpy(a_dst_sim.solver.op_overwrite, a_src_sim.solver.op_overwrite);
  strcpy(a_dst_sim.solver.plot_solution, a_src_sim.solver.plot_solution);
  strcpy(a_dst_sim.solver.model, a_src_sim.solver.model);
  strcpy(a_dst_sim.solver.ib_filename, a_src_sim.solver.ib_filename);

  a_dst_sim.solver.flag_ib = a_src_sim.solver.flag_ib;

#ifdef with_petsc
  a_dst_sim.solver.use_petscTS = a_src_sim.solver.use_petscTS;
#endif

  return(0);
}

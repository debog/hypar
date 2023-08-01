#ifdef with_librom

/*! @file LSROMObject.cpp
 *  @brief Member functions of the class #LSROMObject
 *  @author Ping-Hsuan Tsai
*/

#include <basic.h>
#include <rom_object_ls.h>
#include <simulation_object.h>
#include <physicalmodels/vlasov.h>
//#include <timeintegration.h>
#include <arrayfunctions.h>
#include <io_cpp.h>
#include <mpivars.h>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>


/*! Constructor 
    This function will also look into the file
    \b librom.inp for LS-specific options. 
    Rank 0 reads in the inputs and broadcasts
    them to all the processors.\n\n
    The format of \b librom.inp is as follows:\n

        begin
            <keyword>   <value>
            ...
            <keyword>   <value>
        end

    where the list of keywords and their type (that this function
    will look for) are:\n
    Keyword name           | Type         | Variable                                      | Default value
    ---------------------- | ------------ | --------------------------------------------- | -------------------
    ls_num_win_samples    | int          | #LSROMObject::m_num_window_samples           | INT_MAX
    ls_dirname            | string       | #LSROMObject::m_dirname                      | "LS"
    ls_write_snapshot_mat | bool         | #LSROMObject::m_write_snapshot_mat           | false

    Note: other keywords in this file may be read by other functions.
   
*/
extern "C" int  TimeRHSFunctionExplicit(double*,double*,void*,void*,double);
extern "C" int  CalculateROMDiff(void*,void*);
extern "C" void ResetFilenameIndex(char*, int); /*!< Reset filename index */
extern "C" void IncrementFilenameIndex(char*,int);
extern "C" int  VlasovWriteSpatialField(void*, void*, double*, char*);
extern "C" int FirstDerivativeSecondOrderCentralNoGhosts(double*,double*,int,
                                                         int,int,int*,int,
                                                         int,void*);
extern "C" int SecondDerivativeSecondOrderCentralNoGhosts(double*,double*,int,
                                                          int,int*,int,int,void*);

extern "C" int VlasovAdvection_x(double*,double*,int,void*,double);
extern "C" int HyperbolicFunction_1dir(double*,double*,void*,void*,double,int,
                                       int(*)(double*,double*,int,void*,double),
                                       int(*)(double*,double*,double*,double*,
                                       double*,double*,int,void*,double),int);

LSROMObject::LSROMObject(   const int     a_vec_size,       /*!< vector size */
                            const double  a_dt,             /*!< time step size */
                            const int     a_sampling_freq,  /*!< sampling frequency */
                            const int     a_rdim,           /*!< latent space dimension */
                            const int     a_parametric_id,  /*!< Vector component index */
                            const int     a_rank,           /*!< MPI rank of this process */
                            const int     a_nproc,          /*!< Number of MPI processes */
                            const int     a_sim_idx,        /*!< Simulation index (for ensemble simulations */
                            const int     a_var_idx         /*!< Vector component index */ )
{
  m_rank = a_rank;
  m_nproc = a_nproc;

  /* my implementation */
  m_options.clear();
  m_generator.clear();
  m_projected_init.clear();
  m_snapshots.clear();
  m_romhyperb.clear();
  m_snap.clear();
  m_basis.clear();

  /* precomputation idea */
  m_options_phi.clear();
  m_generator_phi.clear();
  m_projected_init_phi.clear();
  m_snapshots_phi.clear();
  m_basis_phi.clear();
  m_snapshots_e.clear();
  m_basis_e.clear();

  m_ls_is_trained.clear();
  m_intervals.clear();
  m_vec_size = a_vec_size;
  m_dt = a_dt;
  m_sampling_freq = a_sampling_freq;
  m_t_final = -1;
  m_rdim = a_rdim;
  m_rdim_phi = m_rdim;
  m_rdims.clear();
  m_rdims_phi.clear();

  m_sim_idx = a_sim_idx;
  m_var_idx = a_var_idx;
  m_parametric_id = a_parametric_id;

  m_num_window_samples = INT_MAX;
  m_energy_criteria = -1;
	m_nsets = 1;

  char dirname_c_str[_MAX_STRING_SIZE_] = "LS";
  char write_snapshot_mat[_MAX_STRING_SIZE_] = "false";
  char direct_comp_hyperbolic[_MAX_STRING_SIZE_] = "false";
  char solve_phi[_MAX_STRING_SIZE_] = "false";
  char c_err_snap[_MAX_STRING_SIZE_] = "false";

  if (!m_rank) {

    FILE *in;
    in = fopen(_LIBROM_INP_FNAME_,"r");

    if (in) {

      int ferr;
      char word[_MAX_STRING_SIZE_];
      ferr = fscanf(in,"%s", word); if (ferr != 1) return;

      if (std::string(word) == "begin") {
        while (std::string(word) != "end") {
  	      ferr = fscanf(in,"%s",word); if (ferr != 1) return;
          if (std::string(word) == "ls_num_win_samples") {
            ferr = fscanf(in,"%d", &m_num_window_samples); if (ferr != 1) return;
          } else if (std::string(word) == "ls_dirname") {
            ferr = fscanf(in,"%s", dirname_c_str); if (ferr != 1) return;
          } else if (std::string(word) == "ls_write_snapshot_mat") {
            ferr = fscanf(in,"%s", write_snapshot_mat); if (ferr != 1) return;
          } else if (std::string(word) == "ls_direct_comp_hyperbolic") {
            ferr = fscanf(in,"%s", direct_comp_hyperbolic); if (ferr != 1) return;
          } else if (std::string(word) == "ls_rdim_phi") {
            ferr = fscanf(in,"%d", &m_rdim_phi); if (ferr != 1) return;
          } else if (std::string(word) == "ls_solve_phi") {
            ferr = fscanf(in,"%s", solve_phi); if (ferr != 1) return;
          } else if (std::string(word) == "ls_energy_criteria") {
            ferr = fscanf(in,"%le", &m_energy_criteria); if (ferr != 1) return;
          } else if (std::string(word) == "ls_c_err_snap") {
            ferr = fscanf(in,"%s", c_err_snap); if (ferr != 1) return;
          } else if (std::string(word) == "ls_nsets") {
            ferr = fscanf(in,"%d", &m_nsets); if (ferr != 1) return;
          }
          if (ferr != 1) return;
        }
      } else {
        fprintf( stderr, "Error: Illegal format in file \"%s\". Word read is: %s\n",
                 _LIBROM_INP_FNAME_, word);
        return;
      }

      fclose(in);

    }

    /* print useful stuff to screen */
    printf("LSROMObject details:\n");
    printf("  number of samples per window:   %d\n", m_num_window_samples);
    printf("  SVD energy criteria:   %e\n", m_energy_criteria);
    printf("  potential latent space dimension:  %d\n", m_rdim_phi);
    printf("  directory name for LS objects: %s\n", dirname_c_str);
    printf("  write snapshot matrix to file:  %s\n", write_snapshot_mat);
    printf("  directly compute hyperbolic term:  %s\n", direct_comp_hyperbolic);
    printf("  solve potential:  %s\n", solve_phi);
    printf("  compute error for each snapshot:  %s\n", c_err_snap);
    printf("  number of parametric snapshot sets:   %d\n", m_nsets);
    if (m_sim_idx >= 0) {
      printf("  simulation domain:  %d\n", m_sim_idx);
    }
    if (m_var_idx >= 0) {
      printf("  Vector component:  %d\n", m_var_idx);
    }
  }

#ifndef serial
  MPI_Bcast(&m_rdim_phi,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&m_num_window_samples,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(dirname_c_str,_MAX_STRING_SIZE_,MPI_CHAR,0,MPI_COMM_WORLD);
  MPI_Bcast(write_snapshot_mat,_MAX_STRING_SIZE_,MPI_CHAR,0,MPI_COMM_WORLD);
  MPI_Bcast(direct_comp_hyperbolic,_MAX_STRING_SIZE_,MPI_CHAR,0,MPI_COMM_WORLD);
  MPI_Bcast(solve_phi,_MAX_STRING_SIZE_,MPI_CHAR,0,MPI_COMM_WORLD);
  MPI_Bcast(&m_energy_criteria,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Bcast(c_err_snap,_MAX_STRING_SIZE_,MPI_CHAR,0,MPI_COMM_WORLD);
  MPI_Bcast(&m_nsets,1,MPI_INT,0,MPI_COMM_WORLD);
#endif

  m_dirname = std::string( dirname_c_str );
  m_write_snapshot_mat = (std::string(write_snapshot_mat) == "true");
  m_direct_comp_hyperbolic = (std::string(direct_comp_hyperbolic) == "true");
  m_solve_phi = (std::string(solve_phi) == "true");
  m_c_err_snap = (std::string(c_err_snap) == "true");

  if (m_num_window_samples <= m_rdim) {
    printf("ERROR:LSROMObject::LSROMObject() - m_num_window_samples <= m_rdim!!");
  }

  if (m_rdim_phi > m_rdim) {
    printf("ERROR:LSROMObject::LSROMObject() - m_rdim_phi > m_rdim!!");
  }

  m_tic = 0;
  m_curr_win = 0;
}

void LSROMObject::projectInitialSolution(  CAROM::Vector& a_U, /*!< solution vector */
                                           void* a_s)
{
  SimulationObject* sim = (SimulationObject*) a_s;
  HyPar  *solver = (HyPar*) &(sim[0].solver);
  Vlasov *param  = (Vlasov*) solver->physics;
  MPIVariables *mpi = (MPIVariables *) param->m;
  /* Need to modify the code so that is works for multiple windows */
  /* Assumming m_generator[i]->getSpatialBasis() is called before */
  if (m_basis.size() == 0) {
    if (!m_rank) {
      printf("ERROR in LSROMObject::projectInitialSolution() - m_generator is a vector of size 0.\n");
    }
    return;
  }

  CAROM::Vector* m_fomwork;
  m_fomwork = new CAROM::Vector(a_U.getData(), a_U.dim(), false);
  for (int i = 0; i < m_basis.size(); i++) {
    CAROM::Vector* m_working;
    m_working = new CAROM::Vector(m_rdims[i], false);
    m_projected_init.push_back(new CAROM::Vector(m_rdims[i], false));

    m_working = ProjectToRB(m_fomwork, m_basis[i], m_rdims[i]);
    MPISum_double(m_projected_init[i]->getData(), m_working->getData(), m_rdims[i], &mpi->world);
  if (!m_rank) {
    std::cout << "Checking projected initial: ";
    for (int j = 0; j < m_rdims[i]; j++) {
          std::cout << m_projected_init[i]->item(j) << " ";
    }
    std::cout << std::endl;
  }
  }
}

/*! take a sample (solution snapshot) */
void LSROMObject::takeSample(  const CAROM::Vector& a_U, /*!< solution vector */
                               const double a_time, /*!< sample time */
                               void* a_s )
{
  SimulationObject* sim = (SimulationObject*) a_s;

  HyPar  *solver = (HyPar*) &(sim[0].solver);
  Vlasov *param  = (Vlasov*) solver->physics;
  MPIVariables *mpi = (MPIVariables *) param->m;

  std::vector<double> vec_wo_ghosts(param->npts_local_x, 0.0);
  std::vector<double> vec_x_wghosts(param->npts_local_x_wghosts);
  std::vector<int> index(sim[0].solver.ndims);

  if (m_tic == 0) {

    m_options.push_back(new CAROM::Options(m_vec_size, max_num_snapshots, 1, update_right_SV));
    if (m_energy_criteria > 0){
      m_options[m_curr_win]->setSingularValueTol(m_energy_criteria);
    }

    const std::string basisFileName = basisName + std::to_string(m_parametric_id) + "_" + std::to_string(m_tic);
    m_generator.push_back(new CAROM::BasisGenerator(*m_options[m_curr_win], isIncremental, basisFileName));
    m_intervals.push_back( Interval(a_time, m_t_final) );
    m_ls_is_trained.push_back(false);
    m_snap.push_back(0);

    if (!m_rank) {
      printf( "LSROMObject::takeSample() - creating new generator object for sim. domain %d, var %d, t=%f (total: %d).\n",
              m_sim_idx, m_var_idx, m_intervals[m_curr_win].first, m_generator.size());
    }
    bool addSample = m_generator[m_curr_win]->takeSample( a_U.getData(), a_time, m_dt );
//  char buffer[] = "sample";  // Creates a modifiable buffer and copies the string literal
//  OutputlibROMfield(m_generator[m_curr_win]->getSnapshotMatrix()->getColumn(0)->getData(),
//                    sim[0],
//                    buffer);

    if (m_solve_phi) {
      m_options_phi.push_back(new CAROM::Options(param->npts_local_x, max_num_snapshots, 1, update_right_SV));
      if (m_energy_criteria > 0){
        m_options_phi[m_curr_win]->setSingularValueTol(m_energy_criteria);
      }
      const std::string basisFileName_phi = basisName_phi + std::to_string(m_parametric_id) + "_" + std::to_string(m_tic);
      m_generator_phi.push_back(new CAROM::BasisGenerator(*m_options_phi[m_curr_win], isIncremental, basisFileName_phi));
      if (mpi->ip[1] == 0) {
          ArrayCopynD(1,
                      param->potential,
                      vec_wo_ghosts.data(),
                      sim[0].solver.dim_local,
                      sim[0].solver.ghosts,
                      0,
                      index.data(),
                      sim[0].solver.nvars);
      }
      else {
        vec_wo_ghosts = std::vector<double> (vec_wo_ghosts.size(),0.0);
      }
      bool addphiSample = m_generator_phi[m_curr_win]->takeSample( vec_wo_ghosts.data(), a_time, m_dt );

      m_options_e.push_back(new CAROM::Options(param->npts_local_x, max_num_snapshots, 1, update_right_SV));
      m_generator_e.push_back(new CAROM::BasisGenerator(*m_options_e[m_curr_win], isIncremental, basisName));
      if (mpi->ip[1] == 0) {
      ArrayCopynD(1,
                  param->e_field,
                  vec_wo_ghosts.data(),
                  sim[0].solver.dim_local,
                  sim[0].solver.ghosts,
                  0,
                  index.data(),
                  sim[0].solver.nvars);
      }
      else {
        vec_wo_ghosts = std::vector<double> (vec_wo_ghosts.size(),0.0);
      }
      bool addeSample = m_generator_e[m_curr_win]->takeSample( vec_wo_ghosts.data(), a_time, m_dt );
    }
  } else {

    bool addSample = m_generator[m_curr_win]->takeSample( a_U.getData(), a_time, m_dt );
//  char buffer[] = "sample";  // Creates a modifiable buffer and copies the string literal
//  OutputlibROMfield(m_generator[m_curr_win]->getSnapshotMatrix()->getColumn(0)->getData(),
//                    sim[0],
//                    buffer);
    if (m_solve_phi) {
      if (mpi->ip[1] == 0) {
          ArrayCopynD(1,
                      param->potential,
                      vec_wo_ghosts.data(),
                      sim[0].solver.dim_local,
                      sim[0].solver.ghosts,
                      0,
                      index.data(),
                      sim[0].solver.nvars);
      }
      else {
        vec_wo_ghosts = std::vector<double> (vec_wo_ghosts.size(),0.0);
      }
      bool addphiSample = m_generator_phi[m_curr_win]->takeSample( vec_wo_ghosts.data(), a_time, m_dt );

      if (mpi->ip[1] == 0) {
        ArrayCopynD(1,
                    param->e_field,
                    vec_wo_ghosts.data(),
                    sim[0].solver.dim_local,
                    sim[0].solver.ghosts,
                    0,
                    index.data(),
                    sim[0].solver.nvars);
      }
      else {
        vec_wo_ghosts = std::vector<double> (vec_wo_ghosts.size(),0.0);
      }
      bool addeSample = m_generator_e[m_curr_win]->takeSample( vec_wo_ghosts.data(), a_time, m_dt );
    }

    if (m_tic%m_num_window_samples == 0) {

      m_intervals[m_curr_win].second = a_time;
      int ncol = m_generator[m_curr_win]->getSnapshotMatrix()->numColumns();
      if (!m_rank) {
        printf( "LSROMObject::train() - training LS object %d for sim. domain %d, var %d with %d samples.\n",
                m_curr_win, m_sim_idx, m_var_idx, ncol );
      }
      m_ls_is_trained[m_curr_win] = false;
      m_snap.push_back(0);
      m_curr_win++;

      m_options.push_back(new CAROM::Options(m_vec_size, max_num_snapshots, 1, update_right_SV));
      if (m_energy_criteria > 0){
        m_options[m_curr_win]->setSingularValueTol(m_energy_criteria);
      }
      const std::string basisFileName = basisName + std::to_string(m_parametric_id) + "_" + std::to_string(m_curr_win);
      m_generator.push_back(new CAROM::BasisGenerator(*m_options[m_curr_win], isIncremental, basisFileName));

      bool addSample = m_generator[m_curr_win]->takeSample( a_U.getData(), a_time, m_dt );
      if (m_solve_phi) {
        m_options_phi.push_back(new CAROM::Options(param->npts_local_x, max_num_snapshots, 1, update_right_SV));
        if (m_energy_criteria > 0){
          m_options_phi[m_curr_win]->setSingularValueTol(m_energy_criteria);
        }

        const std::string basisFileName_phi = basisName_phi + std::to_string(m_parametric_id) + "_" + std::to_string(m_curr_win);
        m_generator_phi.push_back(new CAROM::BasisGenerator(*m_options_phi[m_curr_win], isIncremental, basisName_phi));

        if (mpi->ip[1] == 0) {
          ArrayCopynD(1,
                      param->potential,
                      vec_wo_ghosts.data(),
                      sim[0].solver.dim_local,
                      sim[0].solver.ghosts,
                      0,
                      index.data(),
                      sim[0].solver.nvars);
        }
        else {
          vec_wo_ghosts = std::vector<double> (vec_wo_ghosts.size(),0.0);
        }
        bool addphiSample = m_generator_phi[m_curr_win]->takeSample( vec_wo_ghosts.data(), a_time, m_dt );

        m_options_e.push_back(new CAROM::Options(param->npts_local_x, max_num_snapshots, 1, update_right_SV));
        m_generator_e.push_back(new CAROM::BasisGenerator(*m_options_e[m_curr_win], isIncremental, basisName));

        if (mpi->ip[1] == 0) {
          ArrayCopynD(1,
                      param->e_field,
                      vec_wo_ghosts.data(),
                      sim[0].solver.dim_local,
                      sim[0].solver.ghosts,
                      0,
                      index.data(),
                      sim[0].solver.nvars);
        }
        else {
          vec_wo_ghosts = std::vector<double> (vec_wo_ghosts.size(),0.0);
        }
        bool addeSample = m_generator_e[m_curr_win]->takeSample( vec_wo_ghosts.data(), a_time, m_dt );
      }

      m_ls_is_trained.push_back(false);
      m_intervals.push_back( Interval(a_time, m_t_final) );

      if (!m_rank) {
        printf( "LSROMObject::takeSample() - creating new generator object for sim. domain %d, var %d, t=%f (total: %d).\n",
                m_sim_idx, m_var_idx, m_intervals[m_curr_win].first, m_generator.size());
      }
    }
  }
  if (!m_rank) {
    printf("checking m_tic %d\n",m_tic);
  }

  m_tic++;
  return;
}

/*! train the LS objects */
void LSROMObject::train(void* a_s)
{
  /* In the training, following quantities are computed:
   * 1. reduced basis (m_generator[i]->getSpatialBasis)
   * 2. reduced hyperbolic operator (constructHyperbolic)
   * 3. Initial condition (projectInitialCondition) */

  SimulationObject* sim = (SimulationObject*) a_s;
  HyPar  *solver = (HyPar*) &(sim[0].solver);
  Vlasov *param  = (Vlasov*) solver->physics;
  MPIVariables *mpi = (MPIVariables *) param->m;

  if (m_generator.size() > 0) {
    for (int i = 0; i < m_generator.size(); i++) {
      if (!m_ls_is_trained[i]) {
        int ncol = m_generator[i]->getSnapshotMatrix()->numColumns();
        if (!m_rank) {
          printf( "LSROMObject::train() - training LS object %d for sim. domain %d, var %d with %d samples.\n", 
                  m_curr_win, m_sim_idx, m_var_idx, ncol );
        }
        if (m_write_snapshot_mat) {
          /* Since here is using gertSnapshotMatrix to write files, need additional process in order to visualize it
             Similar in the case of DMD */
          char idx_string[_MAX_STRING_SIZE_];
          sprintf(idx_string, "%04d", i);
          std::string fname_root(m_dirname + "/snapshot_mat_"+std::string(idx_string));
          m_generator[i]->getSnapshotMatrix()->write(fname_root);
        }

        /* The code below is like m_dmd[i]->train(m_rdim) but in LS-ROM, there is no
           m_ls object, instead, it has m_generator, m_options, m_spatialbasis, m_projected_init */

        /* IMPORTANT!!! m_generator[i]->getSnapshotMatrix() is modified after
         * getSingularValues or (computeSVD) is called, hence need to make a copy of snapshots */
        /* Does everytime invoking getSpatialBasis() do computeSVD again? */
        m_snapshots.push_back(new CAROM::Matrix(*m_generator[i]->getSnapshotMatrix()));

        m_rdims.push_back(m_rdim);
        if (m_generator[i]->getSnapshotMatrix()->numColumns() < m_rdims[i]) {
          throw std::runtime_error("# of snapshots is less than m_rdim");
        }

        m_S  = m_generator[i]->getSingularValues();

        m_basis.push_back(new CAROM::Matrix(
                              m_generator[i]->getSpatialBasis()->getData(),
                              m_generator[i]->getSpatialBasis()->numRows(),
                              m_generator[i]->getSpatialBasis()->numColumns(),
                              false,
                              true));
        if (m_rdims[i] != m_generator[i]->getSpatialBasis()->numColumns()){
          m_rdims[i] = m_generator[i]->getSpatialBasis()->numColumns();
          if (!m_rank) printf("m_rdim %d is reset to %d \n",
                              m_rdim,m_generator[i]->getSpatialBasis()->numColumns());
        }
        if (!m_rank) {
          std::cout << "----------------------------------------\n";
          std::cout << "Time window #" << i << ": # f POD basis : " << m_rdims[i];
          std::cout << std::endl;
        }

        if (i > 0) {
          m_fullscale.push_back(new CAROM::Matrix(m_rdims[i], m_rdims[i-1], false));
          m_matrix = m_basis[i]->transposeMult(m_basis[i-1]);
          MPISum_double(m_fullscale[i-1]->getData(), m_matrix->getData(), m_rdims[i-1]*m_rdims[i], &mpi->world);
        }

        m_projected_init.push_back(new CAROM::Vector(m_rdims[i], false));
        m_romcoef.push_back(new CAROM::Vector(m_rdims[i], false));

        OutputROMBasis(a_s, m_generator[i]->getSpatialBasis(), i);
        CheckSolProjError(a_s,i);

        if ((!m_solve_phi) && (!m_direct_comp_hyperbolic)) {
          m_romhyperb.push_back(new CAROM::Matrix(m_rdims[i],m_rdims[i],false));
          ConstructROMHy(a_s, m_basis[i], i);
          CheckHyProjError(a_s,i);
        }

        if (m_solve_phi) {
          m_snapshots_phi.push_back(new CAROM::Matrix(
                                    *m_generator_phi[i]->getSnapshotMatrix())),
          m_rdims_phi.push_back(m_rdim_phi);
					// m_snapshots_phi is not distributed due to getSnapshotMatrix
          m_S_phi = m_generator_phi[i]->getSingularValues();
					// m_basis_phi is necessary since getSpatialBasis is distributed
          m_basis_phi.push_back(new CAROM::Matrix(
                                m_generator_phi[i]->getSpatialBasis()->getData(),
                                m_generator_phi[i]->getSpatialBasis()->numRows(),
                                m_generator_phi[i]->getSpatialBasis()->numColumns(),
                                false,
                                true));
          if (m_rdims_phi[i] != m_generator_phi[i]->getSpatialBasis()->numColumns()){
            m_rdims_phi[i] = m_generator_phi[i]->getSpatialBasis()->numColumns();
            if (!m_rank) printf("m_rdim_phi %d is reset to %d \n",
                                m_rdim_phi,
                                m_generator_phi[i]->getSpatialBasis()->numColumns());
          }
          if (!m_rank) {
            std::cout << "----------------------------------------\n";
            std::cout << "Time window #" << i << ": # potential POD basis : " << m_rdims_phi[i];
            std::cout << std::endl;
          }
          m_projected_init_phi.push_back(new CAROM::Vector(m_rdims_phi[i], false));

          OutputROMBasisPhi(a_s, m_generator_phi[i]->getSpatialBasis(),i);

          m_romrhs_phi.push_back(new CAROM::Matrix(m_rdims_phi[i], m_rdims[i], false));
          ConstructPotentialROMRhs(a_s,
                                   m_basis[i],
                                   m_basis_phi[i],
                                   i);

          m_romlaplace_phi.push_back(new CAROM::Matrix(m_rdims_phi[i],m_rdims_phi[i],false));
          ConstructPotentialROMLaplace(a_s, m_basis_phi[i], i);

          CheckPotentialProjError(a_s,i);
          CheckLaplaceProjError(a_s,i);
          CheckRhsProjError(a_s,i);

          m_snapshots_e.push_back(new CAROM::Matrix(
                                  *m_generator_e[i]->getSnapshotMatrix()));
          m_basis_e.push_back(new CAROM::Matrix(
                              m_generator_phi[i]->getSpatialBasis()->numRows(),
                              m_rdims_phi[i], false));
          ConstructEBasis(a_s,i);
          CheckEProjError(a_s,i);

          m_romhyperb_x.push_back(new CAROM::Matrix(m_rdims[i], m_rdims[i],false));
          ConstructROMHy_x(a_s, m_basis[i], i);
          m_romhyperb_v.push_back(std::vector<CAROM::Matrix*>());
          ConstructROMHy_v(a_s, m_basis[i], m_basis_e[i], i);
        }

      }
      if (!m_rank) {
        std::cout << "----------------\n";
        std::cout << "Singular Values of f: ";
        for (int i = 0; i < m_S->dim(); i++) {
              std::cout << (m_S->item(i)) << " ";
        }
        std::cout << "\n";
        std::cout << std::endl;
      }
      if (m_solve_phi) {
        if (!m_rank) {
          std::cout << "----------------\n";
          std::cout << "Singular Values for potential: ";
          for (int i = 0; i < m_S_phi->dim(); i++) {
                std::cout << (m_S_phi->item(i)) << " ";
          }
          std::cout << "\n";
          std::cout << std::endl;
        }
      }
    }
  } else {
    printf("ERROR in LSROMObject::train(): m_generator is of size zero!");
  }

  return;
}

/*! compute prediction at given time, the ROM coefficient is stored in m_romcoef  */
const CAROM::Vector* LSROMObject::predict(const double a_t, /*!< time at which to predict solution */
                                          void* a_s )
{
  SimulationObject* sim = (SimulationObject*) a_s;
  HyPar  *solver = (HyPar*) &(sim[0].solver);
  Vlasov *param  = (Vlasov*) solver->physics;
  MPIVariables *mpi = (MPIVariables *) param->m;

  for (int i = 0; i < m_rdims.size(); i++) {
    if (   (a_t >= m_intervals[i].first)
        && (  (a_t < m_intervals[i].second) || (m_intervals[i].second < 0)  ) ){

      if (!m_rank) printf("LS-ROM # %d predicts at time t = %f in interval [%f, %f] \n",
                          i,a_t,m_intervals[i].first,m_intervals[i].second);

      std::vector<int> index(sim[0].solver.ndims);
      std::vector<double> vec_wghosts(sim[0].solver.npoints_local_wghosts*sim[0].solver.nvars);
      std::vector<double> rhs_wghosts(sim[0].solver.npoints_local_wghosts*sim[0].solver.nvars);
      int num_rows = m_basis[i]->numRows();

      CAROM::Vector* m_fomwork;
      m_fomwork = new CAROM::Vector(num_rows,false);
      CAROM::Vector* m_working;
      m_working = new CAROM::Vector(m_rdims[i], false);

      if(std::abs(a_t) < 1e-9) {
        m_working = ProjectToRB(m_snapshots[i]->getColumn(0), m_basis[i], m_rdims[i]);
        MPISum_double(m_projected_init[i]->getData(), m_working->getData(), m_rdims[i], &mpi->world);
        m_romcoef[i] = m_projected_init[i];
        m_fomwork = ReconlibROMfield(m_romcoef[i], m_basis[i], m_rdims[i]);

        ArrayCopynD(sim[0].solver.ndims,
                    m_fomwork->getData(),
                    vec_wghosts.data(),
                    sim[0].solver.dim_local,
                    0,
                    sim[0].solver.ghosts,
                    index.data(),
                    sim[0].solver.nvars);
        sim[0].solver.PostStage( vec_wghosts.data(),
                                 &(sim[0].solver),
                                 &(sim[0].mpi),
                                 0); CHECKERR(ierr);

        /* Setup RK44 parameters */
        TimeExplicitRKInitialize();
        /* Initialize RK44 working variables */
        TimeInitialize(m_rdims[i]);
      } else {
        if ((m_snap[i]==0) && (i>0)) {
          m_fomwork = ReconlibROMfield(m_romcoef[i-1], m_basis[i-1], m_rdims[i-1]);
          m_working = ProjectToRB(m_fomwork, m_basis[i], m_rdims[i]);
          MPISum_double(m_romcoef[i]->getData(), m_working->getData(), m_rdims[i], &mpi->world);
          TimeInitialize(m_rdims[i]);
      } else {
        }
        TimeRK(a_t,a_s,i);
      }

      if ((m_snap[i] < sim[0].solver.n_iter) && (m_snap[i] % m_sampling_freq == 0)){
        int idx;
        idx = m_snap[i]/m_sampling_freq;
        if (!m_rank) printf("idx %d m_snap %d m_sampling_freq %d \n",idx,m_snap[i],m_sampling_freq);
        m_fomwork = ReconlibROMfield(m_romcoef[i], m_basis[i], m_rdims[i]);
        ArrayCopynD(sim[0].solver.ndims,
                    m_fomwork->getData(),
                    vec_wghosts.data(),
                    sim[0].solver.dim_local,
                    0,
                    sim[0].solver.ghosts,
                    index.data(),
                    sim[0].solver.nvars);
        ArrayCopynD(sim[0].solver.ndims,
                    m_snapshots[i]->getColumn(idx)->getData(),
                    rhs_wghosts.data(),
                    sim[0].solver.dim_local,
                    0,
                    sim[0].solver.ghosts,
                    index.data(),
                    sim[0].solver.nvars);
        char buffer[] = "reproderr";  // Creates a modifiable buffer and copies the string literal
        CalSnapROMDiff(&(sim[0].solver),&(sim[0].mpi),rhs_wghosts.data(),vec_wghosts.data(),buffer);
        if (!m_rank) printf("Reconstructive error at # %d snapshot, %.15f %.15f %.15f \n",
                            idx,
                            sim[0].solver.rom_diff_norms[0],
                            sim[0].solver.rom_diff_norms[1],
                            sim[0].solver.rom_diff_norms[2]);
      }
      m_snap[i]++;
      return ReconlibROMfield(m_romcoef[i], m_basis[i], m_rdims[i]);
    }
  }
  printf("ERROR in LSROMObject::predict(): m_generator is of size zero or interval not found!");
  return nullptr;
}

void LSROMObject::save(const std::string& a_fname_root /*!< Filename root */) const
{
    return;
}
void LSROMObject::load(const std::string& a_fname_root /*!< Filename root */)
{
    return;
}

int LSROMObject::TimeInitialize(int a_rdim)
{
  /* Currenty assuming one window only */
  int i;
  /* initialize arrays to NULL, then allocate as necessary */
  m_U.clear();
  m_Udot.clear();

  /* explicit Runge-Kutta methods */
  for (i = 0; i < nstages; i++) {
    m_U.push_back(new CAROM::Vector(a_rdim, false));
    m_Udot.push_back(new CAROM::Vector(a_rdim, false));
  }
  return 0;
}

int LSROMObject::TimeExplicitRKInitialize()
{
  /* Currenty assuming one window only */
  /* Currently only support RK44 */
  nstages = 4;
  A = (double*) calloc (nstages*nstages,sizeof(double));
  b = (double*) calloc (nstages        ,sizeof(double));
  c = (double*) calloc (nstages        ,sizeof(double));
  _ArraySetValue_(A,nstages*nstages,0.0);
  _ArraySetValue_(b,nstages        ,0.0);
  _ArraySetValue_(c,nstages        ,0.0);
  A[4] = 0.5; A[9] = 0.5; A[14] = 1.0;
  c[1] = c[2] = 0.5; c[3] = 1.0;
  b[0] = 1.0/6.0; b[1] = 1.0/3.0; b[2] = 1.0/3.0; b[3] = 1.0/6.0;

  return 0;
}

int LSROMObject::TimeRK(const double a_t, /*!< time at which to predict solution */
                        void* a_s,
                        int idx
                        )
{
  /* Advance the ROM ODE using RK4 scheme */

  SimulationObject* sim = (SimulationObject*) a_s;
  HyPar  *solver = (HyPar*) &(sim[0].solver);
  Vlasov *param  = (Vlasov*) solver->physics;
  MPIVariables *mpi = (MPIVariables *) param->m;

  std::vector<int> index(sim[0].solver.ndims);
  std::vector<double> vec_wghosts(sim[0].solver.npoints_local_wghosts*sim[0].solver.nvars);
  std::vector<double> rhs_wghosts(sim[0].solver.npoints_local_wghosts*sim[0].solver.nvars);

  int ns, stage, i;
  int num_rows = m_basis[idx]->numRows();

  CAROM::Vector* m_fomwork;
  CAROM::Vector* m_rhswork;
  CAROM::Vector* m_romwork;

  m_fomwork = new CAROM::Vector(num_rows,false);
  m_rhswork = new CAROM::Vector(num_rows,false);
  m_romwork = new CAROM::Vector(m_rdims[idx],false);

  CAROM::Vector* m_tmprhs;
  CAROM::Vector* m_tmpsol;
  CAROM::Vector* m_e;
  CAROM::Vector* m_contract1;
  CAROM::Vector* m_contract2;
  if (m_solve_phi) {
    m_contract1 = new CAROM::Vector(m_rdims_phi[idx],false);
    m_contract2 = new CAROM::Vector(m_rdims[idx],false);
  }

    /* Calculate stage values */
  for (stage = 0; stage < nstages; stage++) {
  
    double stagetime = a_t + c[stage]*m_dt;

    _ArrayCopy1D_(  m_romcoef[idx]->getData(),
                    m_U[stage]->getData(),
                    m_rdims[idx] );

    for (i = 0; i < stage; i++) {
      _ArrayAXPY_(  m_Udot[i]->getData(),
                    (m_dt * A[stage*nstages+i]),
                    m_U[stage]->getData(),
                    m_rdims[idx] );
    }

    if (m_direct_comp_hyperbolic) {
//    printf("compute hyperbolic term directly\n");

      m_fomwork = ReconlibROMfield(m_U[stage], m_basis[idx], m_rdims[idx]);

      ArrayCopynD(sim[0].solver.ndims,
                  m_fomwork->getData(),
                  vec_wghosts.data(),
                  sim[0].solver.dim_local,
                  0,
                  sim[0].solver.ghosts,
                  index.data(),
                  sim[0].solver.nvars);

      sim[0].solver.PostStage( vec_wghosts.data(),
                               &(sim[0].solver),
                               &(sim[0].mpi),
                               stagetime); CHECKERR(ierr);

      /* Evaluate F(\phi_j) */
      TimeRHSFunctionExplicit(rhs_wghosts.data(),
                              vec_wghosts.data(),
                              &(sim[0].solver),
                              &(sim[0].mpi),
                              0);



      ArrayCopynD(sim[0].solver.ndims,
                  rhs_wghosts.data(),
                  m_rhswork->getData(),
                  sim[0].solver.dim_local,
                  sim[0].solver.ghosts,
                  0,
                  index.data(),
                  sim[0].solver.nvars);

      m_romwork = ProjectToRB(m_rhswork, m_basis[idx], m_rdims[idx]);
      MPISum_double(m_Udot[stage]->getData(), m_romwork->getData(), m_rdims[idx], &mpi->world);

//    if (!m_rank) {
//      std::cout << "Checking hyperbolic term [directly]: ";
//      for (int j = 0; j < m_rdim; j++) {
//            std::cout << m_Udot[stage]->item(j) << " ";
//      }
//      std::cout << std::endl;
//    }
    }
    else if (m_solve_phi) {
      m_tmprhs = m_romrhs_phi[idx]->mult(m_U[stage]);
      m_tmpsol = m_romlaplace_phi[idx]->mult(m_tmprhs);

      m_romwork = m_romhyperb_x[idx]->mult(m_U[stage]);

//    /* Tensor contraction */
      for (int k = 0; k < m_rdims[idx]; k++) {
        m_contract1 = m_romhyperb_v[idx][k]->mult(m_U[stage]);
        m_contract2->item(k) = m_contract1->inner_product(m_tmpsol);
      }
      m_Udot[stage] = m_romwork->plus(m_contract2);
//
//    m_e = m_basis_e->mult(m_tmpsol);

//    ArrayCopynD(1,
//                m_e->getData(),
//                param->e_field,
//                sim[0].solver.dim_local,
//                0,
//                sim[0].solver.ghosts,
//                index.data(),
//                sim[0].solver.nvars);

//    m_fomwork = ReconlibROMfield(m_U[stage], m_generator[0]->getSpatialBasis(), m_rdim);

//    ArrayCopynD(sim[0].solver.ndims,
//                m_fomwork->getData(),
//                vec_wghosts.data(),
//                sim[0].solver.dim_local,
//                0,
//                sim[0].solver.ghosts,
//                index.data(),
//                sim[0].solver.nvars);
//      solver->ApplyBoundaryConditions(solver,mpi,vec_wghosts.data(),NULL,0);
//      solver->ApplyIBConditions(solver,mpi,vec_wghosts.data(),0);
//      MPIExchangeBoundariesnD(  solver->ndims,
//                                solver->nvars,
//                                solver->dim_local,
//                                solver->ghosts,
//                                mpi,
//                                vec_wghosts.data());

//      HyperbolicFunction_1dir( rhs_wghosts.data(),
//                               vec_wghosts.data(),
//                               solver,
//                               mpi,
//                               0,
//                               1,
//                               solver->FFunction,
//                               solver->Upwind, 1 );

//    ArrayCopynD(sim[0].solver.ndims,
//                rhs_wghosts.data(),
//                m_rhswork->getData(),
//                sim[0].solver.dim_local,
//                sim[0].solver.ghosts,
//                0,
//                index.data(),
//                sim[0].solver.nvars);

//    m_contract2 = ProjectToRB(m_rhswork,m_generator[0]->getSpatialBasis(), m_rdim);
//    m_Udot[stage] = m_romwork->plus(m_contract2);
      *m_Udot[stage] *= -1.0;

      /* Reconstruct potential */
//    m_e = m_basis_e->mult(m_tmpsol);

//    ArrayCopynD(1,
//                m_e->getData(),
//                param->e_field,
//                sim[0].solver.dim_local,
//                0,
//                sim[0].solver.ghosts,
//                index.data(),
//                sim[0].solver.nvars);

//    m_fomwork = ReconlibROMfield(m_U[stage], m_generator[0]->getSpatialBasis(), m_rdim);

//    ArrayCopynD(sim[0].solver.ndims,
//                m_fomwork->getData(),
//                vec_wghosts.data(),
//                sim[0].solver.dim_local,
//                0,
//                sim[0].solver.ghosts,
//                index.data(),
//                sim[0].solver.nvars);

//    /* Evaluate F(\phi_j) */
//    TimeRHSFunctionExplicit(rhs_wghosts.data(),
//                            vec_wghosts.data(),
//                            &(sim[0].solver),
//                            &(sim[0].mpi),
//                            0);

//    ArrayCopynD(sim[0].solver.ndims,
//                rhs_wghosts.data(),
//                m_rhswork->getData(),
//                sim[0].solver.dim_local,
//                sim[0].solver.ghosts,
//                0,
//                index.data(),
//                sim[0].solver.nvars);

//    m_Udot[stage] = ProjectToRB(m_rhswork,m_generator[0]->getSpatialBasis(), m_rdim);
//    if (!m_rank) {
//      std::cout << "Checking hyperbolic term [phi]: ";
//      for (int j = 0; j < m_rdim; j++) {
//            std::cout << m_Udot[stage]->item(j) << " ";
//      }
//      std::cout << std::endl;
//    }
    }
    else {
      m_Udot[stage] = m_romhyperb[idx]->mult(m_U[stage]);
//    if (!m_rank) {
//      std::cout << "Checking hyperbolic term [efficient]: ";
//      for (int j = 0; j < m_rdim; j++) {
//            std::cout << m_Udot[stage]->item(j) << " ";
//      }
//      std::cout << std::endl;
//    }
    }

  }

  /* Step completion */
  for (stage = 0; stage < nstages; stage++) {

    _ArrayAXPY_(  m_Udot[stage]->getData(),
                  (m_dt * b[stage]),
                  m_romcoef[idx]->getData(),
                  m_rdims[idx] );

  }
  if (!m_rank) {
    std::cout << "Checking solved coefficient: ";
    for (int j = 0; j < m_rdims[idx]; j++) {
          std::cout << m_romcoef[idx]->item(j) << " ";
    }
    std::cout << std::endl;
  }
  delete m_fomwork;
  delete m_rhswork;
  delete m_romwork;
  if (m_solve_phi) {
    delete m_tmprhs;
    delete m_tmpsol;
    delete m_e;
    delete m_contract1;
    delete m_contract2;
  }
  return 0;
}

void LSROMObject::OutputlibROMfield(double* vec_data, SimulationObject& a_s, char* filename) {
  WriteArray( a_s.solver.ndims,
              a_s.solver.nvars,
              a_s.solver.dim_global,
              a_s.solver.dim_local,
              0,
              a_s.solver.x,
              vec_data,
              &(a_s.solver),
              &(a_s.mpi),
              filename);
}

CAROM::Vector* LSROMObject::ReconlibROMfield(const CAROM::Vector* a_romcoef, const CAROM::Matrix* a_rombasis, const int a_rdim){

  return a_rombasis->getFirstNColumns(a_rdim)->mult(a_romcoef);
}

CAROM::Vector* LSROMObject::ProjectToRB(CAROM::Vector* a_field, const CAROM::Matrix* a_rombasis, const int a_rdim){
  /* a_field is a serialized solution without ghost points */

  return a_rombasis->getFirstNColumns(a_rdim)->transposeMult(a_field);
}

/*! Check LS objects */
void LSROMObject::check(void* a_s)
{
  SimulationObject* sim = (SimulationObject*) a_s;
        /* Below is just checking if the basis function is orthogonal */
//      projectInitialSolution(*(m_generator[0]->getSpatialBasis()->getColumn(0)));
//      projectInitialSolution(*(m_generator[0]->getSpatialBasis()->getColumn(1)));
//      projectInitialSolution(*(m_generator[0]->getSpatialBasis()->getColumn(2)));
//      projectInitialSolution(*(m_generator[0]->getSpatialBasis()->getColumn(3)));

        /* Compute initial condition */
        /* Can't simply pass *(m_generator[0]->getSnapshotMatrix()->getColumn(0))) since
         * m_generator[0]->getSpatialBasis() is distributed whereas m_generator[0]->getSnapshotMatrix()
         * is not distributed */

//      projectInitialSolution(*(m_snapshots[i]->getColumn(0)));
        if (!m_rank) {
          printf("m_project size %d\n",m_projected_init[0]->dim());
          std::cout << "Checking initial condition: ";
          for (int j = 0; j < m_projected_init[0]->dim(); j++) {
                std::cout << (m_projected_init[0]->item(j)) << " ";
          }
          std::cout << std::endl;
        }
        int num_rows = m_generator[0]->getSpatialBasis()->numRows();

        CAROM::Vector* recon_init;
        recon_init = new CAROM::Vector(num_rows,false);
        recon_init = m_generator[0]->getSpatialBasis()->getFirstNColumns(m_rdim)->mult(m_projected_init[0]);

//      std::vector<int> index(sim[0].solver.ndims);
//      ArrayCopynD(sim[0].solver.ndims,
//                  recon_init->getData(),
//                  sim[0].solver.u_rom_predicted,
//                  sim[0].solver.dim_local,
//                  0,
//                  sim[0].solver.ghosts,
//                  index.data(),
//                  sim[0].solver.nvars);
//        const char* fnam2 = "recon_";
//        std::string fname2 = std::string(fnam2) + std::to_string(0);
//        char* fname_buffer2 = new char[fname2.length() + 1];
//        std::strcpy(fname_buffer2, fname2.c_str());

//        WriteArray( sim[0].solver.ndims,
//                    sim[0].solver.nvars,
//                    sim[0].solver.dim_global,
//                    sim[0].solver.dim_local,
//                    0,
//                    sim[0].solver.x,
//                    recon_init->getData(),
//                    &(sim[0].solver),
//                    &(sim[0].mpi),
//                    fname_buffer2);
//      ArrayCopynD(sim[0].solver.ndims,
//                  m_snapshots->getColumn(0)->getData(),
//                  sim[0].solver.u,
//                  sim[0].solver.dim_local,
//                  0,
//                  sim[0].solver.ghosts,
//                  index.data(),
//                  sim[0].solver.nvars);
//        const char* fnam3 = "snap_col_";
//        std::string fname3 = std::string(fnam3) + std::to_string(0);
//        char* fname_buffer3 = new char[fname3.length() + 1];
//        std::strcpy(fname_buffer3, fname3.c_str());

//        double* vec_data = m_generator[0]->getSnapshotMatrix()->getColumn(0)->getData();
//        WriteArray( sim[0].solver.ndims,
//                    sim[0].solver.nvars,
//                    sim[0].solver.dim_global,
//                    sim[0].solver.dim_local,
//                    0,
//                    sim[0].solver.x,
//                    m_snapshots->getColumn(0)->getData(),
//                    &(sim[0].solver),
//                    &(sim[0].mpi),
//                    fname_buffer3);
//      CalculateROMDiff(  &(sim[0].solver),
//                         &(sim[0].mpi) );

//      projectInitialSolution(*(m_snapshots->getColumn(0)));
//      for (int j = 0; j < 1; j++){
//        /* Extend reduced basis \phi_j with ghost points */
//        std::vector<int> index(sim[0].solver.ndims);
//        ArrayCopynD(sim[0].solver.ndims,
//                    m_snapshots->getColumn(j)->getData(),
//                    vec_wghosts.data(),
//                    sim[0].solver.dim_local,
//                    0,
//                    sim[0].solver.ghosts,
//                    index.data(),
//                    sim[0].solver.nvars);

//        /* Evaluate F(\phi_j) */
//        TimeRHSFunctionExplicit(rhs_wghosts.data(),
//                                vec_wghosts.data(),
//                                &(sim[0].solver),
//                                &(sim[0].mpi),
//                                0);

//        ArrayCopynD(sim[0].solver.ndims,
//                    rhs_wghosts.data(),
//                    sim[0].solver.u,
//                    sim[0].solver.dim_local,
//                    sim[0].solver.ghosts,
//                    sim[0].solver.ghosts,
//                    index.data(),
//                    sim[0].solver.nvars);

//        CAROM::Vector* recon_init;
//        recon_init = new CAROM::Vector(num_rows,false);
//        CAROM::Vector* phi_colwork;
//        phi_colwork = new CAROM::Vector(num_rows,true);
//        /* Remove ghosts point in F(phi_j) */
//        ArrayCopynD(sim[0].solver.ndims,
//                    rhs_wghosts.data(),
//                    phi_colwork->getData(),
//                    sim[0].solver.dim_local,
//                    sim[0].solver.ghosts,
//                    0,
//                    index.data(),
//                    sim[0].solver.nvars);
//        printf("%d \n",phi_colwork->dim());
//        CAROM::Vector* m_working;
//        m_working = m_generator[0]->getSpatialBasis()->getFirstNColumns(m_rdim)->transposeMult(phi_colwork);
//        printf("%d \n",m_working->dim());
//        recon_init = m_generator[0]->getSpatialBasis()->getFirstNColumns(m_rdim)->mult(m_working);
//        ArrayCopynD(sim[0].solver.ndims,
//                    recon_init->getData(),
//                    sim[0].solver.u_rom_predicted,
//                    sim[0].solver.dim_local,
//                    0,
//                    sim[0].solver.ghosts,
//                    index.data(),
//                    sim[0].solver.nvars);
//        recon_init = m_generator[0]->getSpatialBasis()->getFirstNColumns(m_rdim)->mult(m_projected_init[0]);

//        ArrayCopynD(sim[0].solver.ndims,
//                    recon_init->getData(),
//                    vec_wghosts.data(),
//                    sim[0].solver.dim_local,
//                    0,
//                    sim[0].solver.ghosts,
//                    index.data(),
//                    sim[0].solver.nvars);

//        /* Evaluate F(\phi_j) */
//        TimeRHSFunctionExplicit(rhs_wghosts.data(),
//                                vec_wghosts.data(),
//                                &(sim[0].solver),
//                                &(sim[0].mpi),
//                                0);

//        ArrayCopynD(sim[0].solver.ndims,
//                    rhs_wghosts.data(),
//                    sim[0].solver.u_rom_predicted,
//                    sim[0].solver.dim_local,
//                    sim[0].solver.ghosts,
//                    sim[0].solver.ghosts,
//                    index.data(),
//                    sim[0].solver.nvars);

//        CAROM::Vector* m_working;
//        m_working=m_romhyperb->mult(m_projected_init[0]);
//        recon_init = m_generator[0]->getSpatialBasis()->getFirstNColumns(m_rdim)->mult(m_working);
//        ArrayCopynD(sim[0].solver.ndims,
//                    recon_init->getData(),
//                    sim[0].solver.u_rom_predicted,
//                    sim[0].solver.dim_local,
//                    0,
//                    sim[0].solver.ghosts,
//                    index.data(),
//                    sim[0].solver.nvars);
//        CAROM::Vector* m_working;
//        recon_init = phi_hyper->mult(m_projected_init[0]);
//        m_working = m_generator[0]->getSpatialBasis()->getFirstNColumns(m_rdim)->transposeMult(recon_init);
//        recon_init = m_generator[0]->getSpatialBasis()->getFirstNColumns(m_rdim)->mult(m_working);
//        ArrayCopynD(sim[0].solver.ndims,
//                    recon_init->getData(),
//                    sim[0].solver.u_rom_predicted,
//                    sim[0].solver.dim_local,
//                    0,
//                    sim[0].solver.ghosts,
//                    index.data(),
//                    sim[0].solver.nvars);
//      }
}

/*! Construct reduced hyperbolic operator */
void LSROMObject::ConstructROMHy(void* a_s, const CAROM::Matrix* a_rombasis, int idx)
{
  SimulationObject* sim = (SimulationObject*) a_s;
  HyPar  *solver = (HyPar*) &(sim[0].solver);
  Vlasov *param  = (Vlasov*) solver->physics;
  MPIVariables *mpi = (MPIVariables *) param->m;

  int num_rows = a_rombasis->numRows();
  int num_cols = a_rombasis->numColumns();

  /* Compute F(\Phi), where \Phi is the reduced basis matrix */
  CAROM::Matrix* phi_hyper;
  phi_hyper = new CAROM::Matrix(num_rows, m_rdims[idx], false);
  CAROM::Matrix* m_working;
  m_working = new CAROM::Matrix(m_rdims[idx], m_rdims[idx], false);

  CAROM::Vector phi_hyper_col(num_rows,false);
  std::vector<double> vec_wghosts(sim[0].solver.npoints_local_wghosts*sim[0].solver.nvars);
  std::vector<double> rhs_wghosts(sim[0].solver.npoints_local_wghosts*sim[0].solver.nvars);

  for (int j = 0; j < m_rdims[idx]; j++){
    /* Extend reduced basis \phi_j with ghost points */
    std::vector<int> index(sim[0].solver.ndims);
    ArrayCopynD(sim[0].solver.ndims,
                a_rombasis->getColumn(j)->getData(),
                vec_wghosts.data(),
                sim[0].solver.dim_local,
                0,
                sim[0].solver.ghosts,
                index.data(),
                sim[0].solver.nvars);

    sim[0].solver.PostStage( vec_wghosts.data(),
                             &(sim[0].solver),
                             &(sim[0].mpi),
                             0); CHECKERR(ierr);

    /* Evaluate F(\phi_j) */
    TimeRHSFunctionExplicit(rhs_wghosts.data(),
                            vec_wghosts.data(),
                            &(sim[0].solver),
                            &(sim[0].mpi),
                            0);

    /* Remove ghosts point in F(phi_j) */
    ArrayCopynD(sim[0].solver.ndims,
                rhs_wghosts.data(),
                phi_hyper_col.getData(),
                sim[0].solver.dim_local,
                sim[0].solver.ghosts,
                0,
                index.data(),
                sim[0].solver.nvars);

    /* Copy F(phi_j) back to columns of phi_hyper matrix */
    for (int i = 0; i < num_rows; i++) {
      (*phi_hyper)(i, j) = phi_hyper_col.getData()[i];
    }

    char buffer[] = "hyperbasis";  // Creates a modifiable buffer and copies the string literal
    OutputlibROMfield(phi_hyper->getColumn(j)->getData(),
                      sim[0],
                      buffer);
    /* increment the index string, if required */
    if ((!strcmp(sim[0].solver.output_mode,"serial")) && (!strcmp(sim[0].solver.op_overwrite,"no"))) {
        ::IncrementFilenameIndex(sim[0].solver.filename_index,sim[0].solver.index_length);
    }
  }

  // construct hyper_ROM = phi^T phi_hyper
  m_working = a_rombasis->getFirstNColumns(m_rdims[idx])->transposeMult(phi_hyper);
  MPISum_double(m_romhyperb[idx]->getData(), m_working->getData(), m_rdims[idx]*m_rdims[idx], &mpi->world);
  if (!m_rank){
    printf("phi %d %d\n",a_rombasis->numRows(),a_rombasis->numColumns());
    printf("phi_hyper %d %d\n",phi_hyper->numRows(),phi_hyper->numColumns());
    printf("m_romhyperb %d %d\n",m_romhyperb[idx]->numRows(),m_romhyperb[idx]->numColumns());
  }

  ::ResetFilenameIndex( sim[0].solver.filename_index,
                        sim[0].solver.index_length );
  delete phi_hyper;

  return;
}

/*! Dump ROM basis */
void LSROMObject::OutputROMBasis(void* a_s, const CAROM::Matrix* a_rombasis, int idx)
{
  /* m_rdim is written out */
  SimulationObject* sim = (SimulationObject*) a_s;

  for (int j = 0; j < m_rdims[idx]; j++){
    char buffer[] = "basis";  // Creates a modifiable buffer and copies the string literal
    char fbuffer[100];
    sprintf(fbuffer, "%s_%d",buffer,idx);
    OutputlibROMfield(a_rombasis->getColumn(j)->getData(),
                      sim[0],
                      fbuffer);
    /* increment the index string, if required */
    if ((!strcmp(sim[0].solver.output_mode,"serial")) && (!strcmp(sim[0].solver.op_overwrite,"no"))) {
        ::IncrementFilenameIndex(sim[0].solver.filename_index,sim[0].solver.index_length);
    }
  }

  ::ResetFilenameIndex( sim[0].solver.filename_index,
                        sim[0].solver.index_length );
  return;
}

/*! Calculates the L1, L2, & Linf norms of the diff between
 * snapshots and the reproduced solution by libROM.
*/
int LSROMObject::CalSnapROMDiff( void *s, /*!< Solver object of type #HyPar */
                                 void *m,  /*!< MPI object of type #MPIVariables */
                                 double* a_snapshot,
                                 double* a_romsol,
                                 char* filename)

{
  HyPar* solver = (HyPar*) s;
  MPIVariables* mpi = (MPIVariables*) m;
  double sum = 0, global_sum = 0;

  static const double tolerance = 1e-15;

  int size = solver->npoints_local_wghosts * solver->nvars;
  double* u_diff = (double*) calloc (size, sizeof(double));

  /* calculate solution norms (for relative error) */
  double solution_norm[3] = {0.0,0.0,0.0};
  /* L1 */
  sum = ArraySumAbsnD ( solver->nvars,
                        solver->ndims,
                        solver->dim_local,
                        solver->ghosts,
                        solver->index,
                        a_snapshot );
  global_sum = 0; MPISum_double(&global_sum,&sum,1,&mpi->world);
  solution_norm[0] = global_sum/((double)solver->npoints_global);
  /* L2 */
  sum = ArraySumSquarenD  ( solver->nvars,
                            solver->ndims,
                            solver->dim_local,
                            solver->ghosts,
                            solver->index,
                            a_snapshot );
  global_sum = 0; MPISum_double(&global_sum,&sum,1,&mpi->world);
  solution_norm[1] = sqrt(global_sum/((double)solver->npoints_global));
  /* Linf */
  sum = ArrayMaxnD  ( solver->nvars,
                      solver->ndims,
                      solver->dim_local,
                      solver->ghosts,
                      solver->index,
                      a_snapshot  );
  global_sum = 0; MPIMax_double(&global_sum,&sum,1,&mpi->world);
  solution_norm[2] = global_sum;

  /* set u_diff to PDE solution */
  _ArrayCopy1D_(a_snapshot, u_diff, size);
  /* subtract the ROM solutions */
  _ArrayAXPY_(a_romsol, -1.0, u_diff, size);

  /* calculate L1 norm of error */
  sum = ArraySumAbsnD ( solver->nvars,
                        solver->ndims,
                        solver->dim_local,
                        solver->ghosts,
                        solver->index,
                        u_diff  );
  global_sum = 0; MPISum_double(&global_sum,&sum,1,&mpi->world);
  solver->rom_diff_norms[0] = global_sum/((double)solver->npoints_global);

  /* calculate L2 norm of error */
  sum = ArraySumSquarenD  ( solver->nvars,
                            solver->ndims,
                            solver->dim_local,
                            solver->ghosts,
                            solver->index,
                            u_diff  );
  global_sum = 0; MPISum_double(&global_sum,&sum,1,&mpi->world);
  solver->rom_diff_norms[1] = sqrt(global_sum/((double)solver->npoints_global));

  /* calculate Linf norm of error */
  sum = ArrayMaxnD  ( solver->nvars,
                      solver->ndims,
                      solver->dim_local,
                      solver->ghosts,
                      solver->index,
                      u_diff  );
  global_sum = 0; MPIMax_double(&global_sum,&sum,1,&mpi->world);
  solver->rom_diff_norms[2] = global_sum;

  /* decide whether to normalize and report relative diff norms,
    or report absolute diff norms. */
  if (    (solution_norm[0] > tolerance)
      &&  (solution_norm[1] > tolerance)
      &&  (solution_norm[2] > tolerance) ) {
    solver->rom_diff_norms[0] /= solution_norm[0];
    solver->rom_diff_norms[1] /= solution_norm[1];
    solver->rom_diff_norms[2] /= solution_norm[2];
  }
//WriteArray( solver->ndims,
//            solver->nvars,
//            solver->dim_global,
//            solver->dim_local,
//            solver->ghosts,
//            solver->x,
//            u_diff,
//            solver,
//            mpi,
//            filename);

  free(u_diff);
  return 0;
}

/*! Check projection error in solution */
void LSROMObject::CheckSolProjError(void* a_s, int idx)
{
  if (!m_rank) {
    std::cout << "----------------------------------------\n";
    std::cout << "Checking projection error in snapshots: ";
    std::cout << std::endl;
  }

  SimulationObject* sim = (SimulationObject*) a_s;
  HyPar  *solver = (HyPar*) &(sim[0].solver);
  Vlasov *param  = (Vlasov*) solver->physics;
  MPIVariables *mpi = (MPIVariables *) param->m;

  std::vector<double> snap_wghosts(sim[0].solver.npoints_local_wghosts*sim[0].solver.nvars);
  std::vector<double> snapproj_wghosts(sim[0].solver.npoints_local_wghosts*sim[0].solver.nvars);
  std::vector<int> index(sim[0].solver.ndims);

  int num_rows = m_generator[idx]->getSpatialBasis()->numRows();
  int num_cols = m_snapshots[idx]->numColumns();
  CAROM::Vector* recon_init;
  recon_init = new CAROM::Vector(num_rows,false);
  CAROM::Vector* m_working;
  m_working = new CAROM::Vector(m_rdims[idx], false);

  for (int j = 0; j < num_cols; j++){
    m_working = ProjectToRB(m_snapshots[idx]->getColumn(j), m_basis[idx], m_rdims[idx]);
    MPISum_double(m_projected_init[idx]->getData(), m_working->getData(), m_rdims[idx], &mpi->world);
    if (!m_rank) {
      printf("%d m_project size %d\n",idx,m_projected_init[idx]->dim());
      std::cout << "Checking projected coefficients: ";
      for (int k = 0; k < m_projected_init[idx]->dim(); k++) {
            std::cout << (m_projected_init[idx]->item(k)) << " ";
      }
      std::cout << std::endl;
    }
    /* Extend reduced basis \phi_j with ghost points */
    ArrayCopynD(sim[0].solver.ndims,
                m_snapshots[idx]->getColumn(j)->getData(),
                snap_wghosts.data(),
                sim[0].solver.dim_local,
                0,
                sim[0].solver.ghosts,
                index.data(),
                sim[0].solver.nvars);

    recon_init = m_basis[idx]->getFirstNColumns(m_rdims[idx])->mult(m_projected_init[idx]);
    ArrayCopynD(sim[0].solver.ndims,
                recon_init->getData(),
                snapproj_wghosts.data(),
                sim[0].solver.dim_local,
                0,
                sim[0].solver.ghosts,
                index.data(),
                sim[0].solver.nvars);
      char buffer[] = "snapprojerr";  // Creates a modifiable buffer and copies the string literal
//    CalSnapROMDiff(&(sim[0].solver),&(sim[0].mpi),snap_wghosts.data(),snapproj_wghosts.data(),buffer);
      /* increment the index string, if required */
      if ((!strcmp(sim[0].solver.output_mode,"serial")) && (!strcmp(sim[0].solver.op_overwrite,"no"))) {
        ::IncrementFilenameIndex(sim[0].solver.filename_index,sim[0].solver.index_length);
      }
      if (!m_rank) printf("#%d snapshot projection error, %.15f %.15f %.15f \n",
                          j,
                          sim[0].solver.rom_diff_norms[0],
                          sim[0].solver.rom_diff_norms[1],
                          sim[0].solver.rom_diff_norms[2]);
  }
  ::ResetFilenameIndex( sim[0].solver.filename_index,
                        sim[0].solver.index_length );
  delete recon_init;
}


/*! Check projection error in hyperbolic term*/
void LSROMObject::CheckHyProjError(void* a_s, int idx)
{
  if (!m_rank) {
    std::cout << "----------------------------------------\n";
    std::cout << "Time window # " << idx << " Checking projection error in F(snapshots): ";
    std::cout << std::endl;
  }

  SimulationObject* sim = (SimulationObject*) a_s;
  HyPar  *solver = (HyPar*) &(sim[0].solver);
  Vlasov *param  = (Vlasov*) solver->physics;
  MPIVariables *mpi = (MPIVariables *) param->m;

  std::vector<double> vec_wghosts(sim[0].solver.npoints_local_wghosts*sim[0].solver.nvars);
  std::vector<double> rhs_wghosts(sim[0].solver.npoints_local_wghosts*sim[0].solver.nvars);
  std::vector<double> snap_wghosts(sim[0].solver.npoints_local_wghosts*sim[0].solver.nvars);
  std::vector<double> snapproj_wghosts(sim[0].solver.npoints_local_wghosts*sim[0].solver.nvars);
  std::vector<int> index(sim[0].solver.ndims);

  int num_rows = m_generator[idx]->getSpatialBasis()->numRows();
  int num_cols = m_snapshots[idx]->numColumns();
  CAROM::Vector* recon_init;
  CAROM::Vector* phi_colwork;
  recon_init = new CAROM::Vector(num_rows,false);
  phi_colwork = new CAROM::Vector(num_rows,false);
  CAROM::Vector* m_working;
  m_working = new CAROM::Vector(m_rdims[idx], false);

  for (int j = 0; j < num_cols; j++){

    /* Extend snapshot with ghost points */
    ArrayCopynD(sim[0].solver.ndims,
                m_snapshots[idx]->getColumn(j)->getData(),
                vec_wghosts.data(),
                sim[0].solver.dim_local,
                0,
                sim[0].solver.ghosts,
                index.data(),
                sim[0].solver.nvars);
    sim[0].solver.PostStage( vec_wghosts.data(),
                             &(sim[0].solver),
                             &(sim[0].mpi),
                             0); CHECKERR(ierr);

    /* Evaluate F(\phi_j) */
    TimeRHSFunctionExplicit(rhs_wghosts.data(),
                            vec_wghosts.data(),
                            &(sim[0].solver),
                            &(sim[0].mpi),
                            0);

    /* Remove ghosts point in F(phi_j) */
    ArrayCopynD(sim[0].solver.ndims,
                rhs_wghosts.data(),
                phi_colwork->getData(),
                sim[0].solver.dim_local,
                sim[0].solver.ghosts,
                0,
                index.data(),
                sim[0].solver.nvars);

    /* Check 1 */
    /* Project snapshot onto reduced space and reconstruct */
    m_working = m_basis[idx]->getFirstNColumns(m_rdims[idx])->transposeMult(phi_colwork);
    MPISum_double(m_projected_init[idx]->getData(),m_working->getData(),m_rdims[idx],&mpi->world);

    recon_init = m_basis[idx]->getFirstNColumns(m_rdims[idx])->mult(m_projected_init[idx]);
    ArrayCopynD(sim[0].solver.ndims,
                recon_init->getData(),
                snapproj_wghosts.data(),
                sim[0].solver.dim_local,
                0,
                sim[0].solver.ghosts,
                index.data(),
                sim[0].solver.nvars);

    char buffer[] = "hyprojerr";  // Creates a modifiable buffer and copies the string literal
//  CalSnapROMDiff(&(sim[0].solver),&(sim[0].mpi),rhs_wghosts.data(),snapproj_wghosts.data(),buffer);
    if (!m_rank) printf("F(#%d snapshot) projection error (simply project), %.15f %.15f %.15f \n",
                        j,
                        sim[0].solver.rom_diff_norms[0],
                        sim[0].solver.rom_diff_norms[1],
                        sim[0].solver.rom_diff_norms[2]);

    /* Check 3 */
    m_working = ProjectToRB(m_snapshots[idx]->getColumn(j), m_basis[idx], m_rdims[idx]);
    MPISum_double(m_projected_init[idx]->getData(),m_working->getData(), m_rdims[idx], &mpi->world);
    m_working=m_romhyperb[idx]->mult(m_projected_init[idx]);
    recon_init = m_basis[idx]->getFirstNColumns(m_rdims[idx])->mult(m_working);
    ArrayCopynD(sim[0].solver.ndims,
                recon_init->getData(),
                snapproj_wghosts.data(),
                sim[0].solver.dim_local,
                0,
                sim[0].solver.ghosts,
                index.data(),
                sim[0].solver.nvars);

//  CalSnapROMDiff(&(sim[0].solver),&(sim[0].mpi),rhs_wghosts.data(),snapproj_wghosts.data(),buffer);
    if (!m_rank) printf("F(#%d snapshot) projection error (reduced), %.15f %.15f %.15f \n",
                        j,
                        sim[0].solver.rom_diff_norms[0],
                        sim[0].solver.rom_diff_norms[1],
                        sim[0].solver.rom_diff_norms[2]);

    /* increment the index string, if required */
    if ((!strcmp(sim[0].solver.output_mode,"serial")) && (!strcmp(sim[0].solver.op_overwrite,"no"))) {
      ::IncrementFilenameIndex(sim[0].solver.filename_index,sim[0].solver.index_length);
    }
  }
  ::ResetFilenameIndex( sim[0].solver.filename_index,
                        sim[0].solver.index_length );
  delete recon_init;
  delete phi_colwork;
  delete m_working;
}

/*! Construct potential ROM rhs */
void LSROMObject::ConstructPotentialROMRhs(void* a_s,
                                           const CAROM::Matrix* a_rombasis_f,
                                           const CAROM::Matrix* a_rombasis_phi,
                                           int idx)
{
  SimulationObject* sim = (SimulationObject*) a_s;
  HyPar  *solver = (HyPar*) &(sim[0].solver);
  Vlasov *param  = (Vlasov*) solver->physics;
  MPIVariables *mpi = (MPIVariables *) param->m;

  if (param->ndims_x > 1) {
    fprintf(stderr,"Error in ConstructPotentialROMRhs:\n");
    fprintf(stderr,"  Implemented for 1 spatial dimension only.\n");
  }

  int num_rows = a_rombasis_phi->numRows();
  double* int_f= (double*) calloc(num_rows, sizeof(double));
  std::vector<double> int_f_wghosts(param->npts_local_x_wghosts);
  std::vector<int> index(sim[0].solver.ndims);

  /* Integrate f reduced basis over velocity */
  CAROM::Matrix* integral_basis_f;
  integral_basis_f = new CAROM::Matrix(num_rows, m_rdims[idx], false);
  CAROM::Matrix* m_working;
  m_working = new CAROM::Matrix(m_rdims_phi[idx], m_rdims[idx], false);

  for (int j = 0; j < m_rdims[idx]; j++){
    EvaluatePotentialRhs(a_s, a_rombasis_f->getColumn(j), int_f);
    ArrayCopynD(1,
                int_f,
                int_f_wghosts.data(),
                sim[0].solver.dim_local,
                0,
                sim[0].solver.ghosts,
                index.data(),
                sim[0].solver.nvars);
    char buffer[] = "int_basis";
//  ::VlasovWriteSpatialField(&(sim[0].solver),&(sim[0].mpi),int_f_wghosts.data(),buffer);
    /* Copy \int basis_f dv back to columns of integral_basis_f matrix */
    for (int i = 0; i < num_rows; i++) {
      (*integral_basis_f)(i, j) = int_f[i];
    }
    /* increment the index string, if required */
    if ((!strcmp(sim[0].solver.output_mode,"serial")) && (!strcmp(sim[0].solver.op_overwrite,"no"))) {
      ::IncrementFilenameIndex(sim[0].solver.filename_index,sim[0].solver.index_length);
    }
  }
  ::ResetFilenameIndex( sim[0].solver.filename_index,
                        sim[0].solver.index_length );

  // construct rhs = basis_phi^T integral_basis_f
  m_working = a_rombasis_phi->getFirstNColumns(m_rdims_phi[idx])->transposeMult(integral_basis_f);
  MPISum_double(m_romrhs_phi[idx]->getData(),m_working->getData(),m_rdims_phi[idx]*m_rdims[idx],&mpi->world);
  if (!m_rank) {
    printf("f %d %d\n",a_rombasis_f->numRows(),a_rombasis_f->numColumns());
    printf("integral_basis_f %d %d\n",integral_basis_f->numRows(),integral_basis_f->numColumns());
    printf("m_romrhs_phi %d %d\n",m_romrhs_phi[idx]->numRows(),m_romrhs_phi[idx]->numColumns());
    std::cout << "Checking Potential ROM Rhs: \n";
    for (int i = 0; i < m_romrhs_phi[idx]->numRows(); i++) {
      for (int j = 0; j < m_romrhs_phi[idx]->numColumns(); j++) {
          std::cout << (m_romrhs_phi[idx]->item(i,j)) << " ";
      }
    }
    std::cout << std::endl;
  }
  free(int_f);
  delete integral_basis_f;
  delete m_working;
  return;
}

/*! Construct potential ROM laplacian */
void LSROMObject::ConstructPotentialROMLaplace(void* a_s, const CAROM::Matrix* a_rombasis_phi, int idx)
{
  SimulationObject* sim = (SimulationObject*) a_s;
  HyPar  *solver = (HyPar*) &(sim[0].solver);
  Vlasov *param  = (Vlasov*) solver->physics;
  MPIVariables *mpi = (MPIVariables *) param->m;

  if (param->ndims_x > 1) {
    fprintf(stderr,"Error in ConstructPotentialLaplace:\n");
    fprintf(stderr,"  Implemented for 1 spatial dimension only.\n");
  }

  int *dim    = solver->dim_local;
  int  N      = solver->dim_global[0];
  int  ghosts = solver->ghosts;
  int  ndims  = solver->ndims;

  int num_rows = a_rombasis_phi->numRows();
  int num_cols = a_rombasis_phi->numColumns();

  std::vector<double> vec_wghosts(param->npts_local_x_wghosts);
  std::vector<double> rhs_wghosts(param->npts_local_x_wghosts);
  std::vector<int> index(sim[0].solver.ndims);

  /* Integrate f reduced basis over velocity */
  CAROM::Matrix* laplace_phi;
  laplace_phi = new CAROM::Matrix(num_rows, m_rdims_phi[idx], false);
  CAROM::Vector laplace_col(num_rows,false);
  CAROM::Matrix* m_working;
  m_working = new CAROM::Matrix(m_rdims_phi[idx], m_rdims_phi[idx], false);

  for (int j = 0; j < m_rdims_phi[idx]; j++){
    ArrayCopynD(1,
                a_rombasis_phi->getColumn(j)->getData(),
                vec_wghosts.data(),
                sim[0].solver.dim_local,
                0,
                sim[0].solver.ghosts,
                index.data(),
                sim[0].solver.nvars);

     MPIExchangeBoundaries1D(&(sim[0].mpi),vec_wghosts.data(),param->npts_local_x,
                                  sim[0].solver.ghosts,0,sim[0].solver.ndims);
    ::SecondDerivativeSecondOrderCentralNoGhosts(rhs_wghosts.data(),
                                                 vec_wghosts.data(),
                                                 0,
                                                 1,
                                                 &(param->npts_local_x),
                                                 sim[0].solver.ghosts,
                                                 1,
                                                 &(sim[0].mpi));

    ArrayCopynD(1,
                rhs_wghosts.data(),
                laplace_col.getData(),
                sim[0].solver.dim_local,
                sim[0].solver.ghosts,
                0,
                index.data(),
                sim[0].solver.nvars);

    /* Copy \int basis_f dv back to columns of integral_basis_f matrix */
    for (int i = 0; i < num_rows; i++) {
      (*laplace_phi)(i, j) = -1.0*laplace_col.getData()[i]*(solver->dxinv[0])*(solver->dxinv[0]);
    }
  }

  // construct rhs = basis_phi^T integral_basis_f
  m_working = a_rombasis_phi->getFirstNColumns(m_rdims_phi[idx])->transposeMult(laplace_phi);
  MPISum_double(m_romlaplace_phi[idx]->getData(),m_working->getData(),m_rdims_phi[idx]*m_rdims_phi[idx],&mpi->world);
  if (!m_rank) {
    printf("m_romlaplace_phi %d %d\n",m_romlaplace_phi[idx]->numRows(),m_romlaplace_phi[idx]->numColumns());
//  std::cout << "Checking Potential ROM Laplace: \n";
//  for (int i = 0; i < m_romlaplace_phi[idx]->numRows(); i++) {
//    for (int j = 0; j < m_romlaplace_phi[idx]->numColumns(); j++) {
//        std::cout << (m_romlaplace_phi[idx]->item(i,j)) << " ";
//    }
//  }
//  std::cout << std::endl;
  }
  m_romlaplace_phi[idx]->inverse();

  delete laplace_phi;
  delete m_working;
  return;
}

/*! Check projection error in potential solution */
void LSROMObject::CheckPotentialProjError(void* a_s, int idx)
{
  if (!m_rank) {
    std::cout << "------------------------------------------------\n";
    std::cout << "Checking projection error in potential snapshots: ";
    std::cout << std::endl;
  }
  SimulationObject* sim = (SimulationObject*) a_s;
  HyPar  *solver = (HyPar*) &(sim[0].solver);
  Vlasov *param  = (Vlasov*) solver->physics;
  MPIVariables *mpi = (MPIVariables *) param->m;

  std::vector<double> snap_wghosts(param->npts_local_x_wghosts);
  std::vector<double> snapproj_wghosts(param->npts_local_x_wghosts);
  std::vector<int> index(sim[0].solver.ndims);

  int num_rows = m_generator_phi[idx]->getSpatialBasis()->numRows();
  int num_cols = m_snapshots_phi[idx]->numColumns();
  CAROM::Vector* recon_caromvec;
  recon_caromvec= new CAROM::Vector(num_rows,false);
  CAROM::Vector* m_working;
  m_working = new CAROM::Vector(m_rdims_phi[idx], false);

  for (int j = 0; j < num_cols; j++){
    /* Extend reduced basis \phi_j with ghost points */
    ArrayCopynD(1,
                m_snapshots_phi[idx]->getColumn(j)->getData(),
                snap_wghosts.data(),
                sim[0].solver.dim_local,
                0,
                sim[0].solver.ghosts,
                index.data(),
                sim[0].solver.nvars);

    m_working = ProjectToRB(m_snapshots_phi[idx]->getColumn(j),m_basis_phi[idx], m_rdims_phi[idx]);
    MPISum_double(m_projected_init_phi[idx]->getData(),m_working->getData(),m_rdims_phi[idx],&mpi->world);
    if (!m_rank) {
      printf("%d m_project phi size %d\n",idx,m_projected_init_phi[idx]->dim());
      std::cout << "Checking phi projected coefficients: ";
      for (int k = 0; k < m_projected_init_phi[idx]->dim(); k++) {
            std::cout << (m_projected_init_phi[idx]->item(k)) << " ";
      }
      std::cout << std::endl;
    }
    recon_caromvec= m_basis_phi[idx]->getFirstNColumns(m_rdims_phi[idx])->mult(m_projected_init_phi[idx]);
    ArrayCopynD(1,
                recon_caromvec->getData(),
                snapproj_wghosts.data(),
                sim[0].solver.dim_local,
                0,
                sim[0].solver.ghosts,
                index.data(),
                sim[0].solver.nvars);

    char buffer[] = "snapphiprojerr";  // Creates a modifiable buffer and copies the string literal
    CalSnapROMDiff_phi(&(sim[0].solver),&(sim[0].mpi),snap_wghosts.data(),snapproj_wghosts.data(),buffer);
    /* increment the index string, if required */
    if ((!strcmp(sim[0].solver.output_mode,"serial")) && (!strcmp(sim[0].solver.op_overwrite,"no"))) {
      ::IncrementFilenameIndex(sim[0].solver.filename_index,sim[0].solver.index_length);
    }
    if (!m_rank) printf("#%d snapshot projection error phi, %.15f %.15f %.15f \n",
                        j,
                        sim[0].solver.rom_diff_norms[0],
                        sim[0].solver.rom_diff_norms[1],
                        sim[0].solver.rom_diff_norms[2]);
  }
  ::ResetFilenameIndex( sim[0].solver.filename_index,
                        sim[0].solver.index_length );
  delete recon_caromvec;
}

void LSROMObject::projectInitialSolution_phi(  CAROM::Vector& a_U /*!< solution vector */ )
{
  if (m_generator_phi.size() == 0) {
    if (!m_rank) {
      printf("ERROR in LSROMObject::projectInitialSolution_phi() - m_ls is a vector of size 0.\n");
    }
    return;
  }
  for (int i = 0; i < m_generator_phi.size(); i++) {

    m_projected_init_phi.push_back(new CAROM::Vector(m_rdims_phi[i], false));
    m_projected_init_phi[i] = m_generator_phi[i]->getSpatialBasis()->getFirstNColumns(m_rdims_phi[i])->transposeMult(a_U);

    if (!m_rank) {
      std::cout << "Checking projected coefficients phi: ";
      for (int j = 0; j < m_projected_init_phi[i]->dim(); j++) {
            std::cout << (m_projected_init_phi[i]->item(j)) << " ";
      }
      std::cout << std::endl;
    }
  }
}

int LSROMObject::CalSnapROMDiff_phi( void *s, /*!< Solver object of type #HyPar */
                                     void *m,  /*!< MPI object of type #MPIVariables */
                                     double* a_snapshot,
                                     double* a_romsol,
                                     char* filename)
{
  HyPar* solver = (HyPar*) s;
  MPIVariables* mpi = (MPIVariables*) m;
  Vlasov *param  = (Vlasov*) solver->physics;
  double sum = 0, global_sum = 0;

  static const double tolerance = 1e-15;

  int size = param->npts_local_x_wghosts* solver->nvars;
  double* u_diff = (double*) calloc (size, sizeof(double));

  int* dim = (int*)malloc(sizeof(int));
  *dim = param->npts_local_x;

  /* calculate solution norms (for relative error) */
  double solution_norm[3] = {0.0,0.0,0.0};
  /* L1 */
  sum = ArraySumAbsnD ( solver->nvars,
                        1,
                        dim,
                        solver->ghosts,
                        solver->index,
                        a_snapshot );
  global_sum = 0; MPISum_double(&global_sum,&sum,1,&mpi->world);
  solution_norm[0] = global_sum/((double)param->npts_local_x);
  /* L2 */
  sum = ArraySumSquarenD  ( solver->nvars,
                            1,
                            dim,
                            solver->ghosts,
                            solver->index,
                            a_snapshot );
  global_sum = 0; MPISum_double(&global_sum,&sum,1,&mpi->world);
  solution_norm[1] = sqrt(global_sum/((double)param->npts_local_x));
  /* Linf */
  sum = ArrayMaxnD  ( solver->nvars,
                      1,
                      dim,
                      solver->ghosts,
                      solver->index,
                      a_snapshot  );
  global_sum = 0; MPIMax_double(&global_sum,&sum,1,&mpi->world);
  solution_norm[2] = global_sum;

  /* set u_diff to PDE solution */
  _ArrayCopy1D_(a_snapshot, u_diff, size);
  /* subtract the ROM solutions */
  _ArrayAXPY_(a_romsol, -1.0, u_diff, size);

  /* calculate L1 norm of error */
  sum = ArraySumAbsnD ( solver->nvars,
                        1,
                        dim,
                        solver->ghosts,
                        solver->index,
                        u_diff  );
  global_sum = 0; MPISum_double(&global_sum,&sum,1,&mpi->world);
  solver->rom_diff_norms[0] = global_sum/((double)param->npts_local_x);

  /* calculate L2 norm of error */
  sum = ArraySumSquarenD  ( solver->nvars,
                            1,
                            dim,
                            solver->ghosts,
                            solver->index,
                            u_diff  );
  global_sum = 0; MPISum_double(&global_sum,&sum,1,&mpi->world);
  solver->rom_diff_norms[1] = sqrt(global_sum/((double)param->npts_local_x));

  /* calculate Linf norm of error */
  sum = ArrayMaxnD  ( solver->nvars,
                      1,
                      dim,
                      solver->ghosts,
                      solver->index,
                      u_diff  );
  global_sum = 0; MPIMax_double(&global_sum,&sum,1,&mpi->world);
  solver->rom_diff_norms[2] = global_sum;

  /* decide whether to normalize and report relative diff norms,
    or report absolute diff norms. */
  if (    (solution_norm[0] > tolerance)
      &&  (solution_norm[1] > tolerance)
      &&  (solution_norm[2] > tolerance) ) {
    solver->rom_diff_norms[0] /= solution_norm[0];
    solver->rom_diff_norms[1] /= solution_norm[1];
    solver->rom_diff_norms[2] /= solution_norm[2];
  }

  free(u_diff);
  return 0;
}

/*! Check projection error in laplacian */
void LSROMObject::CheckLaplaceProjError(void* a_s, int idx)
{
  if (!m_rank) {
    std::cout << "------------------------------------------------\n";
    std::cout << "Checking projection error in laplacian of phi: ";
    std::cout << std::endl;
  }
  SimulationObject* sim = (SimulationObject*) a_s;
  HyPar  *solver = (HyPar*) &(sim[0].solver);
  Vlasov *param  = (Vlasov*) solver->physics;
  MPIVariables *mpi = (MPIVariables *) param->m;

  std::vector<double> snap_wghosts(param->npts_local_x_wghosts);
  std::vector<double> lap_snap_wghosts(param->npts_local_x_wghosts);
  std::vector<double> snapproj_wghosts(param->npts_local_x_wghosts);
  std::vector<double> lap_snapproj_wghosts(param->npts_local_x_wghosts);
  std::vector<double> recon_wghosts(param->npts_local_x_wghosts);
  std::vector<int> index(sim[0].solver.ndims);

  int num_rows = m_generator_phi[idx]->getSpatialBasis()->numRows();
  CAROM::Vector* recon_field;
  CAROM::Vector* recon_lap;
  CAROM::Vector* field_reduced;
  CAROM::Vector* lap_reduced;
  CAROM::Vector* m_working;
  m_working = new CAROM::Vector(m_rdims_phi[idx], false);

  int num_cols = m_snapshots_phi[idx]->numColumns();
  char buffer[] = "laplace_projerr";  // Creates a modifiable buffer and copies the string literal

  for (int j = 0; j < num_cols; j++){
    /* extend snapshot with ghost points */
    ArrayCopynD(1,
                m_snapshots_phi[idx]->getColumn(j)->getData(),
                snap_wghosts.data(),
                sim[0].solver.dim_local,
                0,
                sim[0].solver.ghosts,
                index.data(),
                sim[0].solver.nvars);

     MPIExchangeBoundaries1D(&(sim[0].mpi),snap_wghosts.data(),param->npts_local_x,
                                  sim[0].solver.ghosts,0,sim[0].solver.ndims);
    ::SecondDerivativeSecondOrderCentralNoGhosts(lap_snap_wghosts.data(),
                                                 snap_wghosts.data(),
                                                 0,
                                                 1,
                                                 &(param->npts_local_x),
                                                 sim[0].solver.ghosts,
                                                 1,
                                                 &(sim[0].mpi));
    _ArrayScale1D_(lap_snap_wghosts,-1.0*(solver->dxinv[0])*(solver->dxinv[0]),param->npts_local_x_wghosts);

    /* Check 1 */
    /* project snapshot onto reduced space */
    m_working = ProjectToRB(m_snapshots_phi[idx]->getColumn(j),m_basis_phi[idx], m_rdims_phi[idx]);
    MPISum_double(m_projected_init_phi[idx]->getData(),m_working->getData(),m_rdims_phi[idx],&mpi->world);
    lap_reduced = m_romlaplace_phi[idx]->mult(m_projected_init_phi[idx]);
    recon_lap = m_generator_phi[idx]->getSpatialBasis()->getFirstNColumns(m_rdims_phi[idx])->mult(lap_reduced);
    ArrayCopynD(1,
                recon_lap->getData(),
                recon_wghosts.data(),
                sim[0].solver.dim_local,
                0,
                sim[0].solver.ghosts,
                index.data(),
                sim[0].solver.nvars);
    CalSnapROMDiff_phi(&(sim[0].solver),&(sim[0].mpi),lap_snap_wghosts.data(),recon_wghosts.data(),buffer);
    if (!m_rank) printf("#%d snapshot projection error in laplacian (reduced), %.15f %.15f %.15f \n",
                        j,
                        sim[0].solver.rom_diff_norms[0],
                        sim[0].solver.rom_diff_norms[1],
                        sim[0].solver.rom_diff_norms[2]);

    /* Check 2 */
    recon_field = m_generator_phi[idx]->getSpatialBasis()->getFirstNColumns(m_rdims_phi[idx])->mult(m_projected_init_phi[idx]);
    ArrayCopynD(1,
                recon_field->getData(),
                snapproj_wghosts.data(),
                sim[0].solver.dim_local,
                0,
                sim[0].solver.ghosts,
                index.data(),
                sim[0].solver.nvars);

     MPIExchangeBoundaries1D(&(sim[0].mpi),snapproj_wghosts.data(),param->npts_local_x,
                                  sim[0].solver.ghosts,0,sim[0].solver.ndims);
    ::SecondDerivativeSecondOrderCentralNoGhosts(lap_snapproj_wghosts.data(),
                                                 snapproj_wghosts.data(),
                                                 0,
                                                 1,
                                                 &(param->npts_local_x),
                                                 sim[0].solver.ghosts,
                                                 1,
                                                 &(sim[0].mpi));
    _ArrayScale1D_(lap_snapproj_wghosts,-1.0*(solver->dxinv[0])*(solver->dxinv[0]),param->npts_local_x_wghosts);
    ArrayCopynD(1,
                lap_snapproj_wghosts.data(),
                recon_field->getData(),
                sim[0].solver.dim_local,
                sim[0].solver.ghosts,
                0,
                index.data(),
                sim[0].solver.nvars);

    field_reduced = m_generator_phi[idx]->getSpatialBasis()->getFirstNColumns(m_rdims_phi[idx])->transposeMult(recon_field);
    recon_lap = m_generator_phi[idx]->getSpatialBasis()->getFirstNColumns(m_rdims_phi[idx])->mult(field_reduced);
    ArrayCopynD(1,
                recon_lap->getData(),
                recon_wghosts.data(),
                sim[0].solver.dim_local,
                0,
                sim[0].solver.ghosts,
                index.data(),
                sim[0].solver.nvars);

    CalSnapROMDiff_phi(&(sim[0].solver),&(sim[0].mpi),lap_snap_wghosts.data(),recon_wghosts.data(),buffer);
    if (!m_rank) printf("#%d snapshot projection error in laplacian (non reduced), %.15f %.15f %.15f \n",
                        j,
                        sim[0].solver.rom_diff_norms[0],
                        sim[0].solver.rom_diff_norms[1],
                        sim[0].solver.rom_diff_norms[2]);
    /* Check 3 */
    ArrayCopynD(1,
                lap_snap_wghosts.data(),
                recon_lap->getData(),
                sim[0].solver.dim_local,
                sim[0].solver.ghosts,
                0,
                index.data(),
                sim[0].solver.nvars);

    lap_reduced = m_generator_phi[idx]->getSpatialBasis()->getFirstNColumns(m_rdims_phi[idx])->transposeMult(recon_lap);
    recon_lap = m_generator_phi[idx]->getSpatialBasis()->getFirstNColumns(m_rdims_phi[idx])->mult(lap_reduced);
    ArrayCopynD(1,
                recon_lap->getData(),
                recon_wghosts.data(),
                sim[0].solver.dim_local,
                0,
                sim[0].solver.ghosts,
                index.data(),
                sim[0].solver.nvars);
    CalSnapROMDiff_phi(&(sim[0].solver),&(sim[0].mpi),lap_snap_wghosts.data(),recon_wghosts.data(),buffer);
    if (!m_rank) printf("#%d snapshot projection error in laplacian (simply project), %.15f %.15f %.15f \n",
                        j,
                        sim[0].solver.rom_diff_norms[0],
                        sim[0].solver.rom_diff_norms[1],
                        sim[0].solver.rom_diff_norms[2]);

    /* increment the index string, if required */
    if ((!strcmp(sim[0].solver.output_mode,"serial")) && (!strcmp(sim[0].solver.op_overwrite,"no"))) {
      ::IncrementFilenameIndex(sim[0].solver.filename_index,sim[0].solver.index_length);
    }
  }
  ::ResetFilenameIndex( sim[0].solver.filename_index,
                        sim[0].solver.index_length );
}

void LSROMObject::EvaluatePotentialRhs(void* a_s, CAROM::Vector* a_vec, double* m_buffer)
{
  SimulationObject* sim = (SimulationObject*) a_s;
  HyPar  *solver = (HyPar*) &(sim[0].solver);
  Vlasov *param  = (Vlasov*) solver->physics;
  MPIVariables *mpi = (MPIVariables *) param->m;

  if (param->ndims_x > 1) {
    fprintf(stderr,"Error in EvaluatePotentialRhs:\n");
    fprintf(stderr,"  Implemented for 1 spatial dimension only.\n");
  }

  int *dim    = solver->dim_local;
  int  N      = solver->dim_global[0];
  int  ghosts = solver->ghosts;
  int  ndims  = solver->ndims;

  int index[ndims], bounds[ndims], bounds_noghost[ndims], offset[ndims];

  std::vector<double> vec_wghosts(sim[0].solver.npoints_local_wghosts*sim[0].solver.nvars);
  std::vector<int> idx(sim[0].solver.ndims);

  // set bounds for array index to include ghost points
  _ArrayCopy1D_(dim,bounds,ndims);
  for (int k = 0; k < ndims; k++) bounds[k] += 2*ghosts;

  // set bounds for array index to NOT include ghost points
  _ArrayCopy1D_(dim,bounds_noghost,ndims);

  // set offset such that index is compatible with ghost point arrangement
  _ArraySetValue_(offset,ndims,-ghosts);

  ArrayCopynD(sim[0].solver.ndims,
              a_vec->getData(),
              vec_wghosts.data(),
              sim[0].solver.dim_local,
              0,
              sim[0].solver.ghosts,
              idx.data(),
              sim[0].solver.nvars);

  double* basis = vec_wghosts.data();
  // First, integrate the particle distribution over velocity.
  // Since the array dimension we want is not unit stride,
  // first manually add up the local part of the array.
  int done = 0; _ArraySetValue_(index,ndims,0);
  _ArraySetValue_(m_buffer,dim[0],0);
  while (!done) {
    //int p; _ArrayIndex1DWO_(ndims,dim,index,offset,ghosts,p);
    int p; _ArrayIndex1D_(ndims,dim,index,ghosts,p);

    // accumulate f at this spatial location
    double dvinv; _GetCoordinate_(1,index[1],dim,ghosts,solver->dxinv,dvinv);
    double x; _GetCoordinate_(0,index[0],dim,ghosts,solver->x,x);
    double v; _GetCoordinate_(1,index[1],dim,ghosts,solver->x,v);

    m_buffer[index[0]] += basis[p] / dvinv;

    _ArrayIncrementIndex_(ndims,bounds_noghost,index,done);
  }
  // Now we can add up globally using MPI reduction 
  for (int i = 0; i < dim[0]; i++) {
    MPISum_double(&m_buffer[i], &m_buffer[i], 1, &mpi->comm[1]);
  }

  // Find the average density over all x
  double average_velocity = 0.0;
  for (int i = 0; i < dim[0]; i++) {
    average_velocity += m_buffer[i];
  }
  MPISum_double(&average_velocity, &average_velocity, 1, &mpi->comm[0]);
  average_velocity /= (double) N;

  /* Copy \int basis_f dv back to columns of integral_basis_f matrix */
  for (int i = 0; i < dim[0]; i++) {
    m_buffer[i] = m_buffer[i] - average_velocity;
  }
  return; 
}

/*! Check projection error in rhs */
void LSROMObject::CheckRhsProjError(void* a_s, int idx)
{
  if (!m_rank) {
    std::cout << "------------------------------------------------\n";
    std::cout << "Checking projection error in rhs : ";
    std::cout << std::endl;
  }
  SimulationObject* sim = (SimulationObject*) a_s;
  HyPar  *solver = (HyPar*) &(sim[0].solver);
  Vlasov *param  = (Vlasov*) solver->physics;
  MPIVariables *mpi = (MPIVariables *) param->m;

  std::vector<double> int_f_wghosts(param->npts_local_x_wghosts);
  std::vector<double> int_basis_wghosts(param->npts_local_x_wghosts);
  std::vector<double> recon_wghosts(param->npts_local_x_wghosts);
  std::vector<int> index(sim[0].solver.ndims);

  int num_rows = m_generator_phi[idx]->getSpatialBasis()->numRows();
  CAROM::Vector* recon_rhs;
  recon_rhs = new CAROM::Vector(num_rows, false);
  CAROM::Vector* int_f_caromvec;
  int_f_caromvec = new CAROM::Vector(num_rows, false);
  CAROM::Vector* m_w1;
  CAROM::Vector* m_w2;
  CAROM::Vector* rhs_reduced;
  m_w1 = new CAROM::Vector(m_rdims_phi[idx], false);
  m_w2 = new CAROM::Vector(m_rdims[idx], false);
  rhs_reduced = new CAROM::Vector(m_rdims_phi[idx], false);

  double* int_f= (double*) calloc(num_rows, sizeof(double));
  int num_cols = m_snapshots[idx]->numColumns();

  for (int j = 0; j < num_cols; j++){

    EvaluatePotentialRhs(a_s, m_snapshots[idx]->getColumn(j), int_f);

    if (mpi->ip[1] == 0) {
      ArrayCopynD(1,int_f,int_f_wghosts.data(),
                  sim[0].solver.dim_local,
                  0,
                  sim[0].solver.ghosts,
                  index.data(),
                  sim[0].solver.nvars);
      }
      else {
        int_f_wghosts = std::vector<double> (int_f_wghosts.size(),0.0);
      }
//  char buffer1[] = "int_f";
//  ::VlasovWriteSpatialField(&(sim[0].solver),&(sim[0].mpi),int_f_wghosts.data(),buffer1);

    ArrayCopynD(1,
                int_f,
                int_f_caromvec->getData(),
                sim[0].solver.dim_local,
                0,
                0,
                index.data(),
                sim[0].solver.nvars);

    /* Check 1 */
    m_w1 = m_basis_phi[idx]->getFirstNColumns(m_rdims_phi[idx])->transposeMult(int_f_caromvec);
    MPISum_double(rhs_reduced->getData(),m_w1->getData(),m_rdims_phi[idx],&mpi->world);
    if (!m_rank) {
      std::cout << "Checking rhs reduced: ";
      for (int k = 0; k < rhs_reduced->dim(); k++) {
            std::cout << (rhs_reduced->item(k)) << " ";
      }
      std::cout << std::endl;
    }
    recon_rhs = m_basis_phi[idx]->getFirstNColumns(m_rdims_phi[idx])->mult(rhs_reduced);
    ArrayCopynD(1,
                recon_rhs->getData(),
                recon_wghosts.data(),
                sim[0].solver.dim_local,
                0,
                sim[0].solver.ghosts,
                index.data(),
                sim[0].solver.nvars);
    char buffer[] = "int_f_projerr";  // Creates a modifiable buffer and copies the string literal
    CalSnapROMDiff_phi(&(sim[0].solver),&(sim[0].mpi),int_f_wghosts.data(),recon_wghosts.data(),buffer);
    if (!m_rank) printf("#%d snapshot projection error in rhs (non reduced), %.15f %.15f %.15f \n",
                        j,
                        sim[0].solver.rom_diff_norms[0],
                        sim[0].solver.rom_diff_norms[1],
                        sim[0].solver.rom_diff_norms[2]);

    /* Check 2 */
//  projectInitialSolution(*(m_snapshots[0]->getColumn(j)));
    m_w2 = ProjectToRB(m_snapshots[idx]->getColumn(j),m_basis[idx], m_rdims[idx]);
    MPISum_double(m_projected_init[idx]->getData(),m_w2->getData(),m_rdims[idx],&mpi->world);
    rhs_reduced = m_romrhs_phi[idx]->mult(m_projected_init[idx]);
    recon_rhs = m_generator_phi[idx]->getSpatialBasis()->getFirstNColumns(m_rdims_phi[idx])->mult(rhs_reduced);
    ArrayCopynD(1,
                recon_rhs->getData(),
                recon_wghosts.data(),
                sim[0].solver.dim_local,
                0,
                sim[0].solver.ghosts,
                index.data(),
                sim[0].solver.nvars);
    CalSnapROMDiff_phi(&(sim[0].solver),&(sim[0].mpi),int_f_wghosts.data(),recon_wghosts.data(),buffer);
    if (!m_rank) printf("#%d snapshot projection error in rhs (reduced), %.15f %.15f %.15f \n",
                        j,
                        sim[0].solver.rom_diff_norms[0],
                        sim[0].solver.rom_diff_norms[1],
                        sim[0].solver.rom_diff_norms[2]);


    /* increment the index string, if required */
    if ((!strcmp(sim[0].solver.output_mode,"serial")) && (!strcmp(sim[0].solver.op_overwrite,"no"))) {
      ::IncrementFilenameIndex(sim[0].solver.filename_index,sim[0].solver.index_length);
    }
  }
  ::ResetFilenameIndex( sim[0].solver.filename_index,
                        sim[0].solver.index_length );
  free(int_f);
  delete recon_rhs;
  delete int_f_caromvec;
  delete rhs_reduced;
}

/*! Check projection error in E solution */
void LSROMObject::CheckEProjError(void* a_s, int idx)
{
  if (!m_rank) {
    std::cout << "------------------------------------------------\n";
    std::cout << "Checking projection error in E snapshots: ";
    std::cout << std::endl;
  }
  SimulationObject* sim = (SimulationObject*) a_s;
  HyPar  *solver = (HyPar*) &(sim[0].solver);
  Vlasov *param  = (Vlasov*) solver->physics;
  MPIVariables *mpi = (MPIVariables *) param->m;

  std::vector<double> snap_wghosts(param->npts_local_x_wghosts);
  std::vector<double> snapproj_wghosts(param->npts_local_x_wghosts);
  std::vector<double> vec_x_wghosts(param->npts_local_x_wghosts*param->ndims_x);
  std::vector<int> index(sim[0].solver.ndims);

  int num_rows = m_generator_phi[idx]->getSpatialBasis()->numRows();
  int num_cols = m_snapshots_e[idx]->numColumns();
  CAROM::Vector* recon_caromvec;
  CAROM::Vector* m_working;
  m_working = new CAROM::Vector(m_rdims_phi[idx], false);

  for (int j = 0; j < num_cols; j++){
    ArrayCopynD(1,
                m_snapshots_e[idx]->getColumn(j)->getData(),
                snap_wghosts.data(),
                sim[0].solver.dim_local,
                0,
                sim[0].solver.ghosts,
                index.data(),
                sim[0].solver.nvars);

//  projectInitialSolution_phi(*(m_snapshots_phi[idx]->getColumn(j)));
    m_working = ProjectToRB(m_snapshots_phi[idx]->getColumn(j),m_basis_phi[idx], m_rdims_phi[idx]);
    MPISum_double(m_projected_init_phi[idx]->getData(),m_working->getData(),m_rdims_phi[idx],&mpi->world);
    recon_caromvec = m_basis_e[idx]->mult(m_projected_init_phi[idx]);
    ArrayCopynD(1,
                recon_caromvec->getData(),
                vec_x_wghosts.data(),
                sim[0].solver.dim_local,
                0,
                sim[0].solver.ghosts,
                index.data(),
                sim[0].solver.nvars);

    char buffer[] = "snapeprojerr";  // Creates a modifiable buffer and copies the string literal
    CalSnapROMDiff_phi(&(sim[0].solver),&(sim[0].mpi),snap_wghosts.data(),vec_x_wghosts.data(),buffer);

    /* increment the index string, if required */
    if ((!strcmp(sim[0].solver.output_mode,"serial")) && (!strcmp(sim[0].solver.op_overwrite,"no"))) {
      ::IncrementFilenameIndex(sim[0].solver.filename_index,sim[0].solver.index_length);
    }
    if (!m_rank) printf("#%d snapshot projection error e, %.15f %.15f %.15f \n",
                        j,
                        sim[0].solver.rom_diff_norms[0],
                        sim[0].solver.rom_diff_norms[1],
                        sim[0].solver.rom_diff_norms[2]);

  }
  ::ResetFilenameIndex( sim[0].solver.filename_index,
                        sim[0].solver.index_length );
  delete recon_caromvec;
}

/*! Dump ROM basis phi */
void LSROMObject::OutputROMBasisPhi(void* a_s, const CAROM::Matrix* a_rombasis, int idx)
{
  SimulationObject* sim = (SimulationObject*) a_s;
  HyPar  *solver = (HyPar*) &(sim[0].solver);
  Vlasov *param  = (Vlasov*) solver->physics;
  MPIVariables *mpi = (MPIVariables *) param->m;

  std::vector<double> vec_wghosts(param->npts_local_x_wghosts);
  std::vector<int> index(sim[0].solver.ndims);

  for (int j = 0; j < a_rombasis->numColumns(); j++){

    ArrayCopynD(1,
                a_rombasis->getColumn(j)->getData(),
                vec_wghosts.data(),
                sim[0].solver.dim_local,
                0,
                sim[0].solver.ghosts,
                index.data(),
                sim[0].solver.nvars);
    char buffer[] = "basis_phi";
    char fbuffer[100];
    sprintf(fbuffer, "%s_%d",buffer,idx);
    ::VlasovWriteSpatialField(&(sim[0].solver),&(sim[0].mpi),vec_wghosts.data(),fbuffer);

    /* increment the index string, if required */
    if ((!strcmp(sim[0].solver.output_mode,"serial")) && (!strcmp(sim[0].solver.op_overwrite,"no"))) {
        ::IncrementFilenameIndex(sim[0].solver.filename_index,sim[0].solver.index_length);
    }
  }

  ::ResetFilenameIndex( sim[0].solver.filename_index,
                        sim[0].solver.index_length );
  return;
}

/*! Construct E basis from phi basis */
void LSROMObject::ConstructEBasis(void* a_s,int idx)
{
  if (!m_rank) {
    std::cout << "------------------------------------------------\n";
    std::cout << "Construct E Basis: ";
    std::cout << std::endl;
  }
  SimulationObject* sim = (SimulationObject*) a_s;
  HyPar  *solver = (HyPar*) &(sim[0].solver);
  Vlasov *param  = (Vlasov*) solver->physics;
  MPIVariables *mpi = (MPIVariables *) param->m;

  std::vector<double> basis_vec_wghosts(param->npts_local_x_wghosts);
  std::vector<double> vec_x_wghosts(param->npts_local_x_wghosts*param->ndims_x);
  std::vector<int> index(sim[0].solver.ndims);

  int num_rows = m_basis_phi[idx]->numRows();
  double* vec_noghosts = (double*) calloc(num_rows, sizeof(double));

  for (int j = 0; j < m_rdims_phi[idx]; j++){

    ArrayCopynD(1,
                m_basis_phi[idx]->getColumn(j)->getData(),
                basis_vec_wghosts.data(),
                sim[0].solver.dim_local,
                0,
                sim[0].solver.ghosts,
                index.data(),
                sim[0].solver.nvars);

     MPIExchangeBoundaries1D(&(sim[0].mpi),basis_vec_wghosts.data(),param->npts_local_x,
                                  sim[0].solver.ghosts,0,sim[0].solver.ndims);
    ::FirstDerivativeSecondOrderCentralNoGhosts(vec_x_wghosts.data(),
                                                basis_vec_wghosts.data(),
                                                0,
                                                1,
                                                1,
                                                &(param->npts_local_x),
                                                sim[0].solver.ghosts,
                                                1,
                                                &(sim[0].mpi));

    if (mpi->ip[1] != 0) vec_x_wghosts = std::vector<double> (vec_x_wghosts.size(),0.0);

    _ArrayScale1D_(vec_x_wghosts,-1.0*(solver->dxinv[0]),param->npts_local_x_wghosts)

    ArrayCopynD(1,
                vec_x_wghosts.data(),
                vec_noghosts,
                sim[0].solver.dim_local,
                sim[0].solver.ghosts,
                0,
                index.data(),
                sim[0].solver.nvars);
    for (int i = 0; i < num_rows; i++) {
      (*m_basis_e[idx])(i, j) = vec_noghosts[i];
    }
  }
  free(vec_noghosts);
}

/*! Construct reduced hyperbolic operator in x direction */
void LSROMObject::ConstructROMHy_x(void* a_s, const CAROM::Matrix* a_rombasis, int idx)
{
  SimulationObject* sim = (SimulationObject*) a_s;
  HyPar  *solver = (HyPar*) &(sim[0].solver);
  Vlasov *param  = (Vlasov*) solver->physics;
  MPIVariables *mpi = (MPIVariables *) param->m;

  int num_rows = a_rombasis->numRows();
  int num_cols = a_rombasis->numColumns();

  /* Compute F(\Phi), where \Phi is the reduced basis matrix */
  CAROM::Matrix* phi_hyper_x;
  phi_hyper_x = new CAROM::Matrix(num_rows, m_rdims[idx], false);
  CAROM::Vector phi_hyper_col(num_rows,false);
  CAROM::Matrix* m_working;
  m_working = new CAROM::Matrix(m_rdims[idx], m_rdims[idx], false);

  std::vector<double> vec_wghosts(sim[0].solver.npoints_local_wghosts*sim[0].solver.nvars);
  std::vector<double> rhs_wghosts(sim[0].solver.npoints_local_wghosts*sim[0].solver.nvars);
  std::vector<int> index(sim[0].solver.ndims);

  for (int j = 0; j < m_rdims[idx]; j++){
    ArrayCopynD(sim[0].solver.ndims,
                a_rombasis->getColumn(j)->getData(),
                vec_wghosts.data(),
                sim[0].solver.dim_local,
                0,
                sim[0].solver.ghosts,
                index.data(),
                sim[0].solver.nvars);

    solver->ApplyBoundaryConditions(solver,mpi,vec_wghosts.data(),NULL,0);
    solver->ApplyIBConditions(solver,mpi,vec_wghosts.data(),0);
    MPIExchangeBoundariesnD(  solver->ndims,
                              solver->nvars,
                              solver->dim_local,
                              solver->ghosts,
                              mpi,
                              vec_wghosts.data());
    HyperbolicFunction_1dir( rhs_wghosts.data(),
                             vec_wghosts.data(),
                             solver,
                             mpi,
                             0,
                             1,
                             solver->FFunction,
                             solver->Upwind, 0 );

    /* Remove ghosts point in F(phi_j) */
    ArrayCopynD(sim[0].solver.ndims,
                rhs_wghosts.data(),
                phi_hyper_col.getData(),
                sim[0].solver.dim_local,
                sim[0].solver.ghosts,
                0,
                index.data(),
                sim[0].solver.nvars);

    /* Copy F(phi_j) back to columns of phi_hyper matrix */
    for (int i = 0; i < num_rows; i++) {
      (*phi_hyper_x)(i, j) = phi_hyper_col.getData()[i];
    }
  }

  // construct hyper_ROM = phi^T phi_hyper
  m_working = a_rombasis->getFirstNColumns(m_rdims[idx])->transposeMult(phi_hyper_x);
  MPISum_double(m_romhyperb_x[idx]->getData(),m_working->getData(),m_rdims[idx]*m_rdims[idx],&mpi->world);
  if (!m_rank) {
    printf("m_romhyperb_x %d %d\n",m_romhyperb_x[idx]->numRows(),m_romhyperb_x[idx]->numColumns());
//  std::cout << "Checking romhyperb_x: \n";
//  for (int i = 0; i < m_romhyperb_x[idx]->numRows(); i++) {
//    for (int j = 0; j < m_romhyperb_x[idx]->numColumns(); j++) {
//        std::cout << (m_romhyperb_x[idx]->item(i,j)) << " ";
//    }
//  }
//  std::cout << std::endl;
  }
  if (!m_rank) printf("m_romhyperb_x %d %d\n",m_romhyperb_x[idx]->numRows(),m_romhyperb_x[idx]->numColumns());

  delete phi_hyper_x;

  return;
}

/*! Construct reduced hyperbolic tensor in v direction */
void LSROMObject::ConstructROMHy_v(void* a_s, const CAROM::Matrix* a_rombasis, const CAROM::Matrix* a_rombasis_phi, int idx)
{
  SimulationObject* sim = (SimulationObject*) a_s;
  HyPar  *solver = (HyPar*) &(sim[0].solver);
  Vlasov *param  = (Vlasov*) solver->physics;
  MPIVariables *mpi = (MPIVariables *) param->m;

  int num_rows = a_rombasis->numRows();
  int num_cols = a_rombasis->numColumns();

  /* Compute F(\Phi), where \Phi is the reduced basis matrix */
  CAROM::Vector phi_hyper_col(num_rows,false);
  CAROM::Matrix* m_working;
  m_working = new CAROM::Matrix(m_rdims_phi[idx], m_rdims[idx], false);

  std::vector<int> index(sim[0].solver.ndims);
  std::vector<double> vec_wghosts(sim[0].solver.npoints_local_wghosts*sim[0].solver.nvars);
  std::vector<double> vec_wo_ghosts(param->npts_local_x);
  std::vector<double> vec_wo_ghosts1(param->npts_local_x);
  std::vector<double> rhs_wghosts(sim[0].solver.npoints_local_wghosts*sim[0].solver.nvars);

  for (int i = 0; i < m_rdims[idx]; i++) {

    m_romhyperb_v[idx].push_back(new CAROM::Matrix(m_rdims_phi[idx], m_rdims[idx], false));
    for (int j = 0; j < m_rdims_phi[idx]; j++) {

      for (int k = 0; k < m_rdims[idx]; k++) {
        ArrayCopynD(sim[0].solver.ndims,
                    a_rombasis->getColumn(k)->getData(),
                    vec_wghosts.data(),
                    sim[0].solver.dim_local,
                    0,
                    sim[0].solver.ghosts,
                    index.data(),
                    sim[0].solver.nvars);

        if (mpi->ip[1] == 0) {
          ArrayCopynD(1,
                      a_rombasis_phi->getColumn(j)->getData(),
                      vec_wo_ghosts.data(),
                      sim[0].solver.dim_local,
                      0,
                      0,
                      index.data(),
                      sim[0].solver.nvars);
        }
        else {
          vec_wo_ghosts = std::vector<double> (vec_wo_ghosts.size(),0.0);
        }
        MPISum_double(vec_wo_ghosts1.data(),vec_wo_ghosts.data(),param->npts_local_x,&mpi->comm[1]);

        memset(param->e_field, 0, sizeof(double) * param->npts_local_x_wghosts);

        ArrayCopynD(1,
                    vec_wo_ghosts1.data(),
                    param->e_field,
                    sim[0].solver.dim_local,
                    0,
                    sim[0].solver.ghosts,
                    index.data(),
                    sim[0].solver.nvars);

        solver->ApplyBoundaryConditions(solver,mpi,vec_wghosts.data(),NULL,0);
        solver->ApplyIBConditions(solver,mpi,vec_wghosts.data(),0);
        MPIExchangeBoundariesnD(  solver->ndims,
                                  solver->nvars,
                                  solver->dim_local,
                                  solver->ghosts,
                                  mpi,
                                  vec_wghosts.data());
        HyperbolicFunction_1dir( rhs_wghosts.data(),
                                 vec_wghosts.data(),
                                 solver,
                                 mpi,
                                 0,
                                 1,
                                 solver->FFunction,
                                 solver->Upwind, 1 );

        ArrayCopynD(sim[0].solver.ndims,
                    rhs_wghosts.data(),
                    phi_hyper_col.getData(),
                    sim[0].solver.dim_local,
                    sim[0].solver.ghosts,
                    0,
                    index.data(),
                    sim[0].solver.nvars);

        (*m_working)(j, k) = a_rombasis->getColumn(i)->inner_product(phi_hyper_col);
      }
    }
    MPISum_double(m_romhyperb_v[idx][i]->getData(),m_working->getData(),m_rdims_phi[idx]*m_rdims[idx],&mpi->world);
//  if (!m_rank) {
//    printf("matrix size %d %d\n",m_romhyperb_v[idx][i]->numRows(),m_romhyperb_v[idx][i]->numColumns());
//  for (int j = 0; j < m_rdims_phi[idx]; j++) {
//    for (int k = 0; k < m_rdims[idx]; k++) {
//      printf("i th %d m_rank %d checking m_romhyperb_v %d %d %f\n",i,m_rank,j,k,((*(m_romhyperb_v[idx][i]))(j, k)));
//    }
//  }
//  }
  }
  return;
}

/*! Write snapshot matrix */
void LSROMObject::writeSnapshot(void* a_s)
{
  if (m_generator.size() > 0) {
    for (int i = 0; i < m_generator.size(); i++) {
      m_generator[i]->writeSnapshot();
      printf("checking snapshot dimension %d %d \n",m_generator[i]->getSnapshotMatrix()->numRows(),m_generator[i]->getSnapshotMatrix()->numColumns());
      delete m_generator[i];
      delete m_options[i];
    }
  }

  if (m_solve_phi) {
    if (m_generator_phi.size() > 0) {
      for (int i = 0; i < m_generator_phi.size(); i++) {
        m_generator_phi[i]->writeSnapshot();
        printf("checking snapshot dimension phi %d %d \n",m_generator_phi[i]->getSnapshotMatrix()->numRows(),m_generator_phi[i]->getSnapshotMatrix()->numColumns());
        delete m_generator_phi[i];
        delete m_options_phi[i];
      }
    }
  }
  outfile_twp.open(outputPath + "/" + std::string(twpfile));
  outfile_twp << m_generator.size();
  if (m_solve_phi) outfile_twp << ", " << m_generator_phi.size();
  outfile_twp.close();
  return;
}

/*! Merge */
void LSROMObject::merge(void* a_s)
{
  SimulationObject* sim = (SimulationObject*) a_s;
  HyPar  *solver = (HyPar*) &(sim[0].solver);
  Vlasov *param  = (Vlasov*) solver->physics;
  MPIVariables *mpi = (MPIVariables *) param->m;

  ReadTimeWindows(a_s);

  for (int sampleWindow = 0; sampleWindow < m_numwindows; ++sampleWindow)
  {
    const std::string basisFileName = basisName + "_" + std::to_string(sampleWindow);
	  std::unique_ptr<CAROM::BasisGenerator> basis_generator;
	  CAROM::Options* options = new CAROM::Options(m_vec_size, max_num_snapshots, 1, update_right_SV);
	  CAROM::BasisGenerator* generator = new CAROM::BasisGenerator(*options, isIncremental, basisFileName);
    for (int paramID=0; paramID<m_nsets; ++paramID)
    {
	    std::string snapshot_filename = basisName + std::to_string(
                                      paramID) + "_" + std::to_string(sampleWindow)
	                                    + "_snapshot";
      generator->loadSamples(snapshot_filename,"snapshot");
    }
    generator->endSamples(); // save the merged basis f
    delete generator;
    delete options;
  }

  if (m_solve_phi) {
	  CAROM::Options* options_phi = new CAROM::Options(param->npts_local_x, max_num_snapshots, 1, update_right_SV);
	  CAROM::BasisGenerator* generator_phi = new CAROM::BasisGenerator(*options_phi, isIncremental, basisName_phi);
    for (int paramID=0; paramID<m_nsets; ++paramID)
    {
      std::string snapshot_filename = basisName_phi + std::to_string(
                                      paramID) + "_" + std::to_string(0)
	                                    + "_snapshot";
      generator_phi->loadSamples(snapshot_filename,"snapshot");
    }
    generator_phi->endSamples(); // save the merged basis f
    delete generator_phi;
    delete options_phi;
  }

	return;
}

/*! Merge */
void LSROMObject::online(void* a_s)
{
  SimulationObject* sim = (SimulationObject*) a_s;
  HyPar  *solver = (HyPar*) &(sim[0].solver);
  Vlasov *param  = (Vlasov*) solver->physics;
  MPIVariables *mpi = (MPIVariables *) param->m;

  ReadTimeWindows(a_s);

  for (int sampleWindow = 0; sampleWindow < m_numwindows; ++sampleWindow)
  {
    const std::string basisFileName = basisName + "_" + std::to_string(sampleWindow);
	  CAROM::BasisReader reader(basisFileName);
    const CAROM::Matrix* spatialbasis;
    if (m_energy_criteria > 0)
    {
      spatialbasis = reader.getSpatialBasis(0.0, 1-m_energy_criteria);
    }
    else
    {
      spatialbasis = reader.getSpatialBasis(0.0, m_rdim);
    }

    const CAROM::Vector* singularf;
    if (!m_rank) {
      singularf = reader.getSingularValues(0);
      std::cout << "----------------\n";
      std::cout << "Time window #" << sampleWindow << ": Singular values of f: ";
      for (int i = 0; i < singularf->dim(); i++) {
            std::cout << (singularf->item(i)) << " ";
      }
      std::cout << "\n";
      std::cout << std::endl;
    }

    m_basis.push_back(new CAROM::Matrix(
                                spatialbasis->getData(),
                                spatialbasis->numRows(),
                                spatialbasis->numColumns(),
                                false,
                                true));

	  int numRowRB, numColumnRB;
    numRowRB = m_basis[sampleWindow]->numRows();
    numColumnRB = m_basis[sampleWindow]->numColumns();
    if (!m_rank) printf("spatial basis dimension is %d x %d\n", numRowRB,
                        numColumnRB);
    if (m_energy_criteria > 0) m_rdim = numColumnRB;
    m_rdims.push_back(m_rdim);
    if (!m_rank) {
      std::cout << "----------------------------------------\n";
      std::cout << "Time window #" << sampleWindow << ": # f POD basis : " << m_rdims[sampleWindow];
      std::cout << std::endl;
    }
    if (sampleWindow > 0) {
      m_fullscale.push_back(new CAROM::Matrix(m_rdims[sampleWindow], m_rdims[sampleWindow-1], false));
      m_matrix = m_basis[sampleWindow]->transposeMult(m_basis[sampleWindow-1]);
      MPISum_double(m_fullscale[sampleWindow-1]->getData(), m_matrix->getData(), m_rdims[sampleWindow-1]*m_rdims[sampleWindow], &mpi->world);
    }
    m_projected_init.push_back(new CAROM::Vector(m_rdims[sampleWindow], false));
    m_romcoef.push_back(new CAROM::Vector(m_rdims[sampleWindow], false));
    m_intervals.push_back( Interval(0, m_t_final) );
    m_snap.push_back(0);

    if ((!m_solve_phi) && (!m_direct_comp_hyperbolic)) {
      m_romhyperb.push_back(new CAROM::Matrix(m_rdims[sampleWindow], m_rdims[sampleWindow], false));
      ConstructROMHy(a_s, m_basis[sampleWindow], sampleWindow);
    }
  }

  if (m_solve_phi) {
	  CAROM::BasisReader reader_phi(basisName_phi);
    const CAROM::Matrix* spatialbasis_phi;
    if (m_energy_criteria > 0)
    {
      spatialbasis_phi = reader_phi.getSpatialBasis(0.0, 1-m_energy_criteria);
    }
    else
    {
      spatialbasis_phi = reader_phi.getSpatialBasis(0.0, m_rdim);
    }

    const CAROM::Vector* singularphi;
    if (!m_rank) {
      singularphi = reader_phi.getSingularValues(0);
      std::cout << "----------------\n";
      std::cout << "Singular Values of potential: ";
      for (int i = 0; i < singularphi->dim(); i++) {
            std::cout << (singularphi->item(i)) << " ";
      }
      std::cout << "\n";
      std::cout << std::endl;
    }

    m_basis_phi.push_back(new CAROM::Matrix(
                                spatialbasis_phi->getData(),
                                spatialbasis_phi->numRows(),
                                spatialbasis_phi->numColumns(),
                                false,
                                true));
	  int numRowRB_phi, numColumnRB_phi;
    numRowRB_phi = m_basis_phi[0]->numRows();
    numColumnRB_phi = m_basis_phi[0]->numColumns();
    if (!m_rank) printf("spatial basis dimension is %d x %d\n", numRowRB_phi,
                        numColumnRB_phi);
    if (m_energy_criteria > 0) m_rdim_phi = numColumnRB_phi;
    m_rdims_phi.push_back(m_rdim_phi);
    if (!m_rank) {
      std::cout << "----------------------------------------\n";
      std::cout << "Time window #" << 0 << ": # Potential POD basis : " << m_rdims_phi[0];
      std::cout << std::endl;
    }
    m_projected_init_phi.push_back(new CAROM::Vector(m_rdims_phi[0], false));
    m_romrhs_phi.push_back(new CAROM::Matrix(m_rdims_phi[0], m_rdims[0], false));
    ConstructPotentialROMRhs(a_s, m_basis[0], m_basis_phi[0], 0);

    m_romlaplace_phi.push_back(new CAROM::Matrix(m_rdims_phi[0], m_rdims_phi[0], false));
    ConstructPotentialROMLaplace(a_s, m_basis_phi[0], 0);
    m_basis_e.push_back(new CAROM::Matrix(m_basis_phi[0]->numRows(),
                        m_rdims_phi[0], false));
    ConstructEBasis(a_s, 0);
    m_romhyperb_x.push_back(new CAROM::Matrix(m_rdims[0], m_rdims[0], false));
    ConstructROMHy_x(a_s, m_basis[0], 0);
    m_romhyperb_v.push_back(std::vector<CAROM::Matrix*>());
    ConstructROMHy_v(a_s, m_basis[0], m_basis_e[0], 0);
  }

}

void LSROMObject::ReadTimeWindows(void* a_s)
{
  if (!m_rank) std::cout << "Reading time window parameters from file " << std::string(twpfile) << std::endl;
  std::ifstream ifs(std::string(twpfile).c_str());
  if (!ifs.is_open())
  {
    std::cout << "Error: invalid file" << std::endl;
  }
  const int nparamRead = m_solve_phi ? 2 : 1;

  std::string line, word;
  int count = 0;
  while (getline(ifs, line))
  {
    std::stringstream s(line);
    std::vector<std::string> row;

    while (getline(s, word, ','))
      row.push_back(word);
    if (row.size() != nparamRead)
    {
      std::cout << "Error: CSV file does not specify " << nparamRead << " parameters" << std::endl;
      ifs.close();
    }
    m_numwindows = stoi(row[0]);
    if (m_solve_phi) m_numwindows_phi = stoi(row[1]);
  }
  ifs.close();
}

#endif

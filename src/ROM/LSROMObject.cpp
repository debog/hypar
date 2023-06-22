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
extern "C" int FirstDerivativeSecondOrderCentralNoGhosts (double*,double*,int,int,void*,void*);
extern "C" int SecondDerivativeSecondOrderCentralNoGhosts (double*,double*,int,void*,void*);
extern "C" int VlasovAdvection_x(double*,double*,int,void*,double);
extern "C" int HyperbolicFunction_1dir (double*,double*,void*,void*,double,int,
                                        int(*)(double*,double*,int,void*,double),
                                        int(*)(double*,double*,double*,double*,double*,
                                        double*,int,void*,double),int);

LSROMObject::LSROMObject(   const int     a_vec_size,       /*!< vector size */
                            const double  a_dt,             /*!< time step size */
                            const int     a_sampling_freq,  /*!< sampling frequency */
                            const int     a_rdim,           /*!< latent space dimension */
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

  /* precomputation idea */
  m_options_phi.clear();
  m_generator_phi.clear();
  m_projected_init_phi.clear();

  m_ls_is_trained.clear();
  m_intervals.clear();
  m_vec_size = a_vec_size;
  m_dt = a_dt;
  m_sampling_freq = a_sampling_freq;
  m_t_final = -1;
  m_rdim = a_rdim;
  m_rdim_phi = m_rdim;

  m_sim_idx = a_sim_idx;
  m_var_idx = a_var_idx;

  m_num_window_samples = INT_MAX;

  char dirname_c_str[_MAX_STRING_SIZE_] = "LS";
  char write_snapshot_mat[_MAX_STRING_SIZE_] = "false";
  char direct_comp_hyperbolic[_MAX_STRING_SIZE_] = "false";
  char solve_phi[_MAX_STRING_SIZE_] = "false";

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
    printf("  potential latent space dimension:  %d\n", m_rdim_phi);
    printf("  directory name for LS objects: %s\n", dirname_c_str);
    printf("  write snapshot matrix to file:  %s\n", write_snapshot_mat);
    printf("  directly compute hyperbolic term:  %s\n", direct_comp_hyperbolic);
    printf("  solve potential:  %s\n", solve_phi);
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
#endif

  m_dirname = std::string( dirname_c_str );
  m_write_snapshot_mat = (std::string(write_snapshot_mat) == "true");
  m_direct_comp_hyperbolic = (std::string(direct_comp_hyperbolic) == "true");
  m_solve_phi = (std::string(solve_phi) == "true");

  if (m_num_window_samples <= m_rdim) {
    printf("ERROR:LSROMObject::LSROMObject() - m_num_window_samples <= m_rdim!!");
  }

  if (m_rdim_phi > m_rdim) {
    printf("ERROR:LSROMObject::LSROMObject() - m_rdim_phi > m_rdim!!");
  }

  m_tic = 0;
  m_curr_win = 0;
  m_snap = 0;
}

void LSROMObject::projectInitialSolution(  CAROM::Vector& a_U /*!< solution vector */ )
{
  /* Need to modify the code so that is works for multiple windows */
  /* Assumming m_generator[i]->getSpatialBasis() is called before */
  if (m_generator.size() == 0) {
    if (!m_rank) {
      printf("ERROR in LSROMObject::projectInitialSolution() - m_ls is a vector of size 0.\n");
    }
    return;
  }
  for (int i = 0; i < m_generator.size(); i++) {

    m_projected_init.push_back(new CAROM::Vector(m_rdim, false));
    m_projected_init[i] = m_generator[i]->getSpatialBasis()->getFirstNColumns(m_rdim)->transposeMult(a_U);

//  if (!m_rank) {
//    printf("m_project size %d\n",m_projected_init[i]->dim());
//    std::cout << "Checking projected coefficients: ";
//    for (int j = 0; j < m_projected_init[i]->dim(); j++) {
//          std::cout << (m_projected_init[i]->item(j)) << " ";
//    }
//    std::cout << std::endl;
//  }

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
  std::vector<double> vec_wghosts(param->npts_local_x*param->ndims_x);
  std::vector<double> vec_x_wghosts(param->npts_local_x_wghosts*param->ndims_x);
  std::vector<int> index(sim[0].solver.ndims);
  if (m_tic == 0) {

    m_options.push_back(new CAROM::Options(m_vec_size, max_num_snapshots, 1, update_right_SV));
    m_generator.push_back(new CAROM::BasisGenerator(*m_options[m_curr_win], isIncremental, basisName));
    m_intervals.push_back( Interval(a_time, m_t_final) );
    m_ls_is_trained.push_back(false);

    if (!m_rank) {
      printf( "LSROMObject::takeSample() - creating new generator object for sim. domain %d, var %d, t=%f (total: %d).\n",
              m_sim_idx, m_var_idx, m_intervals[m_curr_win].first, m_generator.size());
    }

    bool addSample = m_generator[m_curr_win]->takeSample( a_U.getData(), a_time, m_dt );
    char buffer[] = "sample";  // Creates a modifiable buffer and copies the string literal
    OutputlibROMfield(m_generator[0]->getSnapshotMatrix()->getColumn(0)->getData(),
                      sim[0],
                      buffer);

    m_options_phi.push_back(new CAROM::Options(param->npts_local_x, max_num_snapshots, 1, update_right_SV));
    m_generator_phi.push_back(new CAROM::BasisGenerator(*m_options_phi[m_curr_win], isIncremental, basisName));
    ArrayCopynD(1,
                param->potential,
                vec_wghosts.data(),
                sim[0].solver.dim_local,
                sim[0].solver.ghosts,
                0,
                index.data(),
                sim[0].solver.nvars);
    bool addphiSample = m_generator_phi[m_curr_win]->takeSample( vec_wghosts.data(), a_time, m_dt );

    m_options_e.push_back(new CAROM::Options(param->npts_local_x, max_num_snapshots, 1, update_right_SV));
    m_generator_e.push_back(new CAROM::BasisGenerator(*m_options_e[m_curr_win], isIncremental, basisName));
    ArrayCopynD(1,
                param->e_field,
                vec_wghosts.data(),
                sim[0].solver.dim_local,
                sim[0].solver.ghosts,
                0,
                index.data(),
                sim[0].solver.nvars);
    bool addeSample = m_generator_e[m_curr_win]->takeSample( vec_wghosts.data(), a_time, m_dt );
  } else {

    bool addSample = m_generator[m_curr_win]->takeSample( a_U.getData(), a_time, m_dt );
    char buffer[] = "sample";  // Creates a modifiable buffer and copies the string literal
    OutputlibROMfield(m_generator[0]->getSnapshotMatrix()->getColumn(m_tic)->getData(),
                      sim[0],
                      buffer);
    ArrayCopynD(1,
                param->potential,
                vec_wghosts.data(),
                sim[0].solver.dim_local,
                sim[0].solver.ghosts,
                0,
                index.data(),
                sim[0].solver.nvars);
    bool addphiSample = m_generator_phi[m_curr_win]->takeSample( vec_wghosts.data(), a_time, m_dt );

    ArrayCopynD(1,
                param->e_field,
                vec_wghosts.data(),
                sim[0].solver.dim_local,
                sim[0].solver.ghosts,
                0,
                index.data(),
                sim[0].solver.nvars);
    bool addeSample = m_generator_e[m_curr_win]->takeSample( vec_wghosts.data(), a_time, m_dt );
  }
//if (!m_rank) {
//  std::cout << "checking potential field\n";
//  for (int i = 0; i < param->npts_local_x_wghosts; i++) {
//        std::cout << "Index " << i << ": " << param->potential[i] << std::endl;
//  }
//  std::cout << std::endl;
//}
//if (!m_rank) {
//  std::cout << "checking e field\n";
//  for (int i = 0; i < param->npts_local_x_wghosts; i++) {
//        std::cout << "Index " << i << ": " << param->e_field[i] << std::endl;
//  }
//  std::cout << std::endl;
//}
  if (!m_rank) {
    std::cout << "checking e field (no ghosts) \n";
    for (int i = 0; i < vec_wghosts.size(); i++) {
          std::cout << "Index " << i << ": " << vec_wghosts[i] << std::endl;
    }
    std::cout << std::endl;
  }
  ::VlasovWriteSpatialField(&(sim[0].solver),&(sim[0].mpi),param->potential);
  ::FirstDerivativeSecondOrderCentralNoGhosts(vec_x_wghosts.data(),param->potential,0,1,&(sim[0].solver),&(sim[0].mpi));
  MPIExchangeBoundaries1D(mpi, vec_x_wghosts.data(), solver->dim_local[0], solver->ghosts, 0, solver->ndims);
  if (!m_rank) {
    std::cout << "checking grad phi (ghosts)\n";
    for (int i = 0; i < vec_x_wghosts.size(); i++) {
          std::cout << "Index " << i << ": " << -vec_x_wghosts[i]*solver->dxinv[0] << std::endl;
    }
    std::cout << std::endl;
  }
  ArrayCopynD(1,
              vec_x_wghosts.data(),
              vec_wghosts.data(),
              sim[0].solver.dim_local,
              sim[0].solver.ghosts,
              0,
              index.data(),
              sim[0].solver.nvars);
//if (!m_rank) {
//  std::cout << "checking -grad phi (no ghosts)\n";
//  for (int i = 0; i < vec_wghosts.size(); i++) {
//        std::cout << "Index " << i << ": " << -vec_wghosts[i]*solver->dxinv[0] << std::endl;
//  }
//  std::cout << std::endl;
//}
  ::SecondDerivativeSecondOrderCentralNoGhosts(vec_x_wghosts.data(),param->potential,0,&(sim[0].solver),&(sim[0].mpi));
  MPIExchangeBoundaries1D(mpi, vec_x_wghosts.data(), solver->dim_local[0], solver->ghosts, 0, solver->ndims);
//if (!m_rank) {
//  std::cout << "checking - nabla phi (ghosts)\n";
//  for (int i = 0; i < vec_x_wghosts.size(); i++) {
//        std::cout << "Index " << i << ": " << -vec_x_wghosts[i]*(solver->dxinv[0])*(solver->dxinv[0])<< std::endl;
//  }
//  std::cout << std::endl;
//}
//exit (0);
  printf("checking m_tic %d\n",m_tic);

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

  if (m_generator.size() > 0) {
    for (int i = 0; i < m_generator.size(); i++) {
      if (!m_ls_is_trained[i]) {
        int ncol = m_generator[i]->getSnapshotMatrix()->numColumns();
        if (!m_rank) {
          printf( "LSROMObject::train() - training LS object %d for sim. domain %d, var %d with %d samples.\n", 
                  m_curr_win, m_sim_idx, m_var_idx, ncol );
        }
        if (m_write_snapshot_mat) {
          char idx_string[_MAX_STRING_SIZE_];
          sprintf(idx_string, "%04d", i);
          std::string fname_root(m_dirname + "/snapshot_mat_"+std::string(idx_string));
          m_generator[i]->getSnapshotMatrix()->write(fname_root);
          /* Since here is using gertSnapshotMatrix to write files, need additional process in order to visualize it
             Similar in the case of DMD */
        }
        /* The code below is like m_dmd[i]->train(m_rdim) but in LS-ROM, there is no
           m_ls object, instead, it has m_generator, m_options, m_spatialbasis, m_projected_init */

        /* IMPORTANT!!! m_generator[i]->getSnapshotMatrix() is modified after
         * getSingularValues or (computeSVD) is called, hence need to make a copy of snapshots */
        /* Does everytime invoking getSpatialBasis() do computeSVD again? */
        m_snapshots = new CAROM::Matrix(m_generator[i]->getSnapshotMatrix()->getData(),
                                        m_generator[i]->getSnapshotMatrix()->numRows(),
                                        m_generator[i]->getSnapshotMatrix()->numColumns(),
                                        true,
                                        true);
        if (m_generator[i]->getSnapshotMatrix()->numColumns() < m_rdim) {
          throw std::runtime_error("# of snapshots is less than m_rdim");
        }
        m_S  = m_generator[i]->getSingularValues();

        OutputROMBasis(a_s, m_generator[0]->getSpatialBasis());
        ConstructROMHy(a_s, m_generator[0]->getSpatialBasis());
        CheckSolProjError(a_s);
        CheckHyProjError(a_s);

        printf("checking dimension of potential snapshot matrix:%d %d\n",m_generator_phi[i]->getSnapshotMatrix()->numRows(),m_generator_phi[i]->getSnapshotMatrix()->numColumns());
        m_snapshots_phi = new CAROM::Matrix(m_generator_phi[i]->getSnapshotMatrix()->getData(),
                                            m_generator_phi[i]->getSnapshotMatrix()->numRows(),
                                            m_generator_phi[i]->getSnapshotMatrix()->numColumns(),
                                            true,
                                            true);
        m_S_phi = m_generator_phi[i]->getSingularValues();
        OutputROMBasisPhi(a_s, m_generator_phi[0]->getSpatialBasis());
        ConstructPotentialROMRhs(a_s, m_generator[0]->getSpatialBasis(), m_generator_phi[0]->getSpatialBasis());
        ConstructPotentialROMLaplace(a_s, m_generator_phi[0]->getSpatialBasis());

        CheckPotentialProjError(a_s);
        CheckLaplaceProjError(a_s);
        CheckRhsProjError(a_s);
        m_snapshots_e = new CAROM::Matrix(m_generator_e[i]->getSnapshotMatrix()->getData(),
                                          m_generator_e[i]->getSnapshotMatrix()->numRows(),
                                          m_generator_e[i]->getSnapshotMatrix()->numColumns(),
                                          true,
                                          true);
        ConstructEBasis(a_s);
        CheckEProjError(a_s);
//      solver->HyperbolicFunction( solver->hyp,
//                                  m_snapshots->getColumn(0)->getData(),
//                                  solver,
//                                  mpi,
//                                  0,
//                                  1,
//                                  solver->FFunction,
//                                  solver->Upwind );
//      ::VlasovAdvection_x( solver->hyp, m_snapshots->getColumn(0)->getData(),
//                           0, solver, 0);
        ConstructROMHy_x(a_s, m_generator[0]->getSpatialBasis());
        ConstructROMHy_v(a_s, m_generator[0]->getSpatialBasis(), m_basis_e);
        exit(0);

      }
    }
  } else {
    printf("ERROR in LSROMObject::train(): m_generator is of size zero!");
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
  if (!m_rank) {
    std::cout << "----------------\n";
    std::cout << "Singular Values for potential: ";
    for (int i = 0; i < m_S_phi->dim(); i++) {
          std::cout << (m_S_phi->item(i)) << " ";
    }
    std::cout << "\n";
    std::cout << std::endl;
  }

  return;
}

/*! compute prediction at given time, the ROM coefficient is stored in m_romcoef  */
const CAROM::Vector* LSROMObject::predict(const double a_t, /*!< time at which to predict solution */
                                          void* a_s )
{
  SimulationObject* sim = (SimulationObject*) a_s;
  /* Need to modify the code so that is works for multiple windows */
  for (int i = 0; i < m_generator.size(); i++) {
    if (   (a_t >= m_intervals[i].first)
        && (  (a_t < m_intervals[i].second) || (m_intervals[i].second < 0)  ) ){

      if (!m_rank) printf("LSROM predicts at time t = %f in interval [%f, %f] \n",
                          a_t,m_intervals[i].first,m_intervals[i].second);
      std::vector<int> index(sim[0].solver.ndims);
      std::vector<double> vec_wghosts(sim[0].solver.npoints_local_wghosts*sim[0].solver.nvars);
      std::vector<double> rhs_wghosts(sim[0].solver.npoints_local_wghosts*sim[0].solver.nvars);
      int num_rows = m_generator[0]->getSpatialBasis()->numRows();
      CAROM::Vector* m_fomwork;
      CAROM::Vector* m_rhswork;
      m_fomwork = new CAROM::Vector(num_rows,false);
      m_rhswork = new CAROM::Vector(num_rows,true);

      if(std::abs(a_t) < 1e-9) {
        projectInitialSolution(*(m_snapshots->getColumn(0)));
        m_romcoef = new CAROM::Vector(m_projected_init[0]->getData(), m_rdim, false, true);
        m_fomwork = ReconlibROMfield(m_romcoef, m_generator[0]->getSpatialBasis(), m_rdim);
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
        TimeInitialize();
      } else {
        TimeRK(a_t,a_s);
      }
      if ((m_snap < sim[0].solver.n_iter) && (m_snap % m_sampling_freq == 0)){
        int idx;
        idx = m_snap/m_sampling_freq;
        m_fomwork = ReconlibROMfield(m_romcoef, m_generator[0]->getSpatialBasis(), m_rdim);
        ArrayCopynD(sim[0].solver.ndims,
                    m_fomwork->getData(),
                    vec_wghosts.data(),
                    sim[0].solver.dim_local,
                    0,
                    sim[0].solver.ghosts,
                    index.data(),
                    sim[0].solver.nvars);
        ArrayCopynD(sim[0].solver.ndims,
                    m_snapshots->getColumn(idx)->getData(),
                    rhs_wghosts.data(),
                    sim[0].solver.dim_local,
                    0,
                    sim[0].solver.ghosts,
                    index.data(),
                    sim[0].solver.nvars);
        char buffer[] = "reproderr";  // Creates a modifiable buffer and copies the string literal
        CalSnapROMDiff(&(sim[0].solver),&(sim[0].mpi),rhs_wghosts.data(),vec_wghosts.data(),buffer);
        printf("Reproduction error at iter %d, %.15f %.15f %.15f \n",idx,sim[0].solver.rom_diff_norms[0],sim[0].solver.rom_diff_norms[1],sim[0].solver.rom_diff_norms[2]);
      }
      m_snap++;
      return ReconlibROMfield(m_romcoef, m_generator[0]->getSpatialBasis(), m_rdim);
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

int LSROMObject::TimeInitialize()
{
  /* Currenty assuming one window only */
  int i;
  /* initialize arrays to NULL, then allocate as necessary */
  m_U.clear();
  m_Udot.clear();

  /* explicit Runge-Kutta methods */
  for (i = 0; i < nstages; i++) {
    m_U.push_back(new CAROM::Vector(m_rdim,false));
    m_Udot.push_back(new CAROM::Vector(m_rdim,false));
  }

  return(0);
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
  return(0);
}

int LSROMObject::TimeRK(const double a_t, /*!< time at which to predict solution */
                        void* a_s )
{
  /* Currenty assuming one window only */
  /* Advance the ROM ODE using RK4 scheme */

  SimulationObject* sim = (SimulationObject*) a_s;
  HyPar  *solver = (HyPar*) &(sim[0].solver);
  Vlasov *param  = (Vlasov*) solver->physics;
  MPIVariables *mpi = (MPIVariables *) param->m;
  std::vector<int> index(sim[0].solver.ndims);
  std::vector<double> vec_wghosts(sim[0].solver.npoints_local_wghosts*sim[0].solver.nvars);
  std::vector<double> rhs_wghosts(sim[0].solver.npoints_local_wghosts*sim[0].solver.nvars);

  int ns, stage, i;
  int num_rows = m_generator[0]->getSpatialBasis()->numRows();
  CAROM::Vector* m_fomwork;
  CAROM::Vector* m_rhswork;
  CAROM::Vector* m_romwork;
  m_fomwork = new CAROM::Vector(num_rows,false);
  m_rhswork = new CAROM::Vector(num_rows,true);
  m_romwork = new CAROM::Vector(m_rdim,false);
  CAROM::Vector* m_tmprhs;
  CAROM::Vector* m_tmpsol;
  CAROM::Vector* m_e;
  CAROM::Vector* m_contract1;
  CAROM::Vector* m_contract2;
  m_contract2 = new CAROM::Vector(m_rdim,false);

    /* Calculate stage values */
  for (stage = 0; stage < nstages; stage++) {
  
    double stagetime = a_t + c[stage]*m_dt;

    _ArrayCopy1D_(  m_romcoef->getData(),
                    m_U[stage]->getData(),
                    m_rdim );

    for (i = 0; i < stage; i++) {
      _ArrayAXPY_(  m_Udot[i]->getData(),
                    (m_dt * A[stage*nstages+i]),
                    m_U[stage]->getData(),
                    m_rdim );
    }

    if (m_direct_comp_hyperbolic) {
//    printf("compute hyperbolic term directly\n");

      m_fomwork = ReconlibROMfield(m_U[stage], m_generator[0]->getSpatialBasis(), m_rdim);

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

      m_Udot[stage] = ProjectToRB(m_rhswork,m_generator[0]->getSpatialBasis(), m_rdim);
//    if (!m_rank) {
//      std::cout << "Checking hyperbolic term [directly]: ";
//      for (int j = 0; j < m_rdim; j++) {
//            std::cout << m_Udot[stage]->item(j) << " ";
//      }
//      std::cout << std::endl;
//    }
    }
    else if (m_solve_phi) {
      m_tmprhs = m_romrhs_phi->mult(m_U[stage]);
      m_tmpsol = m_romlaplace_phi->mult(m_tmprhs);

      m_romwork = m_romhyperb_x->mult(m_U[stage]);
//    /* Tensor contraction */
      for (int i = 0; i < m_rdim; i++) {
        m_contract1 = m_romhyperb_v[i]->mult(m_U[stage]);
        m_contract2->item(i) = m_contract1->inner_product(m_tmpsol);
      }
      m_Udot[stage] = m_romwork->plus(m_contract2);
      printf("Done tensor contraction\n");
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
      if (!m_rank) {
        std::cout << "Checking hyperbolic term [phi]: ";
        for (int j = 0; j < m_rdim; j++) {
              std::cout << m_Udot[stage]->item(j) << " ";
        }
        std::cout << std::endl;
      }
    }
    else {
      m_Udot[stage] = m_romhyperb->mult(m_U[stage]);
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
                  m_romcoef->getData(),
                  m_rdim );

  }
  if (!m_rank) {
    std::cout << "Checking solved coefficient: ";
    for (int j = 0; j < m_rdim; j++) {
          std::cout << m_romcoef->item(j) << " ";
    }
    std::cout << std::endl;
  }
  delete m_fomwork;
  delete m_rhswork;
  delete m_romwork;
  delete m_contract2;
  return(0);
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

        projectInitialSolution(*(m_snapshots->getColumn(0)));
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
void LSROMObject::ConstructROMHy(void* a_s, const CAROM::Matrix* a_rombasis)
{
  SimulationObject* sim = (SimulationObject*) a_s;

  int num_rows = a_rombasis->numRows();
  int num_cols = a_rombasis->numColumns();

  /* Compute F(\Phi), where \Phi is the reduced basis matrix */
  CAROM::Matrix* phi_hyper;
  phi_hyper = new CAROM::Matrix(num_rows, m_rdim, true);

  CAROM::Vector phi_hyper_col(num_rows,false);
  std::vector<double> vec_wghosts(sim[0].solver.npoints_local_wghosts*sim[0].solver.nvars);
  std::vector<double> rhs_wghosts(sim[0].solver.npoints_local_wghosts*sim[0].solver.nvars);

  for (int j = 0; j < m_rdim; j++){
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
  printf("phi %d %d\n",a_rombasis->numRows(),a_rombasis->numColumns());
  printf("phi_hyper %d %d\n",phi_hyper->numRows(),phi_hyper->numColumns());
  m_romhyperb=a_rombasis->getFirstNColumns(m_rdim)->transposeMult(phi_hyper);
  printf("m_romhyperb %d %d\n",m_romhyperb->numRows(),m_romhyperb->numColumns());
//if (!m_rank) {
//  std::cout << "Checking ROM hyperbolic operator: \n";
//  for (int i = 0; i < m_romhyperb->numRows(); i++) {
//    for (int j = 0; j < m_romhyperb->numColumns(); j++) {
//        std::cout << (m_romhyperb->item(i,j)) << " ";
//    }
//  }
//  std::cout << std::endl;
//}
  ::ResetFilenameIndex( sim[0].solver.filename_index,
                        sim[0].solver.index_length );
  delete phi_hyper;

  return;
}

/*! Dump ROM basis */
void LSROMObject::OutputROMBasis(void* a_s, const CAROM::Matrix* a_rombasis)
{
  /* m_rdim is written out */
  SimulationObject* sim = (SimulationObject*) a_s;

  for (int j = 0; j < m_rdim; j++){
    char buffer[] = "basis";  // Creates a modifiable buffer and copies the string literal
    OutputlibROMfield(a_rombasis->getColumn(j)->getData(),
                      sim[0],
                      buffer);
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
  printf("dim_local %d %d\n",solver->dim_local[0],solver->dim_local[1]);
  printf("%d %d\n",solver->dim_local[0],solver->dim_local[1]);
  printf("npoints_global %d\n",solver->npoints_global);
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
  WriteArray( solver->ndims,
              solver->nvars,
              solver->dim_global,
              solver->dim_local,
              solver->ghosts,
              solver->x,
              u_diff,
              solver,
              mpi,
              filename);

  free(u_diff);
  return 0;
}

/*! Check projection error in solution */
void LSROMObject::CheckSolProjError(void* a_s)
{
  if (!m_rank) {
    std::cout << "----------------------------------------\n";
    std::cout << "Checking projection error in snapshots: ";
    std::cout << std::endl;
  }
  SimulationObject* sim = (SimulationObject*) a_s;
  std::vector<double> snap_wghosts(sim[0].solver.npoints_local_wghosts*sim[0].solver.nvars);
  std::vector<double> snapproj_wghosts(sim[0].solver.npoints_local_wghosts*sim[0].solver.nvars);
  std::vector<int> index(sim[0].solver.ndims);

  int num_rows = m_generator[0]->getSpatialBasis()->numRows();
  int num_cols = m_snapshots->numColumns();
  CAROM::Vector* recon_init;
  recon_init = new CAROM::Vector(num_rows,false);

  for (int j = 0; j < num_cols; j++){
    projectInitialSolution(*(m_snapshots->getColumn(j)));
    /* Extend reduced basis \phi_j with ghost points */
    ArrayCopynD(sim[0].solver.ndims,
                m_snapshots->getColumn(j)->getData(),
                snap_wghosts.data(),
                sim[0].solver.dim_local,
                0,
                sim[0].solver.ghosts,
                index.data(),
                sim[0].solver.nvars);

    recon_init = m_generator[0]->getSpatialBasis()->getFirstNColumns(m_rdim)->mult(m_projected_init[0]);
    ArrayCopynD(sim[0].solver.ndims,
                recon_init->getData(),
                snapproj_wghosts.data(),
                sim[0].solver.dim_local,
                0,
                sim[0].solver.ghosts,
                index.data(),
                sim[0].solver.nvars);
      char buffer[] = "snapprojerr";  // Creates a modifiable buffer and copies the string literal
      CalSnapROMDiff(&(sim[0].solver),&(sim[0].mpi),snap_wghosts.data(),snapproj_wghosts.data(),buffer);
      /* increment the index string, if required */
      if ((!strcmp(sim[0].solver.output_mode,"serial")) && (!strcmp(sim[0].solver.op_overwrite,"no"))) {
        ::IncrementFilenameIndex(sim[0].solver.filename_index,sim[0].solver.index_length);
      }
      printf("#%d snapshot projection error, %.15f %.15f %.15f \n",j,sim[0].solver.rom_diff_norms[0],sim[0].solver.rom_diff_norms[1],sim[0].solver.rom_diff_norms[2]);
  }
  ::ResetFilenameIndex( sim[0].solver.filename_index,
                        sim[0].solver.index_length );
  delete recon_init;
}


/*! Check projection error in hyperbolic term*/
void LSROMObject::CheckHyProjError(void* a_s)
{
  if (!m_rank) {
    std::cout << "----------------------------------------\n";
    std::cout << "Checking projection error in F(snapshots): ";
    std::cout << std::endl;
  }
  SimulationObject* sim = (SimulationObject*) a_s;
  std::vector<double> vec_wghosts(sim[0].solver.npoints_local_wghosts*sim[0].solver.nvars);
  std::vector<double> rhs_wghosts(sim[0].solver.npoints_local_wghosts*sim[0].solver.nvars);
  std::vector<double> snap_wghosts(sim[0].solver.npoints_local_wghosts*sim[0].solver.nvars);
  std::vector<double> snapproj_wghosts(sim[0].solver.npoints_local_wghosts*sim[0].solver.nvars);
  std::vector<int> index(sim[0].solver.ndims);

  int num_rows = m_generator[0]->getSpatialBasis()->numRows();
  int num_cols = m_snapshots->numColumns();
  CAROM::Vector* recon_init;
  CAROM::Vector* phi_colwork;
  recon_init = new CAROM::Vector(num_rows,false);
  phi_colwork = new CAROM::Vector(num_rows,true);

  for (int j = 0; j < num_cols; j++){
    projectInitialSolution(*(m_snapshots->getColumn(j)));
    /* Extend reduced basis \phi_j with ghost points */
    ArrayCopynD(sim[0].solver.ndims,
                m_snapshots->getColumn(j)->getData(),
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

//  /* Remove ghosts point in F(phi_j) */
//  ArrayCopynD(sim[0].solver.ndims,
//              rhs_wghosts.data(),
//              phi_colwork->getData(),
//              sim[0].solver.dim_local,
//              sim[0].solver.ghosts,
//              0,
//              index.data(),
//              sim[0].solver.nvars);
//  CAROM::Vector* m_working;
//  m_working = m_generator[0]->getSpatialBasis()->getFirstNColumns(m_rdim)->transposeMult(phi_colwork);
//  recon_init = m_generator[0]->getSpatialBasis()->getFirstNColumns(m_rdim)->mult(m_working);
//  ArrayCopynD(sim[0].solver.ndims,
//              recon_init->getData(),
//              snapproj_wghosts.data(),
//              sim[0].solver.dim_local,
//              0,
//              sim[0].solver.ghosts,
//              index.data(),
//              sim[0].solver.nvars);
    CAROM::Vector* m_working;
    m_working=m_romhyperb->mult(m_projected_init[0]);
    recon_init = m_generator[0]->getSpatialBasis()->getFirstNColumns(m_rdim)->mult(m_working);
    ArrayCopynD(sim[0].solver.ndims,
                recon_init->getData(),
                snapproj_wghosts.data(),
                sim[0].solver.dim_local,
                0,
                sim[0].solver.ghosts,
                index.data(),
                sim[0].solver.nvars);
//        const char* fnam2 = "recon_";
//        std::string fname2 = std::string(fnam2) + std::to_string(j);
//        char* fname_buffer2 = new char[fname2.length() + 1];
//        std::strcpy(fname_buffer2, fname2.c_str());

//        WriteArray( sim[0].solver.ndims,
//                    sim[0].solver.nvars,
//                    sim[0].solver.dim_global,
//                    sim[0].solver.dim_local,
//                    sim[0].solver.ghosts,
//                    sim[0].solver.x,
//                    snapproj_wghosts.data(),
//                    &(sim[0].solver),
//                    &(sim[0].mpi),
//                    fname_buffer2);
//        const char* fnam3 = "snap_col_";
//        std::string fname3 = std::string(fnam3) + std::to_string(j);
//        char* fname_buffer3 = new char[fname3.length() + 1];
//        std::strcpy(fname_buffer3, fname3.c_str());

//        WriteArray( sim[0].solver.ndims,
//                    sim[0].solver.nvars,
//                    sim[0].solver.dim_global,
//                    sim[0].solver.dim_local,
//                    sim[0].solver.ghosts,
//                    sim[0].solver.x,
//                    rhs_wghosts.data(),
//                    &(sim[0].solver),
//                    &(sim[0].mpi),
//                    fname_buffer3);
    char buffer[] = "hyprojerr";  // Creates a modifiable buffer and copies the string literal
    CalSnapROMDiff(&(sim[0].solver),&(sim[0].mpi),rhs_wghosts.data(),snapproj_wghosts.data(),buffer);
      /* increment the index string, if required */
      if ((!strcmp(sim[0].solver.output_mode,"serial")) && (!strcmp(sim[0].solver.op_overwrite,"no"))) {
        ::IncrementFilenameIndex(sim[0].solver.filename_index,sim[0].solver.index_length);
      }
    printf("F(#%d snapshot) projection error, %.15f %.15f %.15f \n",j,sim[0].solver.rom_diff_norms[0],sim[0].solver.rom_diff_norms[1],sim[0].solver.rom_diff_norms[2]);
  }
  ::ResetFilenameIndex( sim[0].solver.filename_index,
                        sim[0].solver.index_length );
  delete recon_init;
  delete phi_colwork;
}

/*! Construct potential ROM rhs */
void LSROMObject::ConstructPotentialROMRhs(void* a_s, const CAROM::Matrix* a_rombasis_f, const CAROM::Matrix* a_rombasis_phi)
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
  integral_basis_f = new CAROM::Matrix(num_rows, m_rdim, true);

  for (int j = 0; j < m_rdim; j++){
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
    ::VlasovWriteSpatialField(&(sim[0].solver),&(sim[0].mpi),int_f_wghosts.data(),buffer);
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
  m_romrhs_phi=a_rombasis_phi->getFirstNColumns(m_rdim_phi)->transposeMult(integral_basis_f);
  printf("f %d %d\n",a_rombasis_f->numRows(),a_rombasis_f->numColumns());
  printf("integral_basis_f %d %d\n",integral_basis_f->numRows(),integral_basis_f->numColumns());
  printf("m_romrhs_phi %d %d\n",m_romrhs_phi->numRows(),m_romrhs_phi->numColumns());
  if (!m_rank) {
    std::cout << "Checking Potential ROM Rhs: \n";
    for (int i = 0; i < m_romrhs_phi->numRows(); i++) {
      for (int j = 0; j < m_romrhs_phi->numColumns(); j++) {
          std::cout << (m_romrhs_phi->item(i,j)) << " ";
      }
    }
    std::cout << std::endl;
  }
  free(int_f);
  delete integral_basis_f;
  return;
}

/*! Construct potential ROM laplacian */
void LSROMObject::ConstructPotentialROMLaplace(void* a_s, const CAROM::Matrix* a_rombasis_phi)
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
  std::vector<int> idx(sim[0].solver.ndims);

  printf("a_rombasis_phi %d %d\n",a_rombasis_phi->numRows(),a_rombasis_phi->numColumns());

  /* Integrate f reduced basis over velocity */
  CAROM::Matrix* laplace_phi;
  laplace_phi = new CAROM::Matrix(num_rows, m_rdim_phi, true);
  CAROM::Vector laplace_col(num_rows,false);

  for (int j = 0; j < m_rdim_phi; j++){
    ArrayCopynD(1,
                a_rombasis_phi->getColumn(j)->getData(),
                vec_wghosts.data(),
                sim[0].solver.dim_local,
                0,
                sim[0].solver.ghosts,
                idx.data(),
                sim[0].solver.nvars);

    ::SecondDerivativeSecondOrderCentralNoGhosts(rhs_wghosts.data(),vec_wghosts.data(),0,&(sim[0].solver),&(sim[0].mpi));

    ArrayCopynD(1,
                rhs_wghosts.data(),
                laplace_col.getData(),
                sim[0].solver.dim_local,
                sim[0].solver.ghosts,
                0,
                idx.data(),
                sim[0].solver.nvars);

    /* Copy \int basis_f dv back to columns of integral_basis_f matrix */
    for (int i = 0; i < num_rows; i++) {
      (*laplace_phi)(i, j) = -1.0*laplace_col.getData()[i]*(solver->dxinv[0])*(solver->dxinv[0]);
    }
  }

  // construct rhs = basis_phi^T integral_basis_f
  m_romlaplace_phi=a_rombasis_phi->getFirstNColumns(m_rdim_phi)->transposeMult(laplace_phi);
  printf("m_romlaplace_phi %d %d\n",m_romlaplace_phi->numRows(),m_romlaplace_phi->numColumns());
  if (!m_rank) {
    std::cout << "Checking Potential ROM Laplace: \n";
    for (int i = 0; i < m_romlaplace_phi->numRows(); i++) {
      for (int j = 0; j < m_romlaplace_phi->numColumns(); j++) {
          std::cout << (m_romlaplace_phi->item(i,j)) << " ";
      }
    }
    std::cout << std::endl;
  }
  m_romlaplace_phi->inverse();
//if (!m_rank) {
//  std::cout << "Checking Potential ROM inverse Laplace: \n";
//  for (int i = 0; i < m_romlaplace_phi->numRows(); i++) {
//    for (int j = 0; j < m_romlaplace_phi->numColumns(); j++) {
//        std::cout << (m_romlaplace_phi->item(i,j)) << " ";
//    }
//  }
//  std::cout << std::endl;
//}
  return;
}

/*! Check projection error in potential solution */
void LSROMObject::CheckPotentialProjError(void* a_s)
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

  int num_rows = m_generator_phi[0]->getSpatialBasis()->numRows();
  int num_cols = m_snapshots_phi->numColumns();
  CAROM::Vector* recon_caromvec;
  recon_caromvec= new CAROM::Vector(num_rows,false);

  for (int j = 0; j < num_cols; j++){
    /* Extend reduced basis \phi_j with ghost points */
    ArrayCopynD(1,
                m_snapshots_phi->getColumn(j)->getData(),
                snap_wghosts.data(),
                sim[0].solver.dim_local,
                0,
                sim[0].solver.ghosts,
                index.data(),
                sim[0].solver.nvars);

    projectInitialSolution_phi(*(m_snapshots_phi->getColumn(j)));
    recon_caromvec= m_generator_phi[0]->getSpatialBasis()->getFirstNColumns(m_rdim_phi)->mult(m_projected_init_phi[0]);
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
    printf("#%d snapshot projection error phi, %.15f %.15f %.15f \n",j,sim[0].solver.rom_diff_norms[0],sim[0].solver.rom_diff_norms[1],sim[0].solver.rom_diff_norms[2]);
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

    m_projected_init_phi.push_back(new CAROM::Vector(m_rdim_phi, false));
    m_projected_init_phi[i] = m_generator_phi[i]->getSpatialBasis()->getFirstNColumns(m_rdim_phi)->transposeMult(a_U);

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
void LSROMObject::CheckLaplaceProjError(void* a_s)
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

  int num_rows = m_generator_phi[0]->getSpatialBasis()->numRows();
  CAROM::Vector* recon_field;
  CAROM::Vector* recon_lap;
  CAROM::Vector* field_reduced;
  CAROM::Vector* lap_reduced;

  int num_cols = m_snapshots_phi->numColumns();
  char buffer[] = "laplace_projerr";  // Creates a modifiable buffer and copies the string literal

  for (int j = 0; j < num_cols; j++){
    /* extend snapshot with ghost points */
    ArrayCopynD(1,
                m_snapshots_phi->getColumn(j)->getData(),
                snap_wghosts.data(),
                sim[0].solver.dim_local,
                0,
                sim[0].solver.ghosts,
                index.data(),
                sim[0].solver.nvars);

    ::SecondDerivativeSecondOrderCentralNoGhosts(lap_snap_wghosts.data(),snap_wghosts.data(),0,&(sim[0].solver),&(sim[0].mpi));
    _ArrayScale1D_(lap_snap_wghosts,-1.0*(solver->dxinv[0])*(solver->dxinv[0]),param->npts_local_x_wghosts);
//  std::cout << "Checking laplacian of snap: ";
//  for (int i = 0; i < param->npts_local_x_wghosts ; i++) {
//        std::cout << (lap_snap_wghosts.data()[i]) << " ";
//  }
//  std::cout << std::endl;

    /* Check 1 */
    /* project snapshot onto reduced space */
    projectInitialSolution_phi(*(m_snapshots_phi->getColumn(j)));
    lap_reduced = m_romlaplace_phi->mult(m_projected_init_phi[0]);
    recon_lap = m_generator_phi[0]->getSpatialBasis()->getFirstNColumns(m_rdim_phi)->mult(lap_reduced);
    ArrayCopynD(1,
                recon_lap->getData(),
                recon_wghosts.data(),
                sim[0].solver.dim_local,
                0,
                sim[0].solver.ghosts,
                index.data(),
                sim[0].solver.nvars);
    CalSnapROMDiff_phi(&(sim[0].solver),&(sim[0].mpi),lap_snap_wghosts.data(),recon_wghosts.data(),buffer);
    printf("#%d snapshot projection error in laplacian (reduced), %.15f %.15f %.15f \n",j,sim[0].solver.rom_diff_norms[0],sim[0].solver.rom_diff_norms[1],sim[0].solver.rom_diff_norms[2]);

    /* Check 2 */
    recon_field = m_generator_phi[0]->getSpatialBasis()->getFirstNColumns(m_rdim_phi)->mult(m_projected_init_phi[0]);
    ArrayCopynD(1,
                recon_field->getData(),
                snapproj_wghosts.data(),
                sim[0].solver.dim_local,
                0,
                sim[0].solver.ghosts,
                index.data(),
                sim[0].solver.nvars);
    ::SecondDerivativeSecondOrderCentralNoGhosts(lap_snapproj_wghosts.data(),snapproj_wghosts.data(),0,&(sim[0].solver),&(sim[0].mpi));
    _ArrayScale1D_(lap_snapproj_wghosts,-1.0*(solver->dxinv[0])*(solver->dxinv[0]),param->npts_local_x_wghosts);
    ArrayCopynD(1,
                lap_snapproj_wghosts.data(),
                recon_field->getData(),
                sim[0].solver.dim_local,
                sim[0].solver.ghosts,
                0,
                index.data(),
                sim[0].solver.nvars);

    field_reduced = m_generator_phi[0]->getSpatialBasis()->getFirstNColumns(m_rdim_phi)->transposeMult(recon_field);
    recon_lap = m_generator_phi[0]->getSpatialBasis()->getFirstNColumns(m_rdim_phi)->mult(field_reduced);
    ArrayCopynD(1,
                recon_lap->getData(),
                recon_wghosts.data(),
                sim[0].solver.dim_local,
                0,
                sim[0].solver.ghosts,
                index.data(),
                sim[0].solver.nvars);
//    std::cout << "Checking recon_wghosts: ";
//    for (int i = 0; i < param->npts_local_x_wghosts ; i++) {
//          std::cout << (recon_wghosts.data()[i]) << " ";
//    }
//    std::cout << std::endl;
    CalSnapROMDiff_phi(&(sim[0].solver),&(sim[0].mpi),lap_snap_wghosts.data(),recon_wghosts.data(),buffer);
    printf("#%d snapshot projection error in laplacian (non reduced), %.15f %.15f %.15f \n",j,sim[0].solver.rom_diff_norms[0],sim[0].solver.rom_diff_norms[1],sim[0].solver.rom_diff_norms[2]);

    /* Check 3 */
    ArrayCopynD(1,
                lap_snap_wghosts.data(),
                recon_lap->getData(),
                sim[0].solver.dim_local,
                sim[0].solver.ghosts,
                0,
                index.data(),
                sim[0].solver.nvars);

    lap_reduced = m_generator_phi[0]->getSpatialBasis()->getFirstNColumns(m_rdim_phi)->transposeMult(recon_lap);
    recon_lap = m_generator_phi[0]->getSpatialBasis()->getFirstNColumns(m_rdim_phi)->mult(lap_reduced);
    ArrayCopynD(1,
                recon_lap->getData(),
                recon_wghosts.data(),
                sim[0].solver.dim_local,
                0,
                sim[0].solver.ghosts,
                index.data(),
                sim[0].solver.nvars);
    CalSnapROMDiff_phi(&(sim[0].solver),&(sim[0].mpi),lap_snap_wghosts.data(),recon_wghosts.data(),buffer);
    printf("#%d snapshot projection error in laplacian (simply project), %.15f %.15f %.15f \n",j,sim[0].solver.rom_diff_norms[0],sim[0].solver.rom_diff_norms[1],sim[0].solver.rom_diff_norms[2]);

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
void LSROMObject::CheckRhsProjError(void* a_s)
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

  int num_rows = m_generator_phi[0]->getSpatialBasis()->numRows();
  CAROM::Vector* recon_rhs;
  recon_rhs = new CAROM::Vector(num_rows,false);
  CAROM::Vector* int_f_caromvec;
  int_f_caromvec = new CAROM::Vector(num_rows,true);
  CAROM::Vector* rhs_reduced;

  double* int_f= (double*) calloc(num_rows, sizeof(double));
  int num_cols = m_snapshots->numColumns();

  for (int j = 0; j < num_cols; j++){

    EvaluatePotentialRhs(a_s, m_snapshots->getColumn(j), int_f);

    ArrayCopynD(1,
                int_f,
                int_f_wghosts.data(),
                sim[0].solver.dim_local,
                0,
                sim[0].solver.ghosts,
                index.data(),
                sim[0].solver.nvars);
    char buffer1[] = "int_f";
    ::VlasovWriteSpatialField(&(sim[0].solver),&(sim[0].mpi),int_f_wghosts.data(),buffer1);

    ArrayCopynD(1,
                int_f,
                int_f_caromvec->getData(),
                sim[0].solver.dim_local,
                0,
                0,
                index.data(),
                sim[0].solver.nvars);

    /* Check 1 */
    rhs_reduced = m_generator_phi[0]->getSpatialBasis()->getFirstNColumns(m_rdim_phi)->transposeMult(int_f_caromvec);
    if (!m_rank) {
      std::cout << "Checking rhs reduced: ";
      for (int j = 0; j < rhs_reduced->dim(); j++) {
            std::cout << (rhs_reduced->item(j)) << " ";
      }
      std::cout << std::endl;
    }
    recon_rhs = m_generator_phi[0]->getSpatialBasis()->getFirstNColumns(m_rdim_phi)->mult(rhs_reduced);
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
    printf("#%d snapshot projection error in rhs (non reduced), %.15f %.15f %.15f \n",j,sim[0].solver.rom_diff_norms[0],sim[0].solver.rom_diff_norms[1],sim[0].solver.rom_diff_norms[2]);

    /* Check 2 */
    projectInitialSolution(*(m_snapshots->getColumn(j)));
    rhs_reduced = m_romrhs_phi->mult(m_projected_init[0]);
    recon_rhs = m_generator_phi[0]->getSpatialBasis()->getFirstNColumns(m_rdim_phi)->mult(rhs_reduced);
    ArrayCopynD(1,
                recon_rhs->getData(),
                recon_wghosts.data(),
                sim[0].solver.dim_local,
                0,
                sim[0].solver.ghosts,
                index.data(),
                sim[0].solver.nvars);
    CalSnapROMDiff_phi(&(sim[0].solver),&(sim[0].mpi),int_f_wghosts.data(),recon_wghosts.data(),buffer);
    printf("#%d snapshot projection error in rhs (reduced), %.15f %.15f %.15f \n",j,sim[0].solver.rom_diff_norms[0],sim[0].solver.rom_diff_norms[1],sim[0].solver.rom_diff_norms[2]);


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
void LSROMObject::CheckEProjError(void* a_s)
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

  int num_rows = m_generator_phi[0]->getSpatialBasis()->numRows();
  int num_cols = m_snapshots_e->numColumns();
  CAROM::Vector* recon_caromvec;

  for (int j = 0; j < num_cols; j++){
    ArrayCopynD(1,
                m_snapshots_e->getColumn(j)->getData(),
                snap_wghosts.data(),
                sim[0].solver.dim_local,
                0,
                sim[0].solver.ghosts,
                index.data(),
                sim[0].solver.nvars);

    projectInitialSolution_phi(*(m_snapshots_phi->getColumn(j)));
    recon_caromvec = m_basis_e->mult(m_projected_init_phi[0]);
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
    printf("#%d snapshot projection error e, %.15f %.15f %.15f \n",j,sim[0].solver.rom_diff_norms[0],sim[0].solver.rom_diff_norms[1],sim[0].solver.rom_diff_norms[2]);

  }
  ::ResetFilenameIndex( sim[0].solver.filename_index,
                        sim[0].solver.index_length );
  delete recon_caromvec;
}

/*! Dump ROM basis phi */
void LSROMObject::OutputROMBasisPhi(void* a_s, const CAROM::Matrix* a_rombasis)
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
    ::VlasovWriteSpatialField(&(sim[0].solver),&(sim[0].mpi),vec_wghosts.data(),buffer);

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
void LSROMObject::ConstructEBasis(void* a_s)
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

  int num_rows = m_generator_phi[0]->getSpatialBasis()->numRows();
  double* vec_noghosts = (double*) calloc(num_rows, sizeof(double));

  m_basis_e = new CAROM::Matrix(num_rows, m_rdim_phi, true);

  for (int j = 0; j < m_rdim_phi; j++){
    ArrayCopynD(1,
                m_generator_phi[0]->getSpatialBasis()->getColumn(j)->getData(),
                basis_vec_wghosts.data(),
                sim[0].solver.dim_local,
                0,
                sim[0].solver.ghosts,
                index.data(),
                sim[0].solver.nvars);
    ::FirstDerivativeSecondOrderCentralNoGhosts(vec_x_wghosts.data(),basis_vec_wghosts.data(),0,1,&(sim[0].solver),&(sim[0].mpi));
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
      (*m_basis_e)(i, j) = vec_noghosts[i];
    }

  }
  free(vec_noghosts);
}

/*! Construct reduced hyperbolic operator in x direction */
void LSROMObject::ConstructROMHy_x(void* a_s, const CAROM::Matrix* a_rombasis)
{
  SimulationObject* sim = (SimulationObject*) a_s;
  HyPar  *solver = (HyPar*) &(sim[0].solver);
  Vlasov *param  = (Vlasov*) solver->physics;
  MPIVariables *mpi = (MPIVariables *) param->m;

  int num_rows = a_rombasis->numRows();
  int num_cols = a_rombasis->numColumns();

  /* Compute F(\Phi), where \Phi is the reduced basis matrix */
  CAROM::Matrix* phi_hyper_x;
  phi_hyper_x = new CAROM::Matrix(num_rows, m_rdim, true);
  CAROM::Vector phi_hyper_col(num_rows,false);

  std::vector<double> vec_wghosts(sim[0].solver.npoints_local_wghosts*sim[0].solver.nvars);
  std::vector<double> vec1_wghosts(sim[0].solver.npoints_local_wghosts*sim[0].solver.nvars);
  std::vector<double> rhs_wghosts(sim[0].solver.npoints_local_wghosts*sim[0].solver.nvars);
  std::vector<int> index(sim[0].solver.ndims);

  for (int j = 0; j < m_rdim; j++){
    ArrayCopynD(sim[0].solver.ndims,
                a_rombasis->getColumn(j)->getData(),
                vec_wghosts.data(),
                sim[0].solver.dim_local,
                0,
                sim[0].solver.ghosts,
                index.data(),
                sim[0].solver.nvars);
    ArrayCopynD(sim[0].solver.ndims,
                a_rombasis->getColumn(j)->getData(),
                vec1_wghosts.data(),
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
//  printf("checking Hyperbolic\n");
//  solver->ApplyBoundaryConditions(solver,mpi,vec1_wghosts.data(),NULL,0);
//  solver->ApplyIBConditions(solver,mpi,vec1_wghosts.data(),0);
//  MPIExchangeBoundariesnD(  solver->ndims,
//                            solver->nvars,
//                            solver->dim_local,
//                            solver->ghosts,
//                            mpi,
//                            vec1_wghosts.data());
//      solver->HyperbolicFunction( solver->hyp,
//                                  vec1_wghosts.data(),
//                                  solver,
//                                  mpi,
//                                  0,
//                                  1,
//                                  solver->FFunction,
//                                  solver->Upwind );
//  printf("checking rhs_wghosts.data()\n");
//  for (int k = 0; k < rhs_wghosts.size(); k++) {
//    printf("k: %d,  value: %.15f\n",k,rhs_wghosts[k]);
//  }
//  exit (0);

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
  m_romhyperb_x=a_rombasis->getFirstNColumns(m_rdim)->transposeMult(phi_hyper_x);
  printf("m_romhyperb_x %d %d\n",m_romhyperb_x->numRows(),m_romhyperb_x->numColumns());
  delete phi_hyper_x;

  return;
}

/*! Construct reduced hyperbolic tensor in v direction */
void LSROMObject::ConstructROMHy_v(void* a_s, const CAROM::Matrix* a_rombasis, const CAROM::Matrix* a_rombasis_phi)
{
  SimulationObject* sim = (SimulationObject*) a_s;
  HyPar  *solver = (HyPar*) &(sim[0].solver);
  Vlasov *param  = (Vlasov*) solver->physics;
  MPIVariables *mpi = (MPIVariables *) param->m;

  int num_rows = a_rombasis->numRows();
  int num_cols = a_rombasis->numColumns();

  /* Compute F(\Phi), where \Phi is the reduced basis matrix */
  CAROM::Vector phi_hyper_col(num_rows,false);

  std::vector<int> index(sim[0].solver.ndims);
  std::vector<double> vec_wghosts(sim[0].solver.npoints_local_wghosts*sim[0].solver.nvars);
  std::vector<double> rhs_wghosts(sim[0].solver.npoints_local_wghosts*sim[0].solver.nvars);

  for (int i = 0; i < m_rdim; i++) {

    m_romhyperb_v.push_back(new CAROM::Matrix(m_rdim_phi, m_rdim, true));
    for (int j = 0; j < m_rdim_phi; j++) {

      for (int k = 0; k < m_rdim; k++) {
        ArrayCopynD(sim[0].solver.ndims,
                    a_rombasis->getColumn(k)->getData(),
                    vec_wghosts.data(),
                    sim[0].solver.dim_local,
                    0,
                    sim[0].solver.ghosts,
                    index.data(),
                    sim[0].solver.nvars);

        ArrayCopynD(1,
                    a_rombasis_phi->getColumn(j)->getData(),
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

        (*m_romhyperb_v[i])(j, k) = a_rombasis->getColumn(i)->inner_product(phi_hyper_col);
//      printf("checking m_romhyperb_v %d %d %f\n",j,k,(*m_romhyperb_v[i])(j, k));
      }
    }
    printf("matrix size %d %d\n",m_romhyperb_v[i]->numRows(),m_romhyperb_v[i]->numColumns());
  }
    printf("i size %d\n",m_romhyperb_v.size());

  return;
}

#endif

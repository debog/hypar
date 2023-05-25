#ifdef with_librom

/*! @file LSROMObject.cpp
 *  @brief Member functions of the class #LSROMObject
 *  @author Ping-Hsuan Tsai
*/

#include <basic.h>
#include <rom_object_ls.h>
#include <simulation_object.h>
//#include <timeintegration.h>
#include <arrayfunctions.h>
#include <io_cpp.h>
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
extern "C" int TimeRHSFunctionExplicit(double*,double*,void*,void*,double);
extern "C" int CalculateROMDiff(void*,void*);

LSROMObject::LSROMObject( const int     a_vec_size, /*!< vector size */
                            const double  a_dt,       /*!< time step size */
                            const int     a_rdim,     /*!< latent space dimension */
                            const int     a_rank,     /*!< MPI rank of this process */
                            const int     a_nproc,    /*!< Number of MPI processes */
                            const int     a_sim_idx,  /*!< Simulation index (for ensemble simulations */
                            const int     a_var_idx   /*!< Vector component index */ )
{
  m_rank = a_rank;
  m_nproc = a_nproc;

/* my implementation */
  m_options.clear();
  m_generator.clear();
  m_projected_init.clear();

  m_ls_is_trained.clear();
  m_intervals.clear();
  m_vec_size = a_vec_size;
  m_dt = a_dt;
  m_t_final = -1;
  m_rdim = a_rdim;

  m_sim_idx = a_sim_idx;
  m_var_idx = a_var_idx;

  m_num_window_samples = INT_MAX;

  char dirname_c_str[_MAX_STRING_SIZE_] = "LS";
  char write_snapshot_mat[_MAX_STRING_SIZE_] = "false";

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
    printf("  directory name for LS objects: %s\n", dirname_c_str);
    printf("  write snapshot matrix to file:  %s\n", write_snapshot_mat);
    if (m_sim_idx >= 0) {
      printf("  simulation domain:  %d\n", m_sim_idx);
    }
    if (m_var_idx >= 0) {
      printf("  Vector component:  %d\n", m_var_idx);
    }
  }

#ifndef serial
  MPI_Bcast(&m_num_window_samples,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(dirname_c_str,_MAX_STRING_SIZE_,MPI_CHAR,0,MPI_COMM_WORLD);
  MPI_Bcast(write_snapshot_mat,_MAX_STRING_SIZE_,MPI_CHAR,0,MPI_COMM_WORLD);
#endif

  m_dirname = std::string( dirname_c_str );
  m_write_snapshot_mat = (std::string(write_snapshot_mat) == "true");

  if (m_num_window_samples <= m_rdim) {
    printf("ERROR:LSROMObject::LSROMObject() - m_num_window_samples <= m_rdim!!");
  }

  m_tic = 0;
  m_curr_win = 0;
}

void LSROMObject::projectInitialSolution(  CAROM::Vector& a_U /*!< solution vector */ )
{
  /* Assumming m_generator[i]->getSpatialBasis() is called before */
  if (m_generator.size() == 0) {
    if (!m_rank) {
      printf("ERROR in LSROMObject::projectInitialSolution() - m_ls is a vector of size 0.\n");
    }
    return;
  }
  for (int i = 0; i < m_generator.size(); i++) {
    int num_rows = m_generator[i]->getSpatialBasis()->numRows();
    m_projected_init.push_back(new CAROM::Vector(num_rows, false));
    m_projected_init[i] = m_generator[i]->getSpatialBasis()->transposeMult(a_U);
    printf("m_project size %d\n",m_projected_init[i]->dim());
    if (!m_rank) {
      std::cout << "Checking initial condition: ";
      for (int j = 0; j < m_projected_init[i]->dim(); j++) {
            std::cout << (m_projected_init[i]->item(j)) << " ";
      }
      std::cout << std::endl;
    }
  }

  /* Currently only support for one time window
   * Following code is copied from DMDObject for multiple windows and need to
   * be modified accordingly for LS-ROM */
//    m_ls[0]->projectInitialCondition( &a_U );
//    for (int i = 1; i < m_ls.size(); i++) {
//      m_ls[i]->projectInitialCondition( m_ls[i-1]->predict(m_intervals[i].first) );
//    }
//    return;
}

/*! take a sample (solution snapshot) */
void LSROMObject::takeSample(  const CAROM::Vector& a_U, /*!< solution vector */
                                const double a_time, /*!< sample time */ 
                                void* a_s )
{
  SimulationObject* sim = (SimulationObject*) a_s;
  if (m_tic == 0) {

    m_options.push_back(new CAROM::Options(m_vec_size, max_num_snapshots, 1, update_right_SV));
    m_generator.push_back(new CAROM::BasisGenerator(*m_options[m_curr_win], isIncremental, basisName));
//  printf("HEREEEEEE %d \n", isIncremental ? 1 : 0);
//  exit (0);
    m_intervals.push_back( Interval(a_time, m_t_final) );
    m_ls_is_trained.push_back(false);

    if (!m_rank) {
      printf( "LSROMObject::takeSample() - creating new generator object for sim. domain %d, var %d, t=%f (total: %d).\n",
              m_sim_idx, m_var_idx, m_intervals[m_curr_win].first, m_generator.size());
    }
    /* QUESTION: should a_U be centered? */
    bool addSample = m_generator[m_curr_win]->takeSample( a_U.getData(), a_time, m_dt );

    double* vec_data = m_generator[0]->getSnapshotMatrix()->getColumn(0)->getData();

    WriteArray( sim[0].solver.ndims,
                sim[0].solver.nvars,
                sim[0].solver.dim_global,
                sim[0].solver.dim_local,
                0,
                sim[0].solver.x,
                vec_data,
                &(sim[0].solver),
                &(sim[0].mpi),
                "sample_" );
    // In order to use WriteArray function in hypar,
    // Instead of using writeSnapshot(), use getSnapshotMatrix if one likes to visualize the snapshot,
    // Use getspatial basis if one likes to visualize the basis.
    // In both cases, the question is how to use Matrix datatype with WriteArray.
  } else {

    bool addSample = m_generator[m_curr_win]->takeSample( a_U.getData(), a_time, m_dt );

    printf("filename_index %s\n",sim[0].solver.filename_index);
    double* vec_data = m_generator[0]->getSnapshotMatrix()->getColumn(m_tic)->getData();
    WriteArray( sim[0].solver.ndims,
                sim[0].solver.nvars,
                sim[0].solver.dim_global,
                sim[0].solver.dim_local,
                0,
                sim[0].solver.x,
                vec_data,
                &(sim[0].solver),
                &(sim[0].mpi),
                "sample_" );
  }

  m_tic++;
  return;
}

/*! train the LS objects */
void LSROMObject::train(void* a_s)
{
  /* In the training, initial condition for ROM should also be computed */
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

        /* Call SVD on snapshot matrix and also check singular value decomposition
           Does everytime invoking getSpatialBasis() do computeSVD again? */
//      const CAROM::Matrix* snapshots = new CAROM::Matrix(*m_generator[i]->getSnapshotMatrix());
        const CAROM::Matrix* snapshots = new CAROM::Matrix(m_generator[i]->getSnapshotMatrix()->getData(),
                                                           m_generator[i]->getSnapshotMatrix()->numRows(),
                                                           m_generator[i]->getSnapshotMatrix()->numColumns(),
                                                           true,
                                                           true);
        m_snapshots = new CAROM::Matrix(m_generator[i]->getSnapshotMatrix()->getData(),
                                                           m_generator[i]->getSnapshotMatrix()->numRows(),
                                                           m_generator[i]->getSnapshotMatrix()->numColumns(),
                                                           true,
                                                           true);

        /* IMPORTANT!!! m_generator[i]->getSnapshotMatrix() is modified after getSingularValues or (computeSVD) is called */ 
        const CAROM::Vector* sing_vals = m_generator[i]->getSingularValues();
        if (!m_rank) {
          std::cout << "Singular Values: ";
          for (int i = 0; i < sing_vals->dim(); i++) {
                std::cout << (sing_vals->item(i)) << " ";
          }
          std::cout << std::endl;
        }

        /* Compute F(\Phi), where \Phi is the reduced basis matrix */
        int num_rows = m_generator[i]->getSpatialBasis()->numRows();
        int num_cols = m_generator[i]->getSpatialBasis()->numColumns();
        m_rdim = num_cols;

        CAROM::Matrix* phi_hyper;
        phi_hyper = new CAROM::Matrix(num_rows, num_cols, true);
        printf( "Check phi_hyper rows, cols: %d %d \n",num_rows,num_cols);

        CAROM::Vector phi_hyper_col(num_rows,false);
        std::vector<double> vec_wghosts(sim[0].solver.npoints_local_wghosts*sim[0].solver.nvars);
        std::vector<double> rhs_wghosts(sim[0].solver.npoints_local_wghosts*sim[0].solver.nvars);

        for (int j = 0; j < num_cols; j++){
          /* Extend reduced basis \phi_j with ghost points */
          std::vector<int> index(sim[0].solver.ndims);
          ArrayCopynD(sim[0].solver.ndims,
                      m_generator[0]->getSpatialBasis()->getColumn(j)->getData(),
                      vec_wghosts.data(),
                      sim[0].solver.dim_local,
                      0,
                      sim[0].solver.ghosts,
                      index.data(),
                      sim[0].solver.nvars);

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

          const char* fname = "hyper_";
          std::string fname1 = std::string(fname) + std::to_string(j);
          char* fname_buffer = new char[fname1.length() + 1];
          std::strcpy(fname_buffer, fname1.c_str());

          WriteArray( sim[0].solver.ndims,
                      sim[0].solver.nvars,
                      sim[0].solver.dim_global,
                      sim[0].solver.dim_local,
                      0,
                      sim[0].solver.x,
                      phi_hyper->getColumn(j)->getData(),
                      &(sim[0].solver),
                      &(sim[0].mpi),
                      fname_buffer);

          const char* filename = "basis_";
          std::string file_name = std::string(filename) + std::to_string(j);
          char* file_name_buffer = new char[file_name.length() + 1];
          std::strcpy(file_name_buffer, file_name.c_str());

          // ghost is 0, x and y coordinate is not correctly outpost
          WriteArray( sim[0].solver.ndims,
                      sim[0].solver.nvars,
                      sim[0].solver.dim_global,
                      sim[0].solver.dim_local,
                      0,
                      sim[0].solver.x,
                      m_generator[0]->getSpatialBasis()->getColumn(j)->getData(),
                      &(sim[0].solver),
                      &(sim[0].mpi),
                      file_name_buffer);
        }

        // construct hyper_ROM = phi^T phi_hyper
        printf("phi %d %d\n",m_generator[i]->getSpatialBasis()->numRows(),m_generator[i]->getSpatialBasis()->numColumns());
        printf("phi_hyper %d %d\n",phi_hyper->numRows(),phi_hyper->numColumns());
        m_romhyperb=m_generator[0]->getSpatialBasis()->transposeMult(phi_hyper);
        printf("m_romhyperb %d %d\n",m_romhyperb->numRows(),m_romhyperb->numColumns());
        if (!m_rank) {
          std::cout << "Checking ROM operator: \n";
          for (int i = 0; i < m_romhyperb->numRows(); i++) {
            for (int j = 0; j < m_romhyperb->numColumns(); j++) {
                std::cout << (m_romhyperb->item(i,j)) << " ";
            }
          }
          std::cout << std::endl;
        }
        /* Below is just checking if the basis function is orthogonal */
//      projectInitialSolution(*(m_generator[0]->getSpatialBasis()->getColumn(0)));
//      projectInitialSolution(*(m_generator[0]->getSpatialBasis()->getColumn(1)));
//      projectInitialSolution(*(m_generator[0]->getSpatialBasis()->getColumn(2)));
//      projectInitialSolution(*(m_generator[0]->getSpatialBasis()->getColumn(3)));

        /* Compute initial condition */
        /* Can't simply pass *(m_generator[0]->getSnapshotMatrix()->getColumn(0))) since
         * m_generator[0]->getSpatialBasis() is distributed whereas m_generator[0]->getSnapshotMatrix()
         * is not distributed */
//      CAROM::Vector* snap_col;
//      snap_col = new CAROM::Vector(snapshots->getColumn(0)->getData(),
//                                   m_generator[0]->getSnapshotMatrix()->numRows(),true);
//      projectInitialSolution(*snap_col);
        projectInitialSolution(*(m_snapshots->getColumn(0)));

        CAROM::Vector* recon_init;
        recon_init = new CAROM::Vector(num_rows,false);
        recon_init = m_generator[0]->getSpatialBasis()->mult(m_projected_init[0]);

        std::vector<int> index(sim[0].solver.ndims);
        ArrayCopynD(sim[0].solver.ndims,
                    recon_init->getData(),
                    sim[0].solver.u_rom_predicted,
                    sim[0].solver.dim_local,
                    0,
                    sim[0].solver.ghosts,
                    index.data(),
                    sim[0].solver.nvars);
          const char* fnam2 = "recon_";
          std::string fname2 = std::string(fnam2) + std::to_string(0);
          char* fname_buffer2 = new char[fname2.length() + 1];
          std::strcpy(fname_buffer2, fname2.c_str());

          WriteArray( sim[0].solver.ndims,
                      sim[0].solver.nvars,
                      sim[0].solver.dim_global,
                      sim[0].solver.dim_local,
                      0,
                      sim[0].solver.x,
                      recon_init->getData(),
                      &(sim[0].solver),
                      &(sim[0].mpi),
                      fname_buffer2);
        ArrayCopynD(sim[0].solver.ndims,
                    m_snapshots->getColumn(9)->getData(),
                    sim[0].solver.u,
                    sim[0].solver.dim_local,
                    0,
                    sim[0].solver.ghosts,
                    index.data(),
                    sim[0].solver.nvars);
          const char* fnam3 = "snap_col_";
          std::string fname3 = std::string(fnam3) + std::to_string(0);
          char* fname_buffer3 = new char[fname3.length() + 1];
          std::strcpy(fname_buffer3, fname3.c_str());

//        double* vec_data = m_generator[0]->getSnapshotMatrix()->getColumn(0)->getData();
          WriteArray( sim[0].solver.ndims,
                      sim[0].solver.nvars,
                      sim[0].solver.dim_global,
                      sim[0].solver.dim_local,
                      0,
                      sim[0].solver.x,
                      m_snapshots->getColumn(9)->getData(),
                      &(sim[0].solver),
                      &(sim[0].mpi),
                      fname_buffer3);
//      CalculateROMDiff(  &(sim[0].solver),
//                         &(sim[0].mpi) );

      }
    }
  } else {
    printf("ERROR in LSROMObject::train(): m_generator is of size zero!");
  }

  return;
}

/*! compute prediction at given time */
const CAROM::Vector* LSROMObject::predict(const double a_t /*!< time at which to predict solution */ )
{
  int num_rows = m_generator[0]->getSpatialBasis()->numRows();
  CAROM::Vector* recon_init = new CAROM::Vector(num_rows,false);
  if(std::abs(a_t) < 1e-9) {
    printf("a_t %f\n",a_t);
//  projectInitialSolution(*(m_snapshots->getColumn(0)));
//  m_romcoef = new CAROM::Vector(m_projected_init[0]->getData(), m_rdim, false, true);
    m_romcoef =  m_projected_init[0]->getData();
    printf("m_project size %d\n",m_projected_init[0]->dim());
    if (!m_rank) {
      std::cout << "Checking initial condition: ";
      for (int j = 0; j < m_projected_init[0]->dim(); j++) {
            std::cout << (m_projected_init[0]->item(j)) << " ";
      }
      std::cout << std::endl;
    }
//  _ArrayCopy1D_(  m_romcoef,
//                  m_projected_init[0],
//                  m_rdim );
    if (!m_rank) {
      std::cout << "Checking copy initial condition: ";
      for (int j = 0; j < m_rdim; j++) {
            std::cout << m_romcoef[j] << " ";
      }
      std::cout << std::endl;
    }
    recon_init = m_generator[0]->getSpatialBasis()->mult(m_projected_init[0]);
    /* Setup RK44 parameters */
    TimeExplicitRKInitialize();
    printf("checking m_rdim %d\n",m_rdim);
    /* Initialize RK44 working variables */
    TimeInitialize();
  } else {
    TimeRK(a_t);
  }
//    for (int i = 0; i < m_ls.size(); i++) {
//      if (   (a_t >= m_intervals[i].first)
//          && (  (a_t < m_intervals[i].second) || (m_intervals[i].second < 0)  ) ){
//        return m_ls[i]->predict(a_t);
//      }
//    }
//    printf("ERROR in LSROMObject::predict(): m_ls is of size zero or interval not found!");
//    return nullptr;
  return recon_init;
}

void LSROMObject::save(const std::string& a_fname_root /*!< Filename root */) const
{
    return;
}
void LSROMObject::load(const std::string& a_fname_root /*!< Filename root */)
{
    return;
}

void LSROMObject::copyToHyPar(  const CAROM::Vector& a_vec,  /*!< Work vector */
                                    void* a_s, /*!< Array of simulation objects of 
                                                    type #SimulationObject */
                                    int a_idx /*!< Simulation object index */ ) const
{
  SimulationObject* sim = (SimulationObject*) a_s;

  double* u;
  u = sim[a_idx].solver.u;

  std::vector<int> index(sim[a_idx].solver.ndims);

  ArrayCopynD(  sim[a_idx].solver.ndims,
                a_vec.getData(),
                u,
                sim[a_idx].solver.dim_local,
                0,
                sim[a_idx].solver.ghosts,
                index.data(),
                sim[a_idx].solver.nvars );

  return;
}

int LSROMObject::TimeInitialize()
{
  int i;
  /* initialize arrays to NULL, then allocate as necessary */
  m_U.clear();
  m_Udot.clear();

  /* explicit Runge-Kutta methods */
  for (i = 0; i < nstages; i++) {
    m_U.push_back(new CAROM::Vector(m_rdim,false));
    m_Udot.push_back(new CAROM::Vector(m_rdim,false));
  }

//m_U     = (double**) calloc (nstages,sizeof(double*));
//m_Udot  = (double**) calloc (nstages,sizeof(double*));
//for (i = 0; i < nstages; i++) {
//  m_U[i]    = (double*) calloc (m_rdim,sizeof(double));
//  m_Udot[i] = (double*) calloc (m_rdim,sizeof(double));
//}
  return(0);
}

int LSROMObject::TimeExplicitRKInitialize()
{
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
  if (!m_rank) {
    std::cout << "Checking A: ";
    for (int j = 0; j < nstages*nstages; j++) {
         std::cout << A[j] << " ";
    }
    std::cout << std::endl;
    std::cout << "Checking B: ";
    for (int j = 0; j < nstages; j++) {
         std::cout << b[j] << " ";
    }
    std::cout << std::endl;
    std::cout << "Checking C: ";
    for (int j = 0; j < nstages; j++) {
         std::cout << c[j] << " ";
    }
    std::cout << std::endl;
  }
  return(0);
}

int LSROMObject::TimeRK(const double a_t /*!< time at which to predict solution */ )
{
  /* Advance the ROM ODE using RK4 scheme */
  int ns, stage, i;
  printf("checking m_rdim %d\n",m_rdim);
    /* Calculate stage values */
  for (stage = 0; stage < nstages; stage++) {
  
    double stagetime = a_t + c[stage]*m_dt;
    printf("stagetime %f\n",stagetime);
    if (!m_rank) {
      std::cout << "Checking copy initial condition: ";
      for (int j = 0; j < m_rdim; j++) {
            std::cout << m_romcoef[j] << " ";
      }
      std::cout << std::endl;
    }

//  std::cout << "m_romcoef: " << m_romcoef << std::endl;
//  std::cout << "m_U[stage]: " << m_U[stage] << std::endl;
//  std::cout << "m_rdim: " << m_rdim << std::endl;
    _ArrayCopy1D_(  m_romcoef,
                    m_U[stage]->getData(),
                    m_rdim );

    for (i = 0; i < stage; i++) {
      _ArrayAXPY_(  m_Udot[i]->getData(),
                    (m_dt * A[stage*nstages+i]),
                    m_U[stage]->getData(),
                    m_rdim );
    }

    m_Udot[stage] = m_romhyperb->mult(m_U[stage]);
    if (!m_rank) {
      std::cout << "Checking copy initial condition: ";
      for (int j = 0; j < m_rdim; j++) {
            std::cout << m_Udot[stage]->item(j) << " ";
      }
      std::cout << std::endl;
    }
//  TS->RHSFunction( (m_Udot[stage]),
//                   (m_U[stage]),
//                   &(sim[ns].solver),
//                   &(sim[ns].mpi),
//                   stagetime);

//  for (ns = 0; ns < nsims; ns++) {
//    if (sim[ns].solver.PostStage) {
//      sim[ns].solver.PostStage(  (TS->U[stage] + TS->u_offsets[ns]),
//                                 &(sim[ns].solver),
//                                 &(sim[ns].mpi),
//                                 stagetime); CHECKERR(ierr);
//    }
//  }
//
//}
//
///* Step completion */
//for (stage = 0; stage < params->nstages; stage++) {
//
//  for (ns = 0; ns < nsims; ns++) {
//    _ArrayAXPY_(  (TS->Udot[stage] + TS->u_offsets[ns]),
//                  (TS->dt * params->b[stage]),
//                  (sim[ns].solver.u),
//                  (TS->u_sizes[ns]) );
//  }
//
  }
  exit (0);
  return(0);
}

#endif

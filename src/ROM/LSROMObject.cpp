#ifdef with_librom

/*! @file LSROMObject.cpp
 *  @brief Member functions of the class #LSROMObject
 *  @author Ping-Hsuan Tsai
*/

#include <basic.h>
#include <rom_object_ls.h>
#include <simulation_object.h>
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
  m_spatialbasis.clear();

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

/*! take a sample (solution snapshot) */
void LSROMObject::takeSample(  const CAROM::Vector& a_U, /*!< solution vector */
                                const double a_time, /*!< sample time */ 
                                void* a_s )
{
  SimulationObject* sim = (SimulationObject*) a_s;
  if (m_tic == 0) {

    m_options.push_back(new CAROM::Options(m_vec_size, max_num_snapshots, 1, update_right_SV));
    m_generator.push_back(new CAROM::BasisGenerator(*m_options[m_curr_win], isIncremental, basisName));
    m_intervals.push_back( Interval(a_time, m_t_final) );
    m_ls_is_trained.push_back(false);

    if (!m_rank) {
      printf( "LSROMObject::takeSample() - creating new generator object for sim. domain %d, var %d, t=%f (total: %d).\n",
              m_sim_idx, m_var_idx, m_intervals[m_curr_win].first, m_generator.size());
    }
    /* QUESTION: should a_U be centered? */
    bool addSample = m_generator[m_curr_win]->takeSample( a_U.getData(), a_time, m_dt );
    /* Below is printing the singular value */
//  m_S = m_generator[m_curr_win]->getSingularValues();
//  std::cout << "m_S: ";
//  for (int i = 0; i < m_S->dim(); i++) {
//          std::cout << (m_S->item(i)) << " ";
//  }
//  std::cout << std::endl;
//  m_generator[m_curr_win]->writeSnapshot();
//    CAROM::Vector a_tmp(m_generator[0]->getSnapshotMatrix()->getColumn(0),sim[0].solver.dim_global,false);
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
          printf( "inside m_write_snapshot_mat \n");
          char idx_string[_MAX_STRING_SIZE_];
          sprintf(idx_string, "%04d", i);
          std::string fname_root(m_dirname + "/snapshot_mat_"+std::string(idx_string));
          // Need additional process in order to visualize it
          m_generator[i]->getSnapshotMatrix()->write(fname_root);
        }
//      m_generator[i]->train(m_rdim);
//      libROM merge phase code here
        printf( "before g\n");
        const int rom_dim = m_generator[i]->getSpatialBasis()->numColumns();
        printf( "numcolumns %d\n",rom_dim);
        const CAROM::Vector* sing_vals = m_generator[i]->getSingularValues();
        std::cout << "m_S: ";
        for (int i = 0; i < sing_vals->dim(); i++) {
                std::cout << (sing_vals->item(i)) << " ";
        }
        std::cout << std::endl;
//      m_generator[i]->endSamples();
//      const Matrix* mat;
        int num_rows = m_generator[i]->getSpatialBasis()->numRows();
        int num_cols = m_generator[i]->getSpatialBasis()->numColumns();
//      const int d_dim = m_generator[i]->getSpatialBasis()->numRows();
//      CAROM::Matrix* d_basis;
//      d_basis = new CAROM::Matrix(num_rows, num_cols, false);
//      printf( "rows, cols: %d %d \n",num_rows,num_cols);

        for (int j = 0; j < num_cols; j++){
          double* vec_data = m_generator[0]->getSpatialBasis()->getColumn(j)->getData();
          const char* filename = "basis_";
          std::string file_name = std::string(filename) + std::to_string(j);
          char* file_name_buffer = new char[file_name.length() + 1];
          std::strcpy(file_name_buffer, file_name.c_str());

          WriteArray( sim[0].solver.ndims,
                      sim[0].solver.nvars,
                      sim[0].solver.dim_global,
                      sim[0].solver.dim_local,
                      0,
                      sim[0].solver.x,
                      vec_data,
                      &(sim[0].solver),
                      &(sim[0].mpi),
                      file_name_buffer);
        }
      }
    }
  } else {
    printf("ERROR in LSROMObject::train(): m_generator is of size zero!");
  }

  return;
}

void LSROMObject::save(const std::string& a_fname_root /*!< Filename root */) const
{
    return;
}
void LSROMObject::load(const std::string& a_fname_root /*!< Filename root */)
{
    return;
}

#endif

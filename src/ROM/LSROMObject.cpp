#ifdef with_librom

/*! @file LSROMObject.cpp
 *  @brief Member functions of the class #LSROMObject
 *  @author Ping-Hsuan Tsai
*/

#include <basic.h>
#include <rom_object_ls.h>

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
                                const double a_time /*!< sample time */ )
{
  if (m_tic == 0) {

    m_options.push_back(new CAROM::Options(m_vec_size, max_num_snapshots, 1, update_right_SV));
    m_generator.push_back(new CAROM::BasisGenerator(*m_options[m_curr_win], isIncremental, basisName));
    m_intervals.push_back( Interval(a_time, m_t_final) );
    m_ls_is_trained.push_back(false);

    if (!m_rank) {
      printf( "LSROMObject::takeSample() - creating new generator object for sim. domain %d, var %d, t=%f (total: %d).\n",
              m_sim_idx, m_var_idx, m_intervals[m_curr_win].first, m_generator.size());
    }
//  bool addSample = m_generator[m_curr_win]->takeSample( a_U.getData(), a_time, m_dt );

//} else {

//  bool addSample = m_generator[m_curr_win]->takeSample( a_U.getData(), a_time, m_dt );

//  if (m_tic%m_num_window_samples == 0) {
//    printf("checking m_tic %d \n",m_tic);

//    m_intervals[m_curr_win].second = a_time;
//    int ncol = m_ls[m_curr_win]->getSnapshotMatrix()->numColumns();
//    if (!m_rank) {
//      printf( "LSROMObject::train() - training LS object %d for sim. domain %d, var %d with %d samples.\n", 
//              m_curr_win, m_sim_idx, m_var_idx, ncol );
//    }

//    if (m_write_snapshot_mat) {
//      char idx_string[_MAX_STRING_SIZE_];
//      sprintf(idx_string, "%04d", m_curr_win);
//      std::string fname_root(m_dirname + "/snapshot_mat_"+std::string(idx_string));
//      m_ls[m_curr_win]->getSnapshotMatrix()->write(fname_root);
//    }
//    m_ls[m_curr_win]->train(m_rdim);
//    m_ls_is_trained[m_curr_win] = true;

//    m_curr_win++;

//    m_ls.push_back( new CAROM::LS(m_vec_size, m_dt) );
//    m_ls_is_trained.push_back(false);
//    m_intervals.push_back( Interval(a_time, m_t_final) );
//    m_ls[m_curr_win]->takeSample( a_U.getData(), a_time );
//    if (!m_rank) {
//      printf("LSROMObject::takeSample() - creating new LS object for sim. domain %d, var %d, t=%f (total: %d).\n",
//             m_sim_idx, m_var_idx, m_intervals[m_curr_win].first, m_ls.size());
//    }
//  }

  }

  m_tic++;
  return;
}

/*! train the LS objects */
void LSROMObject::train()
{
  /* make sure the number of columns for the last LS isn't less than m_rdim */
  {
    int last_win = m_generator.size() - 1;
    int num_columns( m_generator[last_win]->getSnapshotMatrix()->numColumns() );
    if (num_columns <= m_rdim) {
      m_generator.pop_back();
      m_intervals.pop_back();
      m_intervals[m_intervals.size()-1].second = m_t_final;
      if (!m_rank) {
        printf("LSROMObject::train() - last window LS for sim. domain %d, var %d has %d sample(s) only; deleted it ",
               m_sim_idx,
               m_var_idx,
               num_columns );
        printf("(total: %d).\n", m_generator.size());
      }
    }
  }

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
        }
//      m_generator[i]->train(m_rdim);
//      libROM merge phase code here
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

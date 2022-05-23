#ifdef with_librom

/*! @file libROMInterface.cpp
 *  @brief Functions that implement interface with libROM
 *  @author Debojyoti Ghosh
*/

#include <string.h>
#include <arrayfunctions.h>
#include <simulation_object.h>
#include <librom_interface.h>

/*! Constructor 
    This function will also look into the file
    \b librom.inp for DMD-specific options. 
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
    dmd_num_win_samples    | int          | #DMDROMObject::m_num_window_samples           | INT_MAX
    dmd_dirname            | string       | #DMDROMObject::m_dirname                      | "DMD"
    dmd_write_snapshot_mat | bool         | #DMDROMObject::m_write_snapshot_mat           | false

    Note: other keywords in this file may be read by other functions.
   
*/
DMDROMObject::DMDROMObject( const int     a_vec_size, /*!< vector size */
                            const double  a_dt,       /*!< time step size */
                            const int     a_rdim,     /*!< latent space dimension */
                            const int     a_rank,     /*!< MPI rank of this process */
                            const int     a_nproc     /*!< Number of MPI processes */ )
{
  m_rank = a_rank;
  m_nproc = a_nproc;

  m_dmd.clear();
  m_dmd_is_trained.clear();
  m_intervals.clear();
  m_vec_size = a_vec_size;
  m_dt = a_dt;
  m_t_final = -1;
  m_rdim = a_rdim;

  m_num_window_samples = INT_MAX;

  char dirname_c_str[_MAX_STRING_SIZE_] = "DMD";
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
          if (std::string(word) == "dmd_num_win_samples") {
            ferr = fscanf(in,"%d", &m_num_window_samples); if (ferr != 1) return;
          } else if (std::string(word) == "dmd_dirname") {
            ferr = fscanf(in,"%s", dirname_c_str); if (ferr != 1) return;
          } else if (std::string(word) == "dmd_write_snapshot_mat") {
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
    printf("libROM DMD inputs:\n");
    printf("  number of samples per window:   %d\n", m_num_window_samples);
    printf("  directory name for DMD onjects: %s\n", dirname_c_str);
    printf("  write snapshot matrix to file:  %s\n", write_snapshot_mat);
  }

#ifndef serial
  MPI_Bcast(&m_num_window_samples,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(dirname_c_str,_MAX_STRING_SIZE_,MPI_CHAR,0,MPI_COMM_WORLD);
  MPI_Bcast(write_snapshot_mat,_MAX_STRING_SIZE_,MPI_CHAR,0,MPI_COMM_WORLD);
#endif

  m_dirname = std::string( dirname_c_str );
  m_write_snapshot_mat = (std::string(write_snapshot_mat) == "true");

  if (m_num_window_samples <= m_rdim) {
    printf("ERROR:DMDROMObject::DMDROMObject() - m_num_window_samples <= m_rdim!!");
  }

  m_tic = 0;
  m_curr_win = 0;
}

/*! take a sample (solution snapshot) */
void DMDROMObject::takeSample(  const CAROM::Vector& a_U, /*!< solution vector */
                                const double a_time /*!< sample time */ )
{
  if (m_tic == 0) {

    m_dmd.push_back( new CAROM::DMD(m_vec_size, m_dt) );
    m_dmd_is_trained.push_back(false);
    m_intervals.push_back( Interval(a_time, m_t_final) );

    if (!m_rank) {
      printf("DMDROMObject::takeSample() - creating new DMD object, t=%f (total: %d).\n",
             m_intervals[m_curr_win].first, m_dmd.size());
    }
    m_dmd[m_curr_win]->takeSample( a_U.getData(), a_time );

  } else {

    m_dmd[m_curr_win]->takeSample( a_U.getData(), a_time );

    if (m_tic%m_num_window_samples == 0) {

      m_intervals[m_curr_win].second = a_time;
      int ncol = m_dmd[m_curr_win]->getSnapshotMatrix()->numColumns();
      if (!m_rank) {
        printf("DMDROMObject::train() - training DMD object %d with %d samples.\n", 
                m_curr_win, ncol );
      }

      if (m_write_snapshot_mat) {
        char idx_string[_MAX_STRING_SIZE_];
        sprintf(idx_string, "%04d", m_curr_win);
        std::string fname_root(m_dirname + "/snapshot_mat_"+std::string(idx_string));
        m_dmd[m_curr_win]->getSnapshotMatrix()->write(fname_root);
      }
      m_dmd[m_curr_win]->train(m_rdim);
      m_dmd_is_trained[m_curr_win] = true;

      m_curr_win++;

      m_dmd.push_back( new CAROM::DMD(m_vec_size, m_dt) );
      m_dmd_is_trained.push_back(false);
      m_intervals.push_back( Interval(a_time, m_t_final) );
      m_dmd[m_curr_win]->takeSample( a_U.getData(), a_time );
      if (!m_rank) {
        printf("DMDROMObject::takeSample() - creating new DMD object, t=%f (total: %d).\n",
               m_intervals[m_curr_win].first, m_dmd.size());
      }
    }

  }

  m_tic++;
  return;
}

/*! train the DMD objects */
void DMDROMObject::train()
{
  /* make sure the number of columns for the last DMD isn't less than m_rdim */
  {
    int last_win = m_dmd.size() - 1;
    int num_columns( m_dmd[last_win]->getSnapshotMatrix()->numColumns() );
    if (num_columns <= m_rdim) {
      m_dmd.pop_back();
      m_intervals.pop_back();
      m_intervals[m_intervals.size()-1].second = m_t_final;
      if (!m_rank) {
        printf("DMDROMObject::train() - last window DMD has %d sample(s) only; deleted it ",
               num_columns );
        printf("(total: %d).\n", m_dmd.size());
      }
    }
  }

  if (m_dmd.size() > 0) {
    for (int i = 0; i < m_dmd.size(); i++) {
      if (!m_dmd_is_trained[i]) {
        int ncol = m_dmd[i]->getSnapshotMatrix()->numColumns();
        if (!m_rank) {
          printf("DMDROMObject::train() - training DMD object %d with %d samples.\n", i, ncol );
        }
        if (m_write_snapshot_mat) {
          char idx_string[_MAX_STRING_SIZE_];
          sprintf(idx_string, "%04d", i);
          std::string fname_root(m_dirname + "/snapshot_mat_"+std::string(idx_string));
          m_dmd[i]->getSnapshotMatrix()->write(fname_root);
        }
        m_dmd[i]->train(m_rdim);
      }
    }
  } else {
    printf("ERROR in DMDROMObject::train(): m_dmd is of size zero!");
  }

  return;
}

/*! Save DMD objects to file: the DMD object files will be saved in the subdirectory
 *  with the name #DMDROMObject::m_dirname. They are in a format that libROM can read
 *  from.
 *
 *  \b Note: If the subdirectory #DMDROMObject::m_dirname does not exist, the DMD
 *  objects will not be written (even though the screen output will claim they were
 *  written)!. The code may not report any error, or one may see HDF5 file writing
 *  errors.
*/
void DMDROMObject::save(const std::string& a_fname_root /*!< Filename root */) const
{
  std::string fname_root = m_dirname + "/";
  std::string summary_fname_root = m_dirname + "/";
  std::string header_fname = m_dirname + "/";
  if (a_fname_root == "") {
    fname_root += "dmdobj_";
    summary_fname_root += "dmd_summary_";
    header_fname += "dmd_header.dat";
  } else {
    fname_root += (a_fname_root+"_dmdobj_");
    summary_fname_root += (a_fname_root+"_dmd_summary_");
    header_fname += (a_fname_root+"_dmd_header.dat");
  }

  if (!m_rank) {
    FILE* out;
    out = fopen(header_fname.c_str(), "w");
    fprintf(out, "%d\n", m_dmd.size());
    for (int i = 0; i < m_dmd.size(); i++) {
      fprintf(out, "%1.16e %1.16e\n", m_intervals[i].first, m_intervals[i].second);
    }
    fclose(out);
  }

  for (int i = 0; i < m_dmd.size(); i++) {
    char idx_string[_MAX_STRING_SIZE_];
    sprintf(idx_string, "%04d", i);
    std::string fname = fname_root + std::string(idx_string);
    std::string summary_fname = summary_fname_root + std::string(idx_string);
    if (!m_rank) {
      printf( "  Saving DMD object and summary (%s, %s).\n", 
              fname.c_str(), summary_fname.c_str() );
    }
    m_dmd[i]->save(fname);
    //m_dmd[i]->summary(summary_fname);
  }

  return;
}

/*! Load DMD objects from file: the DMD object files must be in the subdirectory
 *  with the name #DMDROMObject::m_dirname. They must be in a format that libROM can read
 *  from. The number of objects and their parameters must correspond to the simulation
 *  being run for which this object is defined.
*/
void DMDROMObject::load(const std::string& a_fname_root /*!< Filename root */)
{
  std::string fname_root = m_dirname + "/";
  std::string header_fname = m_dirname + "/";
  if (a_fname_root == "") {
    fname_root += "dmdobj_";
    header_fname += "dmd_header.dat";
  } else {
    fname_root += (a_fname_root+"_dmdobj_");
    header_fname += (a_fname_root+"_dmd_header.dat");
  }

  int num_dmds = 0;
  std::vector<double> intervals_start(0), intervals_end(0);
  if (!m_rank) {
    FILE* in;
    in = fopen(header_fname.c_str(), "r");
    fscanf(in, "%d\n", &num_dmds);
    intervals_start.resize(num_dmds);
    intervals_end.resize(num_dmds);
    for (int i = 0; i < num_dmds; i++) {
      fscanf(in, "%lf", &intervals_start[i]);
      fscanf(in, "%lf", &intervals_end[i]);
    }
    fclose(in);
  }
#ifndef serial
  MPI_Bcast(&num_dmds,1,MPI_INT,0,MPI_COMM_WORLD);
  if (m_rank) {
    intervals_start.resize(num_dmds);
    intervals_end.resize(num_dmds);
  }
  MPI_Bcast(intervals_start.data(),num_dmds,MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Bcast(intervals_end.data(),num_dmds,MPI_DOUBLE,0,MPI_COMM_WORLD);
#endif

  for (int i = 0; i < num_dmds; i++) {

    m_intervals.push_back( Interval(intervals_start[i],intervals_end[i]) );

    char idx_string[_MAX_STRING_SIZE_];
    sprintf(idx_string, "%04d", i);
    std::string fname = fname_root + std::string(idx_string);
    if (!m_rank) {
      printf( "  Loading DMD object (%s), time window=[%1.2e,%1.2e].\n", 
              fname.c_str(),
              m_intervals[i].first, m_intervals[i].second );
    }
    m_dmd.push_back( new CAROM::DMD(fname) );
  }

  return;
}

/*! Define the libROM interface
  
    This function also reads libROM-related inputs from the file
    \b librom.inp. Rank 0 reads in the inputs and broadcasts
    them to all the processors.\n\n
    The format of \b librom.inp is as follows:\n

        begin
            <keyword>   <value>
            ...
            <keyword>   <value>
        end

    where the list of keywords and their type that this function will read are:\n
    Keyword name       | Type         | Variable                                      | Default value
    ------------------ | ------------ | --------------------------------------------- | -------------------
    rdim               | int          | #libROMInterface::m_rdim                      | -1
    sampling_frequency | int          | #libROMInterface::m_sampling_freq             | 1
    mode               | string       | #libROMInterface::m_mode                      | "train"
    type               | string       | #libROMInterface::m_rom_type                  | "DMD"
    save_to_file       | string       | #libROMInterface::m_save_ROM                  | "true"

    Note: other keywords in this file may be read by other functions. The default value for \a rdim is invalid, so it \b must be specified.
   
*/
void libROMInterface::define( void*   a_s, /*!< Array of simulation objects of type #SimulationObject */
                              int     a_nsims, /*!< number of simulation objects */
                              int     a_rank,  /*!< MPI rank of this process */
                              int     a_nproc, /*!< Number of MPI processes */
                              double  a_dt     /*!< Time step size */ )
{
  const SimulationObject* sim = (const SimulationObject*) a_s;

  m_rank = a_rank;
  m_nproc = a_nproc;
  m_nsims = a_nsims;

  m_vec_size.resize(m_nsims);
  for (int ns = 0; ns < m_nsims; ns++) {
    m_vec_size[ns] = (sim[ns].solver.npoints_local * sim[ns].solver.nvars);
  }

  char mode_c_str[_MAX_STRING_SIZE_];
  char type_c_str[_MAX_STRING_SIZE_];
  char save_c_str[_MAX_STRING_SIZE_];

  if (!m_rank) {

    FILE *in;
    in = fopen(_LIBROM_INP_FNAME_,"r");

    if (!in) {

      strcpy( mode_c_str, _ROM_MODE_NONE_ );
      strcpy( type_c_str, _ROM_TYPE_NONE_ );
      strcpy( save_c_str, "false" );

    } else {

      strcpy( mode_c_str, _ROM_MODE_TRAIN_ );
      strcpy( type_c_str,_ROM_TYPE_DMD_ );
      strcpy( save_c_str, "true" );

      int ferr;
      char word[_MAX_STRING_SIZE_];
      ferr = fscanf(in,"%s", word); if (ferr != 1) return;

      if (std::string(word) == "begin") {
        while (std::string(word) != "end") {
  	      ferr = fscanf(in,"%s",word); if (ferr != 1) return;
          if (std::string(word) == "rdim") {
            ferr = fscanf(in,"%d", &m_rdim); if (ferr != 1) return;
          } else if (std::string(word) == "sampling_frequency") {
            ferr = fscanf(in,"%d", &m_sampling_freq); if (ferr != 1) return;
          } else if (std::string(word) == "mode") {
            ferr = fscanf(in,"%s", mode_c_str); if (ferr != 1) return;
          } else if (std::string(word) == "type") {
            ferr = fscanf(in,"%s", type_c_str); if (ferr != 1) return;
          } else if (std::string(word) == "save_to_file") {
            ferr = fscanf(in,"%s", save_c_str); if (ferr != 1) return;
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
    printf("libROM inputs and parameters:\n");
    printf("  reduced model dimensionality:  %d\n", m_rdim);
    printf("  sampling frequency:  %d\n", m_sampling_freq);
    printf("  mode: %s\n", mode_c_str);
    printf("  type: %s\n", type_c_str);
    printf("  save to file: %s\n", save_c_str);
  }

#ifndef serial
  MPI_Bcast(&m_rdim,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&m_sampling_freq,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(mode_c_str,_MAX_STRING_SIZE_,MPI_CHAR,0,MPI_COMM_WORLD);
  MPI_Bcast(type_c_str,_MAX_STRING_SIZE_,MPI_CHAR,0,MPI_COMM_WORLD);
  MPI_Bcast(save_c_str,_MAX_STRING_SIZE_,MPI_CHAR,0,MPI_COMM_WORLD);
#endif

  m_mode = std::string( mode_c_str );
  m_rom_type = std::string( type_c_str );

  m_save_ROM = (((std::string(save_c_str) == "true")) && (m_mode == _ROM_MODE_TRAIN_));

  m_rom.resize( m_nsims, nullptr );
  m_U.resize( m_nsims, nullptr );
  if (m_rom_type == _ROM_TYPE_DMD_) {
    for (int ns = 0; ns < m_nsims; ns++) {
      m_rom[ns] = new DMDROMObject( m_vec_size[ns], 
                                    m_sampling_freq*a_dt, 
                                    m_rdim, 
                                    m_rank, 
                                    m_nproc );
    }
  }
  for (int ns = 0; ns < m_nsims; ns++) {
    m_U[ns] = new CAROM::Vector(m_vec_size[ns],true);
  }

  m_is_defined = true;
  return;
}

/*! Copy HyPar solution to the work vectors a_U. Note that the HyPar solution has ghost
 * points but the work vectors a_U is a serialized vector of the solution without ghost
 * points */
void libROMInterface::copyFromHyPar(  std::vector<CAROM::Vector*>& a_U,  /*!< work vector */
                                      void* a_s   /*!< Array of simulation objects of 
                                                       type #SimulationObject */ )
{
  const SimulationObject* sim = (const SimulationObject*) a_s;

  for (int ns = 0; ns < m_nsims; ns++) {
    double* vec = a_U[ns]->getData();
    const double* u = sim[ns].solver.u;

    std::vector<int> index(sim[ns].solver.ndims);

    ArrayCopynD(  sim[ns].solver.ndims,
                  u,
                  vec,
                  sim[ns].solver.dim_local,
                  sim[ns].solver.ghosts,
                  0,
                  index.data(),
                  sim[ns].solver.nvars );
  }

  return;
}

/*! Copy the work vectors a_U to HyPar solution. Note that the HyPar solution has ghost
 * points but the work vectors a_U is a serialized vector of the solution without ghost
 * points.
 *
 * \b Note: if libROMInterface::m_mode is "predict", then the vector will be copied
 * into the main solution variable HyPar::u; otherwise it will be copied into the
 * variable for ROM prediction HyPar::u_rom_predicted */
void libROMInterface::copyToHyPar(  const std::vector<CAROM::Vector*>& a_U,  /*!< Work vector */
                                    void* a_s   /*!< Array of simulation objects of 
                                                     type #SimulationObject */ ) const
{
  SimulationObject* sim = (SimulationObject*) a_s;

  for (int ns = 0; ns < m_nsims; ns++) {
    const double* vec = a_U[ns]->getData();
    double* u;
    if (m_mode == _ROM_MODE_PREDICT_) {
      u = sim[ns].solver.u;
    } else {
      u = sim[ns].solver.u_rom_predicted;
    }

    std::vector<int> index(sim[ns].solver.ndims);

    ArrayCopynD(  sim[ns].solver.ndims,
                  vec,
                  u,
                  sim[ns].solver.dim_local,
                  0,
                  sim[ns].solver.ghosts,
                  index.data(),
                  sim[ns].solver.nvars );
  }

  return;
}

/*! Copy a vector a_vec to HyPar solution. Note that the HyPar solution has ghost
 * points but the vector a_vec is a serialized vector of the solution without ghost
 * points.
 *
 * \b Note: if libROMInterface::m_mode is "predict", then the vector will be copied
 * into the main solution variable HyPar::u; otherwise it will be copied into the
 * variable for ROM prediction HyPar::u_rom_predicted */
void libROMInterface::copyToHyPar(  const CAROM::Vector& a_vec,  /*!< Work vector */
                                    void* a_s, /*!< Array of simulation objects of 
                                                    type #SimulationObject */
                                    int a_idx /*!< Simulation object index */ ) const
{
  SimulationObject* sim = (SimulationObject*) a_s;

  double* u;
  if (m_mode == _ROM_MODE_PREDICT_) {
    u = sim[a_idx].solver.u;
  } else {
    u = sim[a_idx].solver.u_rom_predicted;
  }
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

#endif

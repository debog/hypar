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
    Keyword name       | Type         | Variable                                      | Default value
    ------------------ | ------------ | --------------------------------------------- | -------------------
    dmd_num_win_samples| int          | #DMDROMObject::m_num_window_samples           | INT_MAX

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
  m_intervals.clear();
  m_vec_size = a_vec_size;
  m_dt = a_dt;
  m_rdim = a_rdim;

  m_num_window_samples = INT_MAX;

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
    printf("  number of samples per window:  %d\n", m_num_window_samples);
  }

#ifndef serial
  MPI_Bcast(&m_num_window_samples,1,MPI_INT,0,MPI_COMM_WORLD);
#endif

  if (m_num_window_samples <= m_rdim) {
    printf("ERROR:DMDROMObject::DMDROMObject() - m_num_window_samples <= m_rdim!!");
  }

  m_tic = 0;
  m_curr_win = 0;
}

/*! train the DMD objects */
void DMDROMObject::train()
{
  /* make sure last DMD doesn't have only one column in snapshot matrix;
   * this may happen if total number of iterations is an integer multiple
   * of m_num_window_samples */
  {
    int last_win = m_dmd.size() - 1;
    int num_columns( m_dmd[last_win]->getSnapshotMatrix()->numColumns() );
    if (num_columns <= 1) {
      m_dmd.pop_back();
      m_intervals.pop_back();
      m_intervals[m_intervals.size()-1].second = DBL_MAX;
      if (!m_rank) {
        printf("DMDROMObject::train() - last window DMD has %d sample(s) only; deleted it ",
               num_columns );
        printf("(total: %d).\n", m_dmd.size());
      }
    }
  }

  if (m_dmd.size() > 0) {
    for (int i = 0; i < m_dmd.size(); i++) {
      int ncol = m_dmd[i]->getSnapshotMatrix()->numColumns();
      if (!m_rank) {
        printf("DMDRomObject::train() - training DMD object %d with %d samples.\n", i, ncol );
      }
      m_dmd[i]->train(m_rdim);
    }
  } else {
    printf("ERROR in DMDROMObject::train(): m_dmd is of size zero!");
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
    rdim               | int          | #libROMInterface::m_rdim                      | 10
    sampling_frequency | int          | #libROMInterface::m_sampling_freq             | 1
    mode               | string       | #libROMInterface::m_mode                      | "train"
    type               | string       | #libROMInterface::m_rom_type                  | "DMD"

    Note: other keywords in this file may be read by other functions.
   
*/
void libROMInterface::define( void*   a_s, /*!< Array of simulation objects of type #SimulationObject */
                              int     a_nsims, /*!< number of simulation objects */
                              int     a_rank,  /*!< MPI rank of this process */
                              int     a_nproc, /*!< Number of MPI processes */
                              double  a_dt     /*!< Time step size */ )
{
  SimulationObject* sim = (SimulationObject*) a_s;

  m_rank = a_rank;
  m_nproc = a_nproc;
  m_nsims = a_nsims;

  m_vec_size = 0;
  m_vec_offsets.resize(m_nsims);
  for (int ns = 0; ns < m_nsims; ns++) {
    m_vec_offsets[ns] = m_vec_size;
    m_vec_size += (sim[ns].solver.npoints_local * sim[ns].solver.nvars);
  }

  char mode_c_str[_MAX_STRING_SIZE_];
  char type_c_str[_MAX_STRING_SIZE_];

  if (!m_rank) {

    FILE *in;
    in = fopen(_LIBROM_INP_FNAME_,"r");

    if (!in) {

      strcpy( mode_c_str, "none" );
      strcpy( type_c_str, "none" );

    } else {

      strcpy( mode_c_str, "train" );
      strcpy( type_c_str, "DMD" );

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
    printf("  local vector size: %d\n", m_vec_size);
  }

#ifndef serial
  MPI_Bcast(&m_rdim,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&m_sampling_freq,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(mode_c_str,_MAX_STRING_SIZE_,MPI_CHAR,0,MPI_COMM_WORLD);
  MPI_Bcast(type_c_str,_MAX_STRING_SIZE_,MPI_CHAR,0,MPI_COMM_WORLD);
#endif

  m_mode = std::string( mode_c_str );
  m_rom_type = std::string( type_c_str );

  if (m_mode == "train") {
    if (m_rom_type == _ROM_TYPE_DMD_) {
      m_rom = new DMDROMObject( m_vec_size, a_dt, m_rdim, m_rank, m_nproc );
    }
    m_U = new double[m_vec_size];
  }

  m_is_defined = true;
  return;
}

/*! Copy HyPar solution to the work vector a_U. Note that the HyPar solution has ghost
 * points but the work vector a_U is a serialized vector of the solution without ghost
 * points */
void libROMInterface::copyFromHyPar(  double* a_U,  /*!< work vector */
                                      void* a_s     /*!< Array of simulation objects of 
                                                         type #SimulationObject */ )
{
  SimulationObject* sim = (SimulationObject*) a_s;

  for (int ns = 0; ns < m_nsims; ns++) {
    double* vec = a_U + m_vec_offsets[ns];
    double* u = sim[ns].solver.u;

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

/*! Copy the work vector a_U to HyPar solution. Note that the HyPar solution has ghost
 * points but the work vector a_U is a serialized vector of the solution without ghost
 * points */
void libROMInterface::copyToHyPar(  double* a_U,  /*!< Work vector */
                                    void* a_s     /*!< Array of simulation objects of 
                                                       type #SimulationObject */ )
{
  SimulationObject* sim = (SimulationObject*) a_s;

  for (int ns = 0; ns < m_nsims; ns++) {
    double* vec = a_U + m_vec_offsets[ns];
    double* u = sim[ns].solver.u_rom_predicted;

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

#endif

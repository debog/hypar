#ifdef with_librom

/*! @file libROMInterface.cpp
 *  @brief Functions that implement interface with libROM
 *  @author Debojyoti Ghosh
*/

#include <string.h>
#include <arrayfunctions.h>
#include <simulation_object.h>
#include <rom_object_dmd.h>
#include <librom_interface.h>

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
    component_mode     | string       | #libROMInterface::m_comp_mode                 | "monolithic"
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

  char mode_c_str[_MAX_STRING_SIZE_];
  char comp_mode_c_str[_MAX_STRING_SIZE_];
  char type_c_str[_MAX_STRING_SIZE_];
  char save_c_str[_MAX_STRING_SIZE_];

  if (!m_rank) {

    FILE *in;
    in = fopen(_LIBROM_INP_FNAME_,"r");

    if (!in) {

      strcpy( mode_c_str, _ROM_MODE_NONE_ );
      strcpy( comp_mode_c_str, _ROM_COMP_MODE_MONOLITHIC_ );
      strcpy( type_c_str, _ROM_TYPE_NONE_ );
      strcpy( save_c_str, "false" );

      m_rdim = -1;
      m_sampling_freq = 1;

    } else {

      strcpy( mode_c_str, _ROM_MODE_TRAIN_ );
      strcpy( comp_mode_c_str, _ROM_COMP_MODE_MONOLITHIC_ );
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
          } else if (std::string(word) == "component_mode") {
            ferr = fscanf(in,"%s", comp_mode_c_str); if (ferr != 1) return;
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
    if (std::string(mode_c_str) != _ROM_MODE_NONE_) {
      printf("libROMInterface inputs and parameters:\n");
      printf("  reduced model dimensionality:  %d\n", m_rdim);
      printf("  sampling frequency:  %d\n", m_sampling_freq);
      printf("  mode: %s\n", mode_c_str);
      printf("  component mode: %s\n", comp_mode_c_str);
      printf("  type: %s\n", type_c_str);
      printf("  save to file: %s\n", save_c_str);
    }
  }

#ifndef serial
  MPI_Bcast(&m_rdim,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&m_sampling_freq,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(mode_c_str,_MAX_STRING_SIZE_,MPI_CHAR,0,MPI_COMM_WORLD);
  MPI_Bcast(comp_mode_c_str,_MAX_STRING_SIZE_,MPI_CHAR,0,MPI_COMM_WORLD);
  MPI_Bcast(type_c_str,_MAX_STRING_SIZE_,MPI_CHAR,0,MPI_COMM_WORLD);
  MPI_Bcast(save_c_str,_MAX_STRING_SIZE_,MPI_CHAR,0,MPI_COMM_WORLD);
#endif

  m_mode = std::string( mode_c_str );
  m_comp_mode = std::string( comp_mode_c_str );
  m_rom_type = std::string( type_c_str );

  m_save_ROM = (((std::string(save_c_str) == "true")) && (m_mode == _ROM_MODE_TRAIN_));

  if (m_mode != _ROM_MODE_NONE_) {

    m_vec_size.resize(m_nsims);
    for (int ns = 0; ns < m_nsims; ns++) {
      m_vec_size[ns] = (sim[ns].solver.npoints_local);
      if (m_comp_mode == _ROM_COMP_MODE_MONOLITHIC_) {
        m_vec_size[ns] *= (sim[ns].solver.nvars);
      }
    }

    m_rom.clear();
    m_U.clear();
    if (m_comp_mode == _ROM_COMP_MODE_MONOLITHIC_) {
      for (int ns = 0; ns < m_nsims; ns++) {
        if (m_rom_type == _ROM_TYPE_DMD_) {
          m_rom.push_back(new DMDROMObject( m_vec_size[ns],
                                            m_sampling_freq*a_dt,
                                            m_rdim,
                                            m_rank,
                                            m_nproc,
                                            ns ) );
        }
        m_U.push_back(new CAROM::Vector(m_vec_size[ns],true));
      }
    } else if (m_comp_mode == _ROM_COMP_MODE_COMPONENTWISE_) {
      for (int ns = 0; ns < m_nsims; ns++) {
        m_ncomps.push_back(sim[ns].solver.nvars);
        for (int v = 0; v < sim[ns].solver.nvars; v++) {
          if (m_rom_type == _ROM_TYPE_DMD_) {
            m_rom.push_back(new DMDROMObject( m_vec_size[ns],
                                              m_sampling_freq*a_dt,
                                              m_rdim,
                                              m_rank,
                                              m_nproc,
                                              ns,
                                              v ) );
          }
          m_U.push_back(new CAROM::Vector(m_vec_size[ns],true));
        }
      }
    }

  } else {

    m_rom.clear();
    m_U.clear();
    m_ncomps.clear();

  }

  m_train_wctime = 0;
  m_predict_wctime = 0;

  m_is_defined = true;
  return;
}

/*! Take a sample for training */
void libROMInterface::takeSample( void* a_s,        /*!< Array of simulation objects of
                                                         type #SimulationObject */
                                  const double a_t  /*!< Current simulation time */ )
{
  if (m_U.size() != m_rom.size()) {
    printf( "ERROR in libROMInterface::takeSample(): m_U.size != m_rom.size() on rank %d!!\n",
            m_rank );
  }

  copyFromHyPar( m_U, a_s );
  gettimeofday(&m_train_start, NULL);
  for (int i = 0; i < m_rom.size(); i++) {
    m_rom[i]->takeSample( *(m_U[i]), a_t );
  }
  gettimeofday(&m_train_end, NULL);

  long long walltime;
  walltime = (  (m_train_end.tv_sec*1000000 + m_train_end.tv_usec)
              - (m_train_start.tv_sec*1000000 + m_train_start.tv_usec) );
  m_train_wctime += (double) walltime / 1000000.0;

#ifndef serial
  MPI_Allreduce(  MPI_IN_PLACE,
                  &m_train_wctime,
                  1,
                  MPI_DOUBLE,
                  MPI_MAX,
                  MPI_COMM_WORLD );
#endif
}

/*! Project initial solution for prediction */
void libROMInterface::projectInitialSolution( void* a_s  /*!< Array of simulation objects of
                                                              type #SimulationObject */ )
{
  if (m_U.size() != m_rom.size()) {
    printf("ERROR in libROMInterface::projectInitialSolution(): m_U.size != m_rom.size() on rank %d!\n",
            m_rank );
  }

  copyFromHyPar( m_U, a_s );
  for (int i = 0; i < m_rom.size(); i++) {
    m_rom[i]->projectInitialSolution( *(m_U[i]) );
  }

  return;
}

/*! Train the ROM object */
void libROMInterface::train()
{
  gettimeofday(&m_train_start, NULL);
  for (int i = 0; i < m_rom.size(); i++) {
    m_rom[i]->train();
  }
  gettimeofday(&m_train_end, NULL);

  long long walltime;
  walltime = (  (m_train_end.tv_sec*1000000 + m_train_end.tv_usec)
              - (m_train_start.tv_sec*1000000 + m_train_start.tv_usec) );
  m_train_wctime += (double) walltime / 1000000.0;
#ifndef serial
  MPI_Allreduce(  MPI_IN_PLACE,
                  &m_train_wctime,
                  1,
                  MPI_DOUBLE,
                  MPI_MAX,
                  MPI_COMM_WORLD );
#endif
}

/*! Predict the solution at a given time */
void libROMInterface::predict(void*  a_s, /*!< Array of simulation objects of type #SimulationObject */
                              const double a_t  /*!< time at which to predict solution */) const
{
  m_predict_wctime = 0.0;

  if (m_comp_mode == _ROM_COMP_MODE_MONOLITHIC_) {

    for (int ns = 0; ns < m_nsims; ns++) {
      gettimeofday(&m_predict_start, NULL);
      const CAROM::Vector* const u_predicted = m_rom[ns]->predict(a_t);
      gettimeofday(&m_predict_end, NULL);
      copyToHyPar( *u_predicted, a_s, ns );

      long long walltime;
      walltime = (  (m_predict_end.tv_sec*1000000 + m_predict_end.tv_usec)
                  - (m_predict_start.tv_sec*1000000 + m_predict_start.tv_usec) );
      m_predict_wctime += (double) walltime / 1000000.0;
    }

  } else if (m_comp_mode == _ROM_COMP_MODE_COMPONENTWISE_) {

    int count(0);
    for (int ns = 0; ns < m_nsims; ns++) {
      for (int v = 0; v < m_ncomps[ns]; v++) {
        gettimeofday(&m_predict_start, NULL);
        const CAROM::Vector* const u_predicted = m_rom[count]->predict(a_t);
        gettimeofday(&m_predict_end, NULL);
        copyToHyPar( *u_predicted, a_s, ns, v );

        long long walltime;
        walltime = (  (m_predict_end.tv_sec*1000000 + m_predict_end.tv_usec)
                    - (m_predict_start.tv_sec*1000000 + m_predict_start.tv_usec) );
        m_predict_wctime += (double) walltime / 1000000.0;
        count++;
      }
    }

  }

#ifndef serial
  MPI_Allreduce(  MPI_IN_PLACE,
                  &m_predict_wctime,
                  1,
                  MPI_DOUBLE,
                  MPI_MAX,
                  MPI_COMM_WORLD );
#endif
}

/*! Save ROM object to file */
void libROMInterface::saveROM(const std::string& a_fname_root /*!< filename root */) const
{
  if (m_save_ROM) {

    if (!m_rank) {
      printf("libROMInterface::saveROM() - saving ROM objects.\n");
    }

    if (m_comp_mode == _ROM_COMP_MODE_MONOLITHIC_) {

      for (int ns = 0; ns < m_nsims; ns++) {
        std::string fname_root = a_fname_root;
        if (m_nsims > 1) {
          char idx_string[_MAX_STRING_SIZE_];
          sprintf(idx_string, "sim%03d", ns);
          fname_root += std::string(idx_string);
        }
        m_rom[ns]->save(fname_root);
      }

    } else if (m_comp_mode == _ROM_COMP_MODE_COMPONENTWISE_) {

      int count(0);
      for (int ns = 0; ns < m_nsims; ns++) {
        for (int v = 0; v < m_ncomps[ns]; v++) {
          std::string fname_root = a_fname_root;
          if (m_nsims > 1) {
            char idx_string[_MAX_STRING_SIZE_];
            sprintf(idx_string, "sim%03d", ns);
            fname_root += std::string(idx_string);
          }
          if (m_ncomps[ns] > 1) {
            char idx_string[_MAX_STRING_SIZE_];
            sprintf(idx_string, "var%03d", v);
            fname_root += std::string(idx_string);
          }
          m_rom[count]->save(fname_root);
          count++;
        }
      }

    }

  }

  return;
}

/*! load ROM object from file */
void libROMInterface::loadROM(const std::string& a_fname_root /*!< filename root */)
{
  if (!m_rank) {
    printf("libROMInterface::loadROM() - loading ROM objects.\n");
  }

  if (m_comp_mode == _ROM_COMP_MODE_MONOLITHIC_) {

    for (int ns = 0; ns < m_nsims; ns++) {
      std::string fname_root = a_fname_root;
      if (m_nsims > 1) {
        char idx_string[_MAX_STRING_SIZE_];
        sprintf(idx_string, "sim%03d", ns);
        fname_root += std::string(idx_string);
      }
      m_rom[ns]->load(fname_root);
    }

  } else if (m_comp_mode == _ROM_COMP_MODE_COMPONENTWISE_) {

    int count(0);
    for (int ns = 0; ns < m_nsims; ns++) {
      for (int v = 0; v < m_ncomps[ns]; v++) {
        std::string fname_root = a_fname_root;
        if (m_nsims > 1) {
          char idx_string[_MAX_STRING_SIZE_];
          sprintf(idx_string, "sim%03d", ns);
          fname_root += std::string(idx_string);
        }
        if (m_ncomps[ns] > 1) {
          char idx_string[_MAX_STRING_SIZE_];
          sprintf(idx_string, "var%03d", v);
          fname_root += std::string(idx_string);
        }
        m_rom[count]->load(fname_root);
        count++;
      }
    }

  }

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

  if (m_comp_mode == _ROM_COMP_MODE_MONOLITHIC_) {

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

  } else if (m_comp_mode == _ROM_COMP_MODE_COMPONENTWISE_) {

    int count(0);
    for (int ns = 0; ns < m_nsims; ns++) {
      for (int v = 0; v < sim[ns].solver.nvars; v++) {
        double* vec = a_U[count]->getData();
        const double* u = sim[ns].solver.u;

        std::vector<int> index(sim[ns].solver.ndims);

        ArrayCopynDComponent( sim[ns].solver.ndims,
                              u,
                              vec,
                              sim[ns].solver.dim_local,
                              sim[ns].solver.ghosts,
                              0,
                              index.data(),
                              sim[ns].solver.nvars,
                              1,
                              v,
                              0 );
        count++;
      }
    }

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
  if (m_mode == _ROM_MODE_TRAIN_) {
    u = sim[a_idx].solver.u_rom_predicted;
  } else {
    u = sim[a_idx].solver.u;
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

/*! Copy a vector a_vec to a component of the HyPar solution. Note that the HyPar solution
 * has ghost points but the vector a_vec is a serialized vector of one of the solution components
 * without ghost points.
 *
 * \b Note: if libROMInterface::m_mode is "predict", then the vector will be copied
 * into the main solution variable HyPar::u; otherwise it will be copied into the
 * variable for ROM prediction HyPar::u_rom_predicted */
void libROMInterface::copyToHyPar(  const CAROM::Vector& a_vec,  /*!< Work vector */
                                    void* a_s, /*!< Array of simulation objects of
                                                    type #SimulationObject */
                                    int a_idx /*!< Simulation object index */,
                                    int a_var /*!< Vector component to copy */ ) const
{
  SimulationObject* sim = (SimulationObject*) a_s;

  double* u;
  if (m_mode == _ROM_MODE_TRAIN_) {
    u = sim[a_idx].solver.u_rom_predicted;
  } else {
    u = sim[a_idx].solver.u;
  }
  std::vector<int> index(sim[a_idx].solver.ndims);

  ArrayCopynDComponent(  sim[a_idx].solver.ndims,
                        a_vec.getData(),
                        u,
                        sim[a_idx].solver.dim_local,
                        0,
                        sim[a_idx].solver.ghosts,
                        index.data(),
                        1,
                        sim[a_idx].solver.nvars,
                        0,
                        a_var );

  return;
}

#endif

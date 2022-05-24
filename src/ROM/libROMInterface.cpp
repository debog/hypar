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

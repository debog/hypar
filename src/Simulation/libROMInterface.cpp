/*! @file libROMInterface.cpp
 *  @brief Functions that implement interface with libROM
 *  @author Debojyoti Ghosh
*/

#ifdef with_librom

#include <arrayfunctions.h>
#include <simulation_object.h>
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

    where the list of keywords and their type are:\n
    Keyword name       | Type         | Variable                                      | Default value
    ------------------ | ------------ | --------------------------------------------- | -------------------
    rdim               | int          | #libROMInterface::m_rdim                      | 10
   
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

  if (!m_rank) {

    FILE *in;
    in = fopen(_LIBROM_INP_FNAME_,"r");

    if (!in) {

      fprintf(stderr, "Error in libROMInterface.cpp::Define() -\n");
      fprintf(stderr, "  %s file not found.\n", _LIBROM_INP_FNAME_);
      exit(1);

    } else {

      int ferr;
      char word[_MAX_STRING_SIZE_];
      ferr = fscanf(in,"%s", word); if (ferr != 1) return;

      if (std::string(word) == "begin") {
        while (std::string(word) != "end") {
  	      ferr = fscanf(in,"%s",word); if (ferr != 1) return;
          if (std::string(word) == "rdim") {
            ferr = fscanf(in,"%d",&m_rdim); if (ferr != 1) return;
          } else if (std::string(word) != "end") {
            char useless[_MAX_STRING_SIZE_];
            ferr = fscanf(in,"%s",useless);
            printf("Warning: keyword %s in file \"%s\" with value %s not recognized or extraneous. Ignoring.\n",
                    _LIBROM_INP_FNAME_, word, useless );
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
    printf("  local vector size: %d\n", m_vec_size);
  }

#ifndef serial
  MPI_Bcast(&m_rdim,1,MPI_INT,0,MPI_COMM_WORLD);
#endif

  m_rom = new CAROM::DMD( m_vec_size, a_dt );
  m_U = new double[m_vec_size];

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
    double* u = sim[ns].solver.u;

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

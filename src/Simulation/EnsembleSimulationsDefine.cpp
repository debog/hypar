/*! @file EnsembleSimulationsDefine.cpp
    @brief Define an ensemble simulation object
    @author Debojyoti Ghosh
*/

#include <string>
#include <ensemble_simulations.h>

/*! Define the ensemble simulation object:
    This function also reads sparse grids inputs from the file
    \b simulation.inp. Rank 0 reads in the inputs and broadcasts
    them to all the processors.\n\n
    The format of \b solver.inp is as follows:\n

        begin
            <keyword>   <value>
            <keyword>   <value>
            ...
            <keyword>   <value>
        end

    where the list of keywords and their type are:\n
    Keyword name   | Type    | Variable                       | Default value
    -------------- | ------- | ------------------------------ | ----------------
    nsims          | int     | #EnsembleSimulation::m_nsims   | -1

*/
int EnsembleSimulation::define( int a_rank, /*!< MPI rank of this process */
                                int a_nproc /*!< Total number of MPI ranks */
                              )
{
  if (m_is_defined) {
    fprintf(stderr,"Error: object already defined on rank %d.\n", a_rank);
    return 1;
  }

  m_rank = a_rank;
  m_nproc = a_nproc;

  /* default value */
  m_nsims = -1;

  if (!m_rank) {

    FILE *in;
    in = fopen(_ENSEMBLE_SIM_INP_FNAME_,"r");

    if (!in) {
      fprintf(stderr, "Error in EnsembleSimulations::Define() - %s file not found.\n",
              _ENSEMBLE_SIM_INP_FNAME_);
    } else {

      int ferr;
      char word[_MAX_STRING_SIZE_];
      ferr = fscanf(in,"%s", word); if (ferr != 1) return(1);

      if (std::string(word) == "begin") {

        while (std::string(word) != "end") {

          ferr = fscanf(in,"%s",word); if (ferr != 1) return(1);

          if (std::string(word) == "nsims") {
            ferr = fscanf(in,"%d",&m_nsims); if (ferr != 1) return(1);
          } else if (std::string(word) != "end") {
            char useless[_MAX_STRING_SIZE_];
            ferr = fscanf(in,"%s",useless);
            printf("Warning: keyword %s in file \"%s\" with value %s not recognized or extraneous. Ignoring.\n",
                    _ENSEMBLE_SIM_INP_FNAME_, word, useless );
          }

          if (ferr != 1) return(1);
        }

      } else {
        fprintf(stderr,"Error: Illegal format in file \"%s\". Word read is: %s\n",
                _ENSEMBLE_SIM_INP_FNAME_, word);
        return 1;
      }

      fclose(in);

    }

    if (m_nsims < 1) {
      fprintf(stderr,"Error in InitializeSimulation(): invalid value for nsims (%d)!\n", m_nsims);
      return 1;
    }

    printf("Number of simulation domains: %d\n", m_nsims);

  }

#ifndef serial
  MPI_Bcast(&m_nsims,1,MPI_INT,0,MPI_COMM_WORLD);
#endif

  if (m_nsims < 0) {
    return 1;
  }

  m_sims.resize(m_nsims);
  for (int ns = 0; ns < m_nsims; ns++) {
    m_sims[ns].solver.my_idx = ns;
    m_sims[ns].solver.nsims = m_nsims;
    m_sims[ns].mpi.rank = m_rank;
    m_sims[ns].mpi.nproc = m_nproc;
  }

  if (!m_rank) {
    printf("Allocated simulation object(s).\n");
  }

  m_is_defined = true;
  return 0;
}

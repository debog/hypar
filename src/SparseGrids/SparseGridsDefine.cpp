/*! @file SparseGridsDefine.cpp
    @brief Define a single simulation object
    @author Debojyoti Ghosh, John Loffeld, Lee Ricketson
*/

#include <sparse_grids_simulation.h>

/*! Define the sparse grids simulation object - here, only the 
 * full grid simulation object #SparseGridsSimulation::m_sims_fg 
 * is created. */
int SparseGridsSimulation::define(  int a_rank, /*!< MPI rank of this process */
                                    int a_nproc /*!< Total number of MPI ranks */ 
                                 )
{
  m_rank = a_rank;
  m_nproc = a_nproc;

  if (m_sim_fg != NULL) {
    fprintf(stderr, "Error in SparseGridsSimulation::define() - m_sim_fg not NULL!\n");
    return 1;
  }

  /* default values */
  m_imin = 2;

  if (!m_rank) {

    FILE *in;
    in = fopen(_SPARSEGRIDS_SIM_INP_FNAME_,"r");

    if (!in) {
      fprintf(stderr, "Error in EnsembleSimulations::Define() - %s file not found.\n",
              _SPARSEGRIDS_SIM_INP_FNAME_);
    } else {

      int ferr;
      char word[_MAX_STRING_SIZE_];
      ferr = fscanf(in,"%s", word); if (ferr != 1) return(1);

      if (std::string(word) == "begin") {

        while (std::string(word) != "end") {

  	      ferr = fscanf(in,"%s",word); if (ferr != 1) return(1);

          if (std::string(word) == "log2_imin") {
            ferr = fscanf(in,"%d",&m_imin); if (ferr != 1) return(1);
          } else if (std::string(word) != "end") {
            char useless[_MAX_STRING_SIZE_];
            ferr = fscanf(in,"%s",useless);
            printf("Warning: keyword %s in file \"%s\" with value %s not recognized or extraneous. Ignoring.\n",
                    _SPARSEGRIDS_SIM_INP_FNAME_, word, useless );
          }

          if (ferr != 1) return(1);
        }

      } else {
   		  fprintf(stderr,"Error: Illegal format in file \"%s\". Word read is: \n",
                _SPARSEGRIDS_SIM_INP_FNAME_, word);
        return 1;
      }

      fclose(in);

    }

  }

#ifndef serial
  MPI_Bcast(&m_imin,1,MPI_INT,0,MPI_COMM_WORLD);
#endif

  m_sim_fg = new SimulationObject;
  /* the following are deliberately set to junk values
   * to ensure they are never used */
  m_sim_fg->solver.my_idx = -1;
  m_sim_fg->solver.nsims = -1;
  /* these are okay */
  m_sim_fg->mpi.rank = m_rank;
  m_sim_fg->mpi.nproc = m_nproc;

  if (!m_rank) {
    printf("Allocated full grid simulation object(s).\n");
  }

  m_is_defined = true;
  return 0;
}

/*! @file mpivars_struct.h
    @brief MPI related structure.
    @author Debojyoti Ghosh
 */

#ifndef _MPIVARS_STRUCT_H_
#define _MPIVARS_STRUCT_H_

#ifndef serial
#include <mpi.h>
#endif

/*! \def MPIVariables
 *  \brief Structure of MPI-related variables.
 * This structure contains all the variables needed for parallel computations
 * using the MPI library.
*/

/*! \brief Structure of MPI-related variables.
 *
 * This structure contains all the variables needed for parallel computations
 * using the MPI library.
*/
typedef struct mpi_variables {
  int   m_rank;     /*!< Process rank                                       */
  int   m_nproc;    /*!< Total number of processes                          */
  int   *m_iproc;   /*!< Number of processes along each dimension           */
  int   *m_ip;      /*!< Process rank along each dimension                  */
  int   *m_is,      /*!< Global start index of local domain along each dimension  */
        *m_ie;      /*!< Global end index of local domain along each dimension  */
  int   *m_bcperiodic; /*!< Flag for periodic BCs along any dimension       */

#ifdef serial
  int   m_world; /* Dummy variable */
  int   *m_comm; /* Dummy variable */
#else
  MPI_Comm  m_world;   /*!< Communicator for all processes                  */
  MPI_Comm  *m_comm;   /*!< Sub-communicators                               */
#endif

  int m_N_IORanks;      /*!< Number of IO ranks                             */
  int m_IOParticipant;  /*!< Whether this rank will handle file I/O         */
  int m_CommGroup;      /*!< I/O group this rank is a part of               */
  int m_IORank       ;  /*!< Rank of the process this rank will get I/O from*/
  int m_GroupStartRank; /*!< Starting rank of the IO group                  */
  int m_GroupEndRank;   /*!< Last rank of the IO group                      */
#ifndef serial
  MPI_Comm m_IOWorld;   /*!< Communicator of processes participating in file I/O */
#endif

  double *m_sendbuf, /*!< Buffer to send data */
         *m_recvbuf; /*!< Buffer to receive data */
  int    m_maxbuf;   /*!< Maximum buffer size */

#if defined(HAVE_CUDA)
  int    m_ncalls;
  double m_wctime;
  double m_wctime_total;
  int    *m_cpu_dim;
  double *m_gpu_sendbuf,
         *m_gpu_recvbuf;
#endif
} MPIVariables;

#endif

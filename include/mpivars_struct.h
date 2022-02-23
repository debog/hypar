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
  int   rank;     /*!< Process rank                                       */
  int   nproc;    /*!< Total number of processes                          */
  int   *iproc;   /*!< Number of processes along each dimension           */
  int   *ip;      /*!< Process rank along each dimension                  */
  int   *is,      /*!< Global start index of local domain along each dimension  */
        *ie;      /*!< Global end index of local domain along each dimension  */
  int   *bcperiodic; /*!< Flag for periodic BCs along any dimension       */

#ifdef serial
  int   world; /* Dummy variable */
  int   *comm; /* Dummy variable */
#else
  MPI_Comm  world;   /*!< Communicator for all processes                  */
  MPI_Comm  *comm;   /*!< Sub-communicators                               */
#endif

  int N_IORanks;      /*!< Number of IO ranks                             */
  int IOParticipant;  /*!< Whether this rank will handle file I/O         */
  int CommGroup;      /*!< I/O group this rank is a part of               */
  int IORank       ;  /*!< Rank of the process this rank will get I/O from*/
  int GroupStartRank; /*!< Starting rank of the IO group                  */
  int GroupEndRank;   /*!< Last rank of the IO group                      */
#ifndef serial
  MPI_Comm IOWorld;   /*!< Communicator of processes participating in file I/O */
#endif

  double *sendbuf, /*!< Buffer to send data */
         *recvbuf; /*!< Buffer to receive data */
  int    maxbuf;   /*!< Maximum buffer size */

#if defined(HAVE_CUDA)
  int    ncalls;
  double wctime;
  double wctime_total;
  int    *cpu_dim;
  double *gpu_sendbuf,
         *gpu_recvbuf;
#endif
} MPIVariables;

#endif

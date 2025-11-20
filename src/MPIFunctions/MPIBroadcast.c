/*! @file MPIBroadcast.c
    @brief Functions to broadcast over all MPI ranks
    @author Debojyoti Ghosh
*/

#include <mpivars.h>

/*! Broadcast an array of type \a double to all MPI ranks */
int MPIBroadcast_double(
                          double  *a_x,     /*!< array to broadcast to all ranks */
                          int     a_size,   /*!< size of array to broadcast */
                          int     a_root,   /*!< rank from which to broadcast */
                          void    *comm   /*!< MPI communicator within which to broadcast */
                       )
{
#ifndef serial
  MPI_Bcast(a_x,a_size,MPI_DOUBLE,a_root,*((MPI_Comm*)a_comm));
#endif
  return(0);
}

/*! Broadcast an array of type \a int to all MPI ranks */
int MPIBroadcast_integer(
                          int   *a_x,     /*!< array to broadcast to all ranks */
                          int   a_size,   /*!< size of array to broadcast */
                          int   a_root,   /*!< rank from which to broadcast */
                          void  *comm   /*!< MPI communicator within which to broadcast */
                        )
{
#ifndef serial
  MPI_Bcast(a_x,a_size,MPI_INT,a_root,*((MPI_Comm*)a_comm));
#endif
  return(0);
}

/*! Broadcast an array of type \a char to all MPI ranks */
int MPIBroadcast_character(
                            char  *a_x,   /*!< array to broadcast to all ranks */
                            int   a_size, /*!< size of array to broadcast */
                            int   a_root, /*!< rank from which to broadcast */
                            void  *comm /*!< MPI communicator within which to broadcast */
                          )
{
#ifndef serial
  MPI_Bcast(a_x,a_size,MPI_CHAR,a_root,*((MPI_Comm*)a_comm));
#endif
  return(0);
}

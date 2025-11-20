/*! @file MPIMax.c
    @brief Functions to compute the maximum over all MPI ranks
    @author Debojyoti Ghosh
*/

#include <mpivars.h>

/*!
  Compute the global maximum over all MPI ranks in a given communicator for
  \a int datatype.
  + If \a var is an array of size greater than 1, \a global will be an array
    of the same size with each element as the maximum value of that element
    in \a var on all the MPI ranks in the given communicator.
*/
int MPIMax_integer(
                    int   *a_global, /*!< array to contain the global maximums */
                    int   *a_var,    /*!< the local array */
                    int   a_size,    /*!< size of the local array */
                    void  *comm    /*!< MPI communicator */
                  )
{
#ifdef serial
  int i;
  for (i = 0; i < a_size; i++)  a_global[i] = a_var[i];
#else
  MPI_Allreduce((a_var==a_global?MPI_IN_PLACE:a_var),a_global,a_size,MPI_INT,MPI_MAX,*((MPI_Comm*)a_comm));
#endif
  return(0);
}

/*!
  Compute the global maximum over all MPI ranks in a given communicator for
  \a double datatype.
  + If \a var is an array of size greater than 1, \a global will be an array
    of the same size with each element as the maximum value of that element
    in \a var on all the MPI ranks in the given communicator.
*/
int MPIMax_double(
                    double  *a_global, /*!< array to contain the global maximums */
                    double  *a_var,    /*!< the local array */
                    int     a_size,    /*!< size of the local array */
                    void    *comm    /*!< MPI communicator */
                 )
{
#ifdef serial
  int i;
  for (i = 0; i < a_size; i++)  a_global[i] = a_var[i];
#else
  MPI_Allreduce((a_var==a_global?MPI_IN_PLACE:a_var),a_global,a_size,MPI_DOUBLE,MPI_MAX,*((MPI_Comm*)a_comm));
#endif
  return(0);
}

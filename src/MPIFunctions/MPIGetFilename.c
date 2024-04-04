/*! @file MPIGetFilename.c
    @brief Get filename indexed by MPI rank
    @author Debojyoti Ghosh
*/

#ifndef serial
#include <mpi.h>
#endif
#include <string.h>
#include <basic.h>
#include <mpivars.h>
#include <common.h>

/*!
    Get a string representing a filename indexed by the MPI rank:
    \a filename = \a root.index, where \a index is the string
    corresponding to the MPI rank.
*/
void MPIGetFilename(
                      char  *root,      /*!< filename root */
                      void  *c,         /*!< MPI communicator */
                      char  *filename   /*!< filename */
                   )
{
  char  tail[_MAX_STRING_SIZE_]="";
  int   rank;

#ifndef serial
  MPI_Comm  comm = *((MPI_Comm*)c);
  MPI_Comm_rank(comm,&rank);
#else
  rank = 0;
#endif

  GetStringFromInteger(rank,tail,4);
  strcpy(filename,"");
  strcat(filename,root);
  strcat(filename,"." );
  strcat(filename,tail);

  return;
}

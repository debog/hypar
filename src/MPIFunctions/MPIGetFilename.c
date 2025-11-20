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
                      char  *a_root,      /*!< filename root */
                      void  *a_c,         /*!< MPI communicator */
                      char  *a_filename   /*!< filename */
                   )
{
  char  tail[_MAX_STRING_SIZE_]="";
  int   rank;

#ifndef serial
  MPI_Comm  comm = *((MPI_Comm*)a_c);
  MPI_Comm_rank(comm,&rank);
#else
  rank = 0;
#endif

  GetStringFromInteger(rank,tail,4);
  strcpy(a_filename,"");
  strcat(a_filename,a_root);
  strcat(a_filename,"." );
  strcat(a_filename,tail);

  return;
}

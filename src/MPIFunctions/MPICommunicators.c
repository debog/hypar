/*! @file MPICommunicators.c
    @brief Functions to create and destroy MPI subcommunicators
    @author Debojyoti Ghosh
*/

#include <stdlib.h>
#include <arrayfunctions.h>
#include <mpivars.h>

/*!
  Create subcommunicators from MPI_WORLD, where each subcommunicator contains
  MPI ranks along a spatial dimension. Consider a two-dimensional problem,
  partitioned on 21 MPI ranks as follows:
  @image html mpi_ranks.png
  @image latex mpi_ranks.eps width=0.9\textwidth

  This function will create 10 subcommunicators with the following ranks:
  + 0,1,2,3,4,5,6
  + 7,8,9,10,11,12,13
  + 14,15,16,17,18,19,20
  + 0,7,14
  + 1,8,15
  + 2,9,16
  + 3,10,17
  + 4,11,18
  + 5,12,19
  + 6,13,20

  These subcommunicators are useful for parallel computations along
  grid lines. For example, a compact finite-difference scheme solves
  implicit systems along grid lines in every spatial dimension. Thus,
  the subcommunicator may be passed on to the parallel systems solver
  instead of MPI_WORLD.
*/
int MPICreateCommunicators(
                            int   a_ndims,  /*!< Number of spatial dimensions */
                            void  *a_m      /*!< MPI object of type #MPIVariables */
                          )
{
  MPIVariables *mpi = (MPIVariables*) a_m;
#ifdef serial
  mpi->m_comm = NULL;
#else
  int          i,n,color,key;
  int          *ip,*iproc;

  mpi->m_comm = (MPI_Comm*) calloc (a_ndims, sizeof(MPI_Comm));
  if (a_ndims == 1) MPI_Comm_dup(mpi->m_world,mpi->m_comm);
  else {
    ip    = (int*) calloc (a_ndims-1,sizeof(int));
    iproc = (int*) calloc (a_ndims-1,sizeof(int));
    for (n=0; n<a_ndims; n++) {
      int tick=0;
      for (i=0; i<a_ndims; i++) {
        if (i != n) {
          ip[tick]    = mpi->m_ip[i];
          iproc[tick] = mpi->m_iproc[i];
          tick++;
        }
      }
      _ArrayIndex1D_(a_ndims-1,iproc,ip,0,color);
      key   = mpi->m_ip[n];
      MPI_Comm_split(mpi->m_world,color,key,&mpi->m_comm[n]);
    }
    free(ip);
    free(iproc);
  }
#endif
  return(0);
}

/*!
  Free the subcommunicators created in MPICreateCommunicators().
*/
int MPIFreeCommunicators(
                          int   a_ndims,  /*!< Number of spatial dimensions */
                          void  *a_m      /*!< MPI object of type #MPIVariables */
                        )
{
#ifndef serial
  MPIVariables *mpi = (MPIVariables*) a_m;
  int          n;
  for (n=0; n<a_ndims; n++) MPI_Comm_free(&mpi->m_comm[n]);
  free(mpi->m_comm);
  if (mpi->m_IOParticipant) MPI_Comm_free(&mpi->m_IOWorld);
  MPI_Comm_free(&mpi->m_world);
#endif
  return(0);
}

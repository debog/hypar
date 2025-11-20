/*! @file MPIExchangeBoundaries1D.c
    @brief Exchange data and fill ghost points for a 1D array
    @author Debojyoti Ghosh
*/

#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mpivars.h>

/*!
  Exchange the data across MPI ranks and fill in ghost points for a 1D array. In a multidimensional
  simulation, a 1D array is an array of data along one of the spatial dimensions, i.e. its an array
  storing a variable that varies in only one of the spatial dimension. For example, for a
  2D problem on a Cartesian grid (with spatial dimensions x and y), the array of x-coordinates is
  a 1D array along x, and the array of y-coordinates is a 1D array along y. Thus, the size of the
  1D array is equal to the size of the domain along the spatial dimension corresponding to that array.

  Consider a two-dimensional problem, partitioned on 21 MPI ranks as follows:
  @image html mpi_ranks.png
  @image latex mpi_ranks.eps width=0.9\textwidth
  and consider rank 9.

  If the argument \a dir is specified as 0, and thus we are dealing with a 1D array
  along dimension 0, then
  + Rank 9 will exchange data with ranks 8 and 10, and fill in its ghost points.

  If \a dir is specified as 1, and thus we are dealing with a 1D array along dimension
  1, then
  + Rank 9 will exchange data with ranks 2 and 16, and fill in its ghost points.
*/
int MPIExchangeBoundaries1D(
                              void    *a_m,     /*!< MPI object of type MPIVariables */
                              double  *a_x,     /*!< The 1D array for which to exchange data */
                              int     a_N,      /*!< Size of the array */
                              int     a_ghosts, /*!< Number of ghost points */
                              int     a_dir,    /*!< Spatial dimension corresponding to the 1D array */
                              int     ndims   /*!< Number of spatial dimensions in the simulation */
                           )
{
#ifndef serial
  MPIVariables  *mpi = (MPIVariables*) a_m;
  int           i;

  int *ip     = mpi->m_ip;
  int *iproc  = mpi->m_iproc;
  int non      = 0; /* number of neighbours */

  int neighbor_rank[2] = {-1,-1};
  int nip[a_ndims];

  /* each process has 2 neighbors (except at physical boundaries)       */
  /* calculate the rank of these neighbors (-1 -> none)                 */
  _ArrayCopy1D_(ip,nip,a_ndims);  nip[a_dir]--;
  if (ip[a_dir] == 0)             neighbor_rank[0] = -1;
  else                          neighbor_rank[0] = MPIRank1D(a_ndims,iproc,nip);
  _ArrayCopy1D_(ip,nip,a_ndims);  nip[a_dir]++;
  if (ip[a_dir] == (iproc[a_dir]-1))neighbor_rank[1] = -1;
  else                          neighbor_rank[1] = MPIRank1D(a_ndims,iproc,nip);

  /* Allocate send and receive buffers */
  double sendbuf[2][a_ghosts], recvbuf[2][a_ghosts];

  /* count number of neighbors and copy data to send buffers */
  non = 0;
  if (neighbor_rank[0] != -1) {
    non++;
    for (i = 0; i < a_ghosts; i++) sendbuf[0][i] = a_x[i+a_ghosts];
  }
  if (neighbor_rank[1] != -1) {
    non++;
    for (i = 0; i < a_ghosts; i++) sendbuf[1][i] = a_x[i+a_N];
  }
  MPI_Request requests[2*non];
  MPI_Status  statuses[2*non];

  /* exchange the data */
  int tick = 0;
  if (neighbor_rank[0]!= -1) {
    MPI_Irecv(recvbuf[0],a_ghosts,MPI_DOUBLE,neighbor_rank[0],1631,mpi->m_world,&requests[tick]);
    MPI_Isend(sendbuf[0],a_ghosts,MPI_DOUBLE,neighbor_rank[0],1631,mpi->m_world,&requests[tick+non]);
    tick++;
  }
  if (neighbor_rank[1] != -1) {
    MPI_Irecv(recvbuf[1],a_ghosts,MPI_DOUBLE,neighbor_rank[1],1631,mpi->m_world,&requests[tick]);
    MPI_Isend(sendbuf[1],a_ghosts,MPI_DOUBLE,neighbor_rank[1],1631,mpi->m_world,&requests[tick+non]);
    tick++;
  }

  /* Wait till data transfer is done */
  MPI_Waitall(2*non,requests,statuses);

  /* copy received data to ghost points */
  if (neighbor_rank[0] != -1) for (i = 0; i < a_ghosts; i++) a_x[i]          = recvbuf[0][i];
  if (neighbor_rank[1] != -1) for (i = 0; i < a_ghosts; i++) a_x[i+a_N+a_ghosts] = recvbuf[1][i];

#endif
  return(0);
}


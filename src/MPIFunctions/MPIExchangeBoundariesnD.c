/*! @file MPIExchangeBoundariesnD.c
    @brief Exchange data and fill in ghost points for an n-dimensional array
    @author Debojyoti Ghosh
*/

#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mpivars.h>

/*!
  Exchange data across MPI ranks, and fill in ghost points for an n-dimensional array
  (where \a n is the total number of spatial dimensions). If any of the physical boundaries
  are periodic, this function also exchanges data and fills in the ghost points for these
  boundaries.

  The n-dimensional array must be stored in the memory as a single-index array, with the following order of mapping:
  + Number of variables (vector components)
  + Spatial dimension 0
  + Spatial dimension 1
  + ...
  + Spatial dimensions \a ndims-1

  For example, consider a 2D simulation (\a ndims = 2), of size \f$7 \times 3\f$, with \f$4\f$ vector components
  (\a nvars = 4). The following figure shows the layout (without the ghost points):
  @image html layout.png
  @image latex layout.eps width=0.9\textwidth

  The bold numbers in parentheses represent the 2D indices. The numbers below them are the indices of the array
  that correspond to that 2D location. Thus, elements 40,41,42, and 43 in the array are the 1st, 2nd, 3rd, and
  4th vector components at location (1,3).

  If \f${\bf i}\left[{\rm ndims}\right]\f$ is an integer array representing an n-dimensional index
  (for example, \f$\left(5,4\right)\f$ in 2D, \f$\left(3,5,2\right)\f$ in 3D), and the number of vector
  components is \a nvars, then:
  + #_ArrayIndex1D_ computes the index \f$p\f$ in the array corresponding to \f${\bf i}\f$. In the above example,
    \f${\bf i} = \left(1,3\right) \rightarrow p = 10\f$.
  + \a var[nvars*p+v] accesses the \a v-th component of the n-dimensional array \a var at location \f${\bf i}\f$.
    In the above example, to access the 3rd vector component at location \f$\left(1,3\right)\f$, we have \f$p=10\f$,
    so \a var [4*10+2] = \a var [42].
*/
int MPIExchangeBoundariesnD(
                              int     a_ndims,  /*!< Number of spatial dimensions */
                              int     a_nvars,  /*!< Number of variables (vector components) at each grid location */
                              int     *a_dim,   /*!< Integer array whose elements are the local size along each spatial dimension */
                              int     a_ghosts, /*!< Number of ghost points */
                              void    *a_m,     /*!< MPI object of type #MPIVariables */
                              double  *var    /*!< The array for which to exchange data and fill in ghost points */
                           )
{
#ifndef serial
  MPIVariables  *mpi = (MPIVariables*) a_m;
  int           d;

  int *ip     = mpi->m_ip;
  int *iproc  = mpi->m_iproc;
  int *bcflag = mpi->m_bcperiodic;

  int neighbor_rank[2*a_ndims], nip[a_ndims], index[a_ndims], bounds[a_ndims], offset[a_ndims];
  MPI_Request rcvreq[2*a_ndims], sndreq[2*a_ndims];
  for (d=0; d<2*a_ndims; d++) rcvreq[d] = sndreq[d] = MPI_REQUEST_NULL;

  /* each process has 2*a_ndims neighbors (except at non-periodic physical boundaries)  */
  /* calculate the rank of these neighbors (-1 -> none)                               */
  for (d = 0; d < a_ndims; d++) {
    _ArrayCopy1D_(ip,nip,a_ndims);
    if (ip[d] == 0) nip[d] = iproc[d]-1;
    else            nip[d]--;
    if ((ip[d] == 0) && (!bcflag[d])) neighbor_rank[2*d]   = -1;
    else                              neighbor_rank[2*d]   = MPIRank1D(a_ndims,iproc,nip);
    _ArrayCopy1D_(ip,nip,a_ndims);
    if (ip[d] == (iproc[d]-1)) nip[d] = 0;
    else                       nip[d]++;
    if ((ip[d] == (iproc[d]-1)) && (!bcflag[d]))  neighbor_rank[2*d+1] = -1;
    else                                          neighbor_rank[2*d+1] = MPIRank1D(a_ndims,iproc,nip);
  }

  /* calculate dimensions of each of the send-receive regions */
  double *sendbuf = mpi->m_sendbuf;
  double *recvbuf = mpi->m_recvbuf;
  int    stride   = mpi->m_maxbuf;
  int    bufdim[a_ndims];
  for (d = 0; d < a_ndims; d++) {
    bufdim[d] = 1;
    int i;
    for (i = 0; i < a_ndims; i++) {
      if (i == d) bufdim[d] *= a_ghosts;
      else        bufdim[d] *= a_dim[i];
    }
  }

  /* post the receive requests */
  for (d = 0; d < a_ndims; d++) {
    if (neighbor_rank[2*d  ] != -1) {
      MPI_Irecv(&recvbuf[2*d*stride],bufdim[d]*a_nvars,MPI_DOUBLE,neighbor_rank[2*d  ],1630,
                mpi->m_world,&rcvreq[2*d]);
    }
    if (neighbor_rank[2*d+1] != -1) {
      MPI_Irecv(&recvbuf[(2*d+1)*stride],bufdim[d]*a_nvars,MPI_DOUBLE,neighbor_rank[2*d+1],1631,
                mpi->m_world,&rcvreq[2*d+1]);
    }
  }

  /* count number of neighbors and copy data to send buffers */
  for (d = 0; d < a_ndims; d++) {
    _ArrayCopy1D_(a_dim,bounds,a_ndims); bounds[d] = a_ghosts;
    if (neighbor_rank[2*d] != -1) {
      _ArraySetValue_(offset,a_ndims,0);
      int done = 0; _ArraySetValue_(index,a_ndims,0);
      while (!done) {
        int p1; _ArrayIndex1DWO_(a_ndims,a_dim,index,offset,a_ghosts,p1);
        int p2; _ArrayIndex1D_(a_ndims,bounds,index,0,p2);
        _ArrayCopy1D_((a_var+a_nvars*p1),(sendbuf+2*d*stride+a_nvars*p2),a_nvars);
        _ArrayIncrementIndex_(a_ndims,bounds,index,done);
      }
    }
    if (neighbor_rank[2*d+1] != -1) {
      _ArraySetValue_(offset,a_ndims,0);offset[d] = a_dim[d]-a_ghosts;
      int done = 0; _ArraySetValue_(index,a_ndims,0);
      while (!done) {
        int p1; _ArrayIndex1DWO_(a_ndims,a_dim,index,offset,a_ghosts,p1);
        int p2; _ArrayIndex1D_(a_ndims,bounds,index,0,p2);
        _ArrayCopy1D_((a_var+a_nvars*p1),(sendbuf+(2*d+1)*stride+a_nvars*p2),a_nvars);
        _ArrayIncrementIndex_(a_ndims,bounds,index,done);
      }
    }
  }

  /* send the data */
  for (d = 0; d < a_ndims; d++) {
    if (neighbor_rank[2*d  ] != -1) {
      MPI_Isend(&sendbuf[2*d*stride],bufdim[d]*a_nvars,MPI_DOUBLE,neighbor_rank[2*d  ],1631,
                mpi->m_world,&sndreq[2*d]);
    }
    if (neighbor_rank[2*d+1] != -1) {
      MPI_Isend(&sendbuf[(2*d+1)*stride],bufdim[d]*a_nvars,MPI_DOUBLE,neighbor_rank[2*d+1],1630,
                mpi->m_world,&sndreq[2*d+1]);
    }
  }

  /* Wait till data is done received */
  MPI_Status status_arr[2*a_ndims];
  MPI_Waitall(2*a_ndims,rcvreq,status_arr);
  /* copy received data to ghost points */
  for (d = 0; d < a_ndims; d++) {
    _ArrayCopy1D_(a_dim,bounds,a_ndims); bounds[d] = a_ghosts;
    if (neighbor_rank[2*d] != -1) {
      _ArraySetValue_(offset,a_ndims,0); offset[d] = -a_ghosts;
      int done = 0; _ArraySetValue_(index,a_ndims,0);
      while (!done) {
        int p1; _ArrayIndex1DWO_(a_ndims,a_dim,index,offset,a_ghosts,p1);
        int p2; _ArrayIndex1D_(a_ndims,bounds,index,0,p2);
        _ArrayCopy1D_((recvbuf+2*d*stride+a_nvars*p2),(a_var+a_nvars*p1),a_nvars);
        _ArrayIncrementIndex_(a_ndims,bounds,index,done);
      }
    }
    if (neighbor_rank[2*d+1] != -1) {
      _ArraySetValue_(offset,a_ndims,0); offset[d] = a_dim[d];
      int done = 0; _ArraySetValue_(index,a_ndims,0);
      while (!done) {
        int p1; _ArrayIndex1DWO_(a_ndims,a_dim,index,offset,a_ghosts,p1);
        int p2; _ArrayIndex1D_(a_ndims,bounds,index,0,p2);
        _ArrayCopy1D_((recvbuf+(2*d+1)*stride+a_nvars*p2),(a_var+a_nvars*p1),a_nvars);
        _ArrayIncrementIndex_(a_ndims,bounds,index,done);
      }
    }
  }
  /* Wait till send requests are complete before freeing memory */
  MPI_Waitall(2*a_ndims,sndreq,status_arr);

#endif
  return(0);
}


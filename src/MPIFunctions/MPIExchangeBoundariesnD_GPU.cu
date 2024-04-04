/*! @file MPIExchangeBoundariesnD_GPU.cu
    @brief Exchange data and fill ghost points for a 1D array
    @author Youngdae Kim
*/

#include <basic_gpu.h>
#include <arrayfunctions_gpu.h>
#include <mpivars.h>

#ifdef CUDA_VAR_ORDERDING_AOS

/*! Kernel function for gpuMPIExchangeBoundariesnD() */
__global__
void gpuFillSendBufLeft_kernel(
  int npoints_grid,
  int d,
  int ndims,
  int nvars,
  int ghosts,
  int stride,
  const int    * __restrict__ dim,
  const double * __restrict__ var,
  double       * __restrict__ sendbuf
)
{
  int p = blockDim.x * blockIdx.x + threadIdx.x;

  if (p < npoints_grid) {
    int index[GPU_MAX_NDIMS], bounds[GPU_MAX_NDIMS], offset[GPU_MAX_NDIMS];
    int p1, p2;

    _ArrayCopy1D_(dim,bounds,ndims); bounds[d] = ghosts;
    _ArraySetValue_(offset,ndims,0);
    _ArrayIndexnD_(ndims,p,bounds,index,0);
    _ArrayIndex1DWO_(ndims,dim,index,offset,ghosts,p1);
    _ArrayIndex1D_(ndims,bounds,index,0,p2);
    _ArrayCopy1D_((var+nvars*p1),(sendbuf+2*d*stride+nvars*p2),nvars);
  }

  return;
}

/*! Kernel function for gpuMPIExchangeBoundariesnD() */
__global__
void gpuFillSendBufRight_kernel(
  int npoints_grid,
  int d,
  int ndims,
  int nvars,
  int ghosts,
  int stride,
  const int    * __restrict__ dim,
  const double * __restrict__ var,
  double       * __restrict__ sendbuf
)
{
  int p = blockDim.x * blockIdx.x + threadIdx.x;

  if (p < npoints_grid) {
    int index[GPU_MAX_NDIMS], bounds[GPU_MAX_NDIMS], offset[GPU_MAX_NDIMS];
    int p1, p2;

    _ArrayCopy1D_(dim,bounds,ndims); bounds[d] = ghosts;
    _ArraySetValue_(offset,ndims,0); offset[d] = dim[d]-ghosts;
    _ArrayIndexnD_(ndims,p,bounds,index,0);
    _ArrayIndex1DWO_(ndims,dim,index,offset,ghosts,p1);
    _ArrayIndex1D_(ndims,bounds,index,0,p2);
    _ArrayCopy1D_((var+nvars*p1),(sendbuf+(2*d+1)*stride+nvars*p2),nvars);
  }

  return;
}

/*! Kernel function for gpuMPIExchangeBoundariesnD() */
__global__
void gpuFillRecvBufLeft_kernel(
  int npoints_grid,
  int d,
  int ndims,
  int nvars,
  int ghosts,
  int stride,
  const int    * __restrict__ dim,
  const double * __restrict__ recvbuf,
  double * __restrict__ var
)
{
  int p = blockDim.x * blockIdx.x + threadIdx.x;

  if (p < npoints_grid) {
    int index[GPU_MAX_NDIMS], bounds[GPU_MAX_NDIMS], offset[GPU_MAX_NDIMS];
    int p1, p2;

    _ArrayCopy1D_(dim,bounds,ndims); bounds[d] = ghosts;
    _ArraySetValue_(offset,ndims,0); offset[d] = -ghosts;
    _ArrayIndexnD_(ndims,p,bounds,index,0);
    _ArrayIndex1DWO_(ndims,dim,index,offset,ghosts,p1);
    _ArrayIndex1D_(ndims,bounds,index,0,p2);
    _ArrayCopy1D_((recvbuf+2*d*stride+nvars*p2),(var+nvars*p1),nvars);
  }

  return;
}

/*! Kernel function for gpuMPIExchangeBoundariesnD() */
__global__
void gpuFillRecvBufRight_kernel(
  int npoints_grid,
  int d,
  int ndims,
  int nvars,
  int ghosts,
  int stride,
  const int    * __restrict__ dim,
  const double * __restrict__ recvbuf,
  double       * __restrict__ var
)
{
  int p = blockDim.x * blockIdx.x + threadIdx.x;

  if (p < npoints_grid) {
    int index[GPU_MAX_NDIMS], bounds[GPU_MAX_NDIMS], offset[GPU_MAX_NDIMS];
    int p1, p2;

    _ArrayCopy1D_(dim,bounds,ndims); bounds[d] = ghosts;
    _ArraySetValue_(offset,ndims,0); offset[d] = dim[d];
    _ArrayIndexnD_(ndims,p,bounds,index,0);
    _ArrayIndex1DWO_(ndims,dim,index,offset,ghosts,p1);
    _ArrayIndex1D_(ndims,bounds,index,0,p2);
    _ArrayCopy1D_((recvbuf+(2*d+1)*stride+nvars*p2),(var+nvars*p1),nvars);
  }

  return;
}

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
extern "C" int gpuMPIExchangeBoundariesnD(
  int    ndims,
  int    nvars,
  const int * __restrict__ dim,
  int    ghosts,
  void   *m,
  double    * __restrict__ var
)
{
#ifndef serial
  MPIVariables  *mpi = (MPIVariables*) m;
  int           d;

  int *ip      = mpi->ip;
  int *iproc   = mpi->iproc;
  int *bcflag  = mpi->bcperiodic;
  int *cpu_dim = mpi->cpu_dim;

  int neighbor_rank[2*ndims], nip[ndims], bounds[ndims];
  MPI_Request rcvreq[2*ndims], sndreq[2*ndims];
  for (d=0; d<2*ndims; d++) rcvreq[d] = sndreq[d] = MPI_REQUEST_NULL;

  /* each process has 2*ndims neighbors (except at non-periodic physical boundaries)  */
  /* calculate the rank of these neighbors (-1 -> none)                               */
  for (d = 0; d < ndims; d++) {
    _ArrayCopy1D_(ip,nip,ndims);
    if (ip[d] == 0) nip[d] = iproc[d]-1;
    else            nip[d]--;
    if ((ip[d] == 0) && (!bcflag[d])) neighbor_rank[2*d]   = -1;
    else                              neighbor_rank[2*d]   = MPIRank1D(ndims,iproc,nip);
    _ArrayCopy1D_(ip,nip,ndims);
    if (ip[d] == (iproc[d]-1)) nip[d] = 0;
    else                       nip[d]++;
    if ((ip[d] == (iproc[d]-1)) && (!bcflag[d]))  neighbor_rank[2*d+1] = -1;
    else                                          neighbor_rank[2*d+1] = MPIRank1D(ndims,iproc,nip);
  }

  /* calculate dimensions of each of the send-receive regions */
  double *sendbuf = mpi->gpu_sendbuf;
  double *recvbuf = mpi->gpu_recvbuf;
  int    stride   = mpi->maxbuf;
  int    bufdim[ndims];
  for (d = 0; d < ndims; d++) {
    bufdim[d] = 1;
    int i;
    for (i = 0; i < ndims; i++) {
      if (i == d) bufdim[d] *= ghosts;
      else        bufdim[d] *= cpu_dim[i];
    }
  }

  /* post the receive requests */
  for (d = 0; d < ndims; d++) {
    if (neighbor_rank[2*d  ] != -1) {
      MPI_Irecv(&recvbuf[2*d*stride],bufdim[d]*nvars,MPI_DOUBLE,neighbor_rank[2*d  ],1630,
                mpi->world,&rcvreq[2*d]);
    }
    if (neighbor_rank[2*d+1] != -1) {
      MPI_Irecv(&recvbuf[(2*d+1)*stride],bufdim[d]*nvars,MPI_DOUBLE,neighbor_rank[2*d+1],1631,
                mpi->world,&rcvreq[2*d+1]);
    }
  }

  /* count number of neighbors and copy data to send buffers */
  for (d = 0; d < ndims; d++) {
    _ArrayCopy1D_(cpu_dim,bounds,ndims); bounds[d] = ghosts;
    if (neighbor_rank[2*d] != -1) {
      int npoints_grid = 1; for (int i=0; i<ndims; i++) npoints_grid *= bounds[i];
      int nblocks = (npoints_grid-1)/GPU_THREADS_PER_BLOCK + 1;
      gpuFillSendBufLeft_kernel<<<nblocks, GPU_THREADS_PER_BLOCK>>>(npoints_grid,d,ndims,nvars,ghosts,stride,dim,var,sendbuf);
    }
    if (neighbor_rank[2*d+1] != -1) {
      int npoints_grid = 1; for (int i=0; i<ndims; i++) npoints_grid *= bounds[i];
      int nblocks = (npoints_grid-1)/GPU_THREADS_PER_BLOCK + 1;
      gpuFillSendBufRight_kernel<<<nblocks, GPU_THREADS_PER_BLOCK>>>(npoints_grid,d,ndims,nvars,ghosts,stride,dim,var,sendbuf);
    }
  }
  cudaDeviceSynchronize();

  /* send the data */
  for (d = 0; d < ndims; d++) {
    if (neighbor_rank[2*d  ] != -1) {
      MPI_Isend(&sendbuf[2*d*stride],bufdim[d]*nvars,MPI_DOUBLE,neighbor_rank[2*d  ],1631,
                mpi->world,&sndreq[2*d]);
    }
    if (neighbor_rank[2*d+1] != -1) {
      MPI_Isend(&sendbuf[(2*d+1)*stride],bufdim[d]*nvars,MPI_DOUBLE,neighbor_rank[2*d+1],1630,
                mpi->world,&sndreq[2*d+1]);
    }
  }

  /* Wait till data is done received */
  MPI_Status status_arr[2*ndims];
  MPI_Waitall(2*ndims,rcvreq,status_arr);
  /* copy received data to ghost points */
  for (d = 0; d < ndims; d++) {
    _ArrayCopy1D_(cpu_dim,bounds,ndims); bounds[d] = ghosts;
    if (neighbor_rank[2*d] != -1) {
      int npoints_grid = 1; for (int i=0; i<ndims; i++) npoints_grid *= bounds[i];
      int nblocks = (npoints_grid-1)/GPU_THREADS_PER_BLOCK + 1;
      gpuFillRecvBufLeft_kernel<<<nblocks, GPU_THREADS_PER_BLOCK>>>(npoints_grid,d,ndims,nvars,ghosts,stride,dim,recvbuf,var);
    }
    if (neighbor_rank[2*d+1] != -1) {
      int npoints_grid = 1; for (int i=0; i<ndims; i++) npoints_grid *= bounds[i];
      int nblocks = (npoints_grid-1)/GPU_THREADS_PER_BLOCK + 1;
      gpuFillRecvBufRight_kernel<<<nblocks, GPU_THREADS_PER_BLOCK>>>(npoints_grid,d,ndims,nvars,ghosts,stride,dim,recvbuf,var);
    }
  }
  cudaDeviceSynchronize();
  /* Wait till send requests are complete before freeing memory */
  MPI_Waitall(2*ndims,sndreq,status_arr);
#endif
  return 0;
}

#else

/*! Kernel function for gpuMPIExchangeBoundariesnD() */
__global__
void gpuFillSendBufLeft_kernel(
  int npoints_grid,
  int npoints_local_wghosts,
  int d,
  int ndims,
  int nvars,
  int ghosts,
  int stride,
  const int    * __restrict__ dim,
  const double * __restrict__ var,
  double       * __restrict__ sendbuf
)
{
  int p = blockDim.x * blockIdx.x + threadIdx.x;

  if (p < npoints_grid) {
    int index[GPU_MAX_NDIMS], bounds[GPU_MAX_NDIMS], offset[GPU_MAX_NDIMS];
    int p1, p2;

    _ArrayCopy1D_(dim,bounds,ndims); bounds[d] = ghosts;
    _ArraySetValue_(offset,ndims,0);
    _ArrayIndexnD_(ndims,p,bounds,index,0);
    _ArrayIndex1DWO_(ndims,dim,index,offset,ghosts,p1);
    _ArrayIndex1D_(ndims,bounds,index,0,p2);
    for (int v=0; v<nvars; v++) {
      sendbuf[2*d*stride+nvars*p2+v] = var[p1+v*npoints_local_wghosts];
    }
    //_ArrayCopy1D_((var+nvars*p1),(sendbuf+2*d*stride+nvars*p2),nvars);
  }

  return;
}

/*! Kernel function for gpuMPIExchangeBoundariesnD() */
__global__
void gpuFillSendBufRight_kernel(
  int npoints_grid,
  int npoints_local_wghosts,
  int d,
  int ndims,
  int nvars,
  int ghosts,
  int stride,
  const int    * __restrict__ dim,
  const double * __restrict__ var,
  double       * __restrict__ sendbuf
)
{
  int p = blockDim.x * blockIdx.x + threadIdx.x;

  if (p < npoints_grid) {
    int index[GPU_MAX_NDIMS], bounds[GPU_MAX_NDIMS], offset[GPU_MAX_NDIMS];
    int p1, p2;

    _ArrayCopy1D_(dim,bounds,ndims); bounds[d] = ghosts;
    _ArraySetValue_(offset,ndims,0); offset[d] = dim[d]-ghosts;
    _ArrayIndexnD_(ndims,p,bounds,index,0);
    _ArrayIndex1DWO_(ndims,dim,index,offset,ghosts,p1);
    _ArrayIndex1D_(ndims,bounds,index,0,p2);
    for (int v=0; v<nvars; v++) {
      sendbuf[(2*d+1)*stride+nvars*p2+v] = var[p1+v*npoints_local_wghosts];
    }
    //_ArrayCopy1D_((var+nvars*p1),(sendbuf+(2*d+1)*stride+nvars*p2),nvars);
  }

  return;
}

/*! Kernel function for gpuMPIExchangeBoundariesnD() */
__global__
void gpuFillRecvBufLeft_kernel(
  int npoints_grid,
  int npoints_local_wghosts,
  int d,
  int ndims,
  int nvars,
  int ghosts,
  int stride,
  const int    * __restrict__ dim,
  const double * __restrict__ recvbuf,
  double * __restrict__ var
)
{
  int p = blockDim.x * blockIdx.x + threadIdx.x;

  if (p < npoints_grid) {
    int index[GPU_MAX_NDIMS], bounds[GPU_MAX_NDIMS], offset[GPU_MAX_NDIMS];
    int p1, p2;

    _ArrayCopy1D_(dim,bounds,ndims); bounds[d] = ghosts;
    _ArraySetValue_(offset,ndims,0); offset[d] = -ghosts;
    _ArrayIndexnD_(ndims,p,bounds,index,0);
    _ArrayIndex1DWO_(ndims,dim,index,offset,ghosts,p1);
    _ArrayIndex1D_(ndims,bounds,index,0,p2);
    for (int v=0; v<nvars; v++) {
      var[p1+v*npoints_local_wghosts] = recvbuf[2*d*stride+nvars*p2+v];
    }
    //_ArrayCopy1D_((recvbuf+2*d*stride+nvars*p2),(var+nvars*p1),nvars);
  }

  return;
}

/*! Kernel function for gpuMPIExchangeBoundariesnD() */
__global__
void gpuFillRecvBufRight_kernel(
  int npoints_grid,
  int npoints_local_wghosts,
  int d,
  int ndims,
  int nvars,
  int ghosts,
  int stride,
  const int    * __restrict__ dim,
  const double * __restrict__ recvbuf,
  double       * __restrict__ var
)
{
  int p = blockDim.x * blockIdx.x + threadIdx.x;

  if (p < npoints_grid) {
    int index[GPU_MAX_NDIMS], bounds[GPU_MAX_NDIMS], offset[GPU_MAX_NDIMS];
    int p1, p2;

    _ArrayCopy1D_(dim,bounds,ndims); bounds[d] = ghosts;
    _ArraySetValue_(offset,ndims,0); offset[d] = dim[d];
    _ArrayIndexnD_(ndims,p,bounds,index,0);
    _ArrayIndex1DWO_(ndims,dim,index,offset,ghosts,p1);
    _ArrayIndex1D_(ndims,bounds,index,0,p2);
    for (int v=0; v<nvars; v++) {
      var[p1+v*npoints_local_wghosts] = recvbuf[(2*d+1)*stride+nvars*p2+v];
    }
    //_ArrayCopy1D_((recvbuf+(2*d+1)*stride+nvars*p2),(var+nvars*p1),nvars);
  }

  return;
}

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
extern "C" int gpuMPIExchangeBoundariesnD(
  int    ndims,
  int    nvars,
  const int * __restrict__ dim,
  int    ghosts,
  void   *m,
  double    * __restrict__ var
)
{
#ifndef serial

  MPIVariables  *mpi = (MPIVariables*) m;
  int           d;

  int *ip      = mpi->ip;
  int *iproc   = mpi->iproc;
  int *bcflag  = mpi->bcperiodic;
  int *cpu_dim = mpi->cpu_dim;
  int size = 1; for (d=0; d<ndims; d++) size *= (cpu_dim[d]+2*ghosts);

  int neighbor_rank[2*ndims], nip[ndims], bounds[ndims];
  MPI_Request rcvreq[2*ndims], sndreq[2*ndims];
  for (d=0; d<2*ndims; d++) rcvreq[d] = sndreq[d] = MPI_REQUEST_NULL;

  /* each process has 2*ndims neighbors (except at non-periodic physical boundaries)  */
  /* calculate the rank of these neighbors (-1 -> none)                               */
  for (d = 0; d < ndims; d++) {
    _ArrayCopy1D_(ip,nip,ndims);
    if (ip[d] == 0) nip[d] = iproc[d]-1;
    else            nip[d]--;
    if ((ip[d] == 0) && (!bcflag[d])) neighbor_rank[2*d]   = -1;
    else                              neighbor_rank[2*d]   = MPIRank1D(ndims,iproc,nip);
    _ArrayCopy1D_(ip,nip,ndims);
    if (ip[d] == (iproc[d]-1)) nip[d] = 0;
    else                       nip[d]++;
    if ((ip[d] == (iproc[d]-1)) && (!bcflag[d]))  neighbor_rank[2*d+1] = -1;
    else                                          neighbor_rank[2*d+1] = MPIRank1D(ndims,iproc,nip);
  }

  /* calculate dimensions of each of the send-receive regions */
  double *sendbuf = mpi->gpu_sendbuf;
  double *recvbuf = mpi->gpu_recvbuf;
  int    stride   = mpi->maxbuf;
  int    bufdim[ndims];
  for (d = 0; d < ndims; d++) {
    bufdim[d] = 1;
    int i;
    for (i = 0; i < ndims; i++) {
      if (i == d) bufdim[d] *= ghosts;
      else        bufdim[d] *= cpu_dim[i];
    }
  }

  /* post the receive requests */
  for (d = 0; d < ndims; d++) {
    if (neighbor_rank[2*d  ] != -1) {
      MPI_Irecv(&recvbuf[2*d*stride],bufdim[d]*nvars,MPI_DOUBLE,neighbor_rank[2*d  ],1630,
                mpi->world,&rcvreq[2*d]);
    }
    if (neighbor_rank[2*d+1] != -1) {
      MPI_Irecv(&recvbuf[(2*d+1)*stride],bufdim[d]*nvars,MPI_DOUBLE,neighbor_rank[2*d+1],1631,
                mpi->world,&rcvreq[2*d+1]);
    }
  }

  /* count number of neighbors and copy data to send buffers */
  for (d = 0; d < ndims; d++) {
    _ArrayCopy1D_(cpu_dim,bounds,ndims); bounds[d] = ghosts;
    if (neighbor_rank[2*d] != -1) {
      int npoints_grid = 1; for (int i=0; i<ndims; i++) npoints_grid *= bounds[i];
      int nblocks = (npoints_grid-1)/GPU_THREADS_PER_BLOCK + 1;
      gpuFillSendBufLeft_kernel<<<nblocks, GPU_THREADS_PER_BLOCK>>>(
        npoints_grid,size,d,ndims,nvars,ghosts,stride,dim,var,sendbuf);
    }
    if (neighbor_rank[2*d+1] != -1) {
      int npoints_grid = 1; for (int i=0; i<ndims; i++) npoints_grid *= bounds[i];
      int nblocks = (npoints_grid-1)/GPU_THREADS_PER_BLOCK + 1;
      gpuFillSendBufRight_kernel<<<nblocks, GPU_THREADS_PER_BLOCK>>>(
        npoints_grid,size,d,ndims,nvars,ghosts,stride,dim,var,sendbuf);
    }
  }
  cudaDeviceSynchronize();

  /* send the data */
  for (d = 0; d < ndims; d++) {
    if (neighbor_rank[2*d  ] != -1) {
      MPI_Isend(&sendbuf[2*d*stride],bufdim[d]*nvars,MPI_DOUBLE,neighbor_rank[2*d  ],1631,
                mpi->world,&sndreq[2*d]);
    }
    if (neighbor_rank[2*d+1] != -1) {
      MPI_Isend(&sendbuf[(2*d+1)*stride],bufdim[d]*nvars,MPI_DOUBLE,neighbor_rank[2*d+1],1630,
                mpi->world,&sndreq[2*d+1]);
    }
  }

  /* Wait till data is done received */
  MPI_Waitall(2*ndims,rcvreq,MPI_STATUS_IGNORE);

  /* copy received data to ghost points */
  for (d = 0; d < ndims; d++) {
    _ArrayCopy1D_(cpu_dim,bounds,ndims); bounds[d] = ghosts;
    if (neighbor_rank[2*d] != -1) {
      int npoints_grid = 1; for (int i=0; i<ndims; i++) npoints_grid *= bounds[i];
      int nblocks = (npoints_grid-1)/GPU_THREADS_PER_BLOCK + 1;
      gpuFillRecvBufLeft_kernel<<<nblocks, GPU_THREADS_PER_BLOCK>>>(
        npoints_grid,size,d,ndims,nvars,ghosts,stride,dim,recvbuf,var);
    }
    if (neighbor_rank[2*d+1] != -1) {
      int npoints_grid = 1; for (int i=0; i<ndims; i++) npoints_grid *= bounds[i];
      int nblocks = (npoints_grid-1)/GPU_THREADS_PER_BLOCK + 1;
      gpuFillRecvBufRight_kernel<<<nblocks, GPU_THREADS_PER_BLOCK>>>(
        npoints_grid,size,d,ndims,nvars,ghosts,stride,dim,recvbuf,var);
    }
  }
  cudaDeviceSynchronize();

  /* Wait till send requests are complete before freeing memory */
  MPI_Waitall(2*ndims,sndreq,MPI_STATUS_IGNORE);

#endif
  return 0;
}

#endif

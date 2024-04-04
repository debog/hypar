/*! @file SparseGridsInterpolate.cpp
    @brief Functions to interpolate from one grid to another
    @author Debojyoti Ghosh, John Loffeld, Lee Ricketson
*/

#include <cmath>
#include <basic.h>
#include <arrayfunctions.h>
#include <std_vec_ops.h>
#include <sparse_grids_simulation.h>

/*! Interpolate grid coordinates from one grid to another. Note that along each dimension,
    the ratio of the number of grid points in the source grid and that in the
    destination grid must be an integer power of 2 (negative or positive).

    Since the source and destination grids may have different processor layouts,
    the interpolation is carried out on rank 0 by allocating global arrays on it.
    Therefore, this function is *not* scalable.
*/
void SparseGridsSimulation::interpolateGrid ( SimulationObject* const       a_dst,  /*!< Destination object */
                                              const SimulationObject* const a_src   /*!< Source object */
                                            )
{
  /* get the source and destination grid dimensions */
  GridDimensions dim_src, dim_dst;
  StdVecOps::copyFrom(dim_src, a_src->solver.dim_global, m_ndims);
  StdVecOps::copyFrom(dim_dst, a_dst->solver.dim_global, m_ndims);

  /* gather the source on rank 0 */
  double* xg_src = NULL;
  if (!m_rank) allocateGridArrays(dim_src, &xg_src);
  {
    int offset_global, offset_local;
    offset_global = offset_local = 0;
    for (int d=0; d<m_ndims; d++) {
      MPIGatherArray1D(  (void*) &(a_src->mpi),
                         (a_src->mpi.rank ? NULL : &xg_src[offset_global]),
                         &(a_src->solver.x[offset_local+a_src->solver.ghosts]),
                         a_src->mpi.is[d],
                         a_src->mpi.ie[d],
                         a_src->solver.dim_local[d],
                         0 );
      offset_global += a_src->solver.dim_global[d];
      offset_local  += a_src->solver.dim_local [d] + 2*a_src->solver.ghosts;
    }
  }

  /* now do the interpolation, dimension-by-dimension */
  double *xg_dst = NULL;
  if (!m_rank) {

    GridDimensions dim_to(m_ndims,0);
    GridDimensions dim_from(m_ndims,0);

    double* x_from;
    double* x_to;

    dim_to = dim_src;
    x_to = xg_src;
    x_from = NULL;

    for (int dir = 0; dir < m_ndims; dir++) {

      dim_from = dim_to;
      dim_to[dir] = dim_dst[dir];

      if (dim_from[dir] == dim_to[dir]) continue;

      double fac = (dim_to[dir] > dim_from[dir] ?
                        (double)dim_to[dir]/(double)dim_from[dir]
                      : (double)dim_from[dir]/(double)dim_to[dir] );
      if (!isPowerOfTwo((int)fac)) {
        fprintf(stderr,"Error in SparseGridsSimulation::interpolate() - \n");
        fprintf(stderr,"  refinement/coarsening factor not a power of 2!\n");
        exit(1);
      }

      if (x_from != NULL) free(x_from);
      x_from = x_to;

      allocateGridArrays(dim_to, &x_to);
      if (dim_to[dir] < dim_from[dir]) {
        coarsenGrid1D(dim_from, dim_to, x_from, x_to, dir);
      } else {
        refineGrid1D(dim_from, dim_to, x_from, x_to, dir);
      }

    }

    /* dim_to should be equal to dim_dst now */
    for (int d = 0; d < m_ndims; d++) {
      if (dim_to[d] != dim_dst[d]) {
        fprintf(stderr,"Error in SparseGridsSimulation::interpolate() - \n");
        fprintf(stderr,"  dim_to[%d] != dim_dst[%d]!\n", d, d);
        exit(1);
      }
    }

    if (x_from != NULL) free(x_from);
    xg_dst = x_to;

  }

  /* partition the destination */
  int offset_global = 0, offset_local = 0;
  for (int d=0; d<m_ndims; d++) {
    MPIPartitionArray1D(  (void*) &(a_dst->mpi),
                          (a_dst->mpi.rank ? NULL : &xg_dst[offset_global]),
                          &(a_dst->solver.x[offset_local+a_dst->solver.ghosts]),
                          a_dst->mpi.is[d],
                          a_dst->mpi.ie[d],
                          a_dst->solver.dim_local[d],
                          0);
    offset_global += a_dst->solver.dim_global[d];
    offset_local  += a_dst->solver.dim_local [d] + 2*a_dst->solver.ghosts;
  }

  if (!m_rank) {
    free(xg_dst);
  }

  /* exchange MPI-boundary values of x between processors */
  {
    int offset = 0;
    for (int d = 0; d < m_ndims; d++) {
      MPIExchangeBoundaries1D(  (void*) &(a_dst->mpi),
                                &(a_dst->solver.x[offset]),
                                a_dst->solver.dim_local[d],
                                a_dst->solver.ghosts,
                                d,
                                m_ndims);
      offset += (a_dst->solver.dim_local[d] + 2*a_dst->solver.ghosts);
    }
  }
  /* fill in ghost values of x at physical boundaries by extrapolation */
  {
    int offset = 0;
    for (int d = 0; d < m_ndims; d++) {
      double* X       = &(a_dst->solver.x[offset]);
      int*    dim     = a_dst->solver.dim_local;
      int     ghosts  = a_dst->solver.ghosts;
      if (a_dst->mpi.ip[d] == 0) {
        /* fill left boundary along this dimension */
        for (int i = 0; i < ghosts; i++) {
          int delta = ghosts - i;
          X[i] = X[ghosts] + ((double) delta) * (X[ghosts]-X[ghosts+1]);
        }
      }
      if (a_dst->mpi.ip[d] == a_dst->mpi.iproc[d]-1) {
        /* fill right boundary along this dimension */
        for (int i = dim[d]+ghosts; i < dim[d]+2*ghosts; i++) {
          int delta = i - (dim[d]+ghosts-1);
          X[i] =  X[dim[d]+ghosts-1]
                  + ((double) delta) * (X[dim[d]+ghosts-1]-X[dim[d]+ghosts-2]);
        }
      }
      offset  += (dim[d] + 2*ghosts);
    }
  }
  return;
}

/*! Coarsen along a given dimension - the source and destination must have the same sizes
 *  along all other dimensions.
 *
 *  Note that the grid must *not* have any ghost points!
*/
void SparseGridsSimulation::coarsenGrid1D(  const GridDimensions& a_dim_src,  /*!< Grid size of source data */
                                            const GridDimensions& a_dim_dst,  /*!< Grid size of destination data */
                                            const double* const   a_x_src,    /*!< Source grid */
                                            double* const         a_x_dst,    /*!< Destination grid */
                                            const int             a_dir       /*!< Dimension along which to coarsen */
                                         )
{
  for (int d = 0; d < m_ndims; d++) {
    if ((d != a_dir) && (a_dim_src[d] != a_dim_dst[d])) {
      fprintf(stderr, "Error in SparseGridsSimulation::coarsenGrid1D() -\n");
      fprintf(stderr, " a_dim_src[%d] != a_dim_dst[%d]\n", d, d);
      exit(1);
    }
  }

  int n_src = a_dim_src[a_dir];
  int n_dst = a_dim_dst[a_dir];
  if (n_dst > n_src) {
    fprintf(stderr, "Error in SparseGridsSimulation::coarsenGrid1D() -\n");
    fprintf(stderr, " destination grid is finer than source grid along a_dir!\n");
    exit(1);
  }

  double fac = ((double) n_src) / ((double) n_dst);
  int stride = (int) fac;
  if (std::abs(((double)stride)-fac) > _MACHINE_ZERO_) {
    fprintf(stderr, "Error in SparseGridsSimulation::coarsenGrid1D() -\n");
    fprintf(stderr, "  non-integer coarsening factor!\n");
    exit(1);
  }

  const double* x_src = a_x_src;
  double* x_dst = a_x_dst;

  for (int d = 0; d < m_ndims; d++) {

    if (d == a_dir) {

      for (int i = 0; i < a_dim_dst[d]; i++) {
        double avg = 0;
        for (int j=stride*i; j<stride*(i+1); j++) {
          if (j >= n_src) {
            fprintf(stderr, "Error in SparseGridsSimulation::coarsenGrid1D() -\n");
            fprintf(stderr, "  j >= n_src!\n");
            exit(1);
          }
          avg += x_src[j];
        }
        avg /= (double) stride;
        x_dst[i] = avg;
      }

    } else {

      _ArrayCopy1D_(x_src, x_dst, a_dim_dst[d]);

    }

    x_src += a_dim_src[d];
    x_dst += a_dim_dst[d];
  }

  return;
}

/*! Refine along a given dimension - the source and destination must have the same sizes
 *  along all other dimensions.
 *
 *  Note that the grid must *not* have any ghost points!
*/
void SparseGridsSimulation::refineGrid1D( const GridDimensions& a_dim_src,  /*!< Grid size of source data */
                                          const GridDimensions& a_dim_dst,  /*!< Grid size of destination data */
                                          const double* const   a_x_src,    /*!< Source grid */
                                          double* const         a_x_dst,    /*!< Destination grid */
                                          const int             a_dir       /*!< Dimension along which to refine */
                                        )
{
  for (int d = 0; d < m_ndims; d++) {
    if ((d != a_dir) && (a_dim_src[d] != a_dim_dst[d])) {
      fprintf(stderr, "Error in SparseGridsSimulation::coarsenGrid1D() -\n");
      fprintf(stderr, " a_dim_src[%d] != a_dim_dst[%d]\n", d, d);
      exit(1);
    }
  }

  int n_src = a_dim_src[a_dir];
  int n_dst = a_dim_dst[a_dir];
  if (n_dst < n_src) {
    fprintf(stderr, "Error in SparseGridsSimulation::refineGrid1D() -\n");
    fprintf(stderr, "  destination grid is coarser than source grid along a_dir!\n");
    exit(1);
  }

  double fac = ((double) n_dst) / ((double) n_src);
  int stride = (int) fac;
  if (std::abs(((double)stride)-fac) > _MACHINE_ZERO_) {
    fprintf(stderr, "Error in SparseGridsSimulation::refineGrid1D() -\n");
    fprintf(stderr, "  non-integer refinement factor!\n");
    exit(1);
  }

  int offset = 0;
  for (int d = 0; d < a_dir; d++) offset += a_dim_src[d];
  const double* x_src = a_x_src + offset;
  double* x_dst = a_x_dst + offset;

  fprintf(stderr, "Error in SparseGridsSimulation::refineGrid1D() -\n");
  fprintf(stderr, "  NOT YET IMPLEMENTED! Why do you need this?\n");
  exit(1);

  return;
}


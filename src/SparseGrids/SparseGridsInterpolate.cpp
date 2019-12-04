/*! @file SparseGridsInterpolate.cpp
    @brief Functions to interpolate from one grid to another
    @author Debojyoti Ghosh, John Loffeld, Lee Ricketson
*/

#include <arrayfunctions.h>
#include <std_vec_ops.h>
#include <sparse_grids_simulation.h>

/*! Interpolate data from one grid to another. Note that along each dimension,
    the ratio of the number of grid points in the source grid and that in the 
    destination grid must be an integer power of 2 (negative or positive).
*/
void SparseGridsSimulation::interpolate(  SimulationObject* const       a_dst,  /*!< Destination object */
                                          const SimulationObject* const a_src   /*!< Source object */
                                      )
{
  /* get the source and destination grid dimensions */
  GridDimensions dim_dst;
  StdVecOps::copyFrom(dim_dst, a_dst->solver.dim_global, m_ndims);

  /* get number of vector components */
  int nvars = a_src->solver.nvars;
  if (nvars != a_dst->solver.nvars) {
    fprintf(stderr, "Error in SparseGridsSimulation::interpolate(): unequal nvars\n");
    exit(1);
  }

  double *ug_dst = NULL;
  interpolate(dim_dst, &ug_dst, a_src); 

  /* partition the destination */
  MPIPartitionArraynD( m_ndims,
                       (void*) &(a_dst->mpi),
                       (a_dst->mpi.rank ? NULL : ug_dst),
                       a_dst->solver.u,
                       a_dst->solver.dim_global,
                       a_dst->solver.dim_local,
                       a_dst->solver.ghosts,
                       nvars); 

  if (!m_rank) {
    delete[] ug_dst;
  }

  return;
}

/*! Interpolate data from one grid to another of a desired resolution. 
    Note that along each dimension, the ratio of the number of grid points 
    in the source grid and that in the destination grid must be an integer 
    power of 2 (negative or positive).

    The incoming pointer must be NULL. After this function is executed, it
    will point to a chunk of memory with the interpolated solution. This is
    the *global* solution.

    Since the source and destination grids may have different processor layouts,
    the interpolation is carried out on rank 0 by allocating global arrays on it.
    Therefore, this function is *not* scalable.
*/
void SparseGridsSimulation::interpolate(const GridDimensions& a_dim_dst, /*!< grid dimensions to interpolate to */
                                        double** const a_u_dst, /*!< pointer to array containing interpolated data */
                                        const SimulationObject* const a_src /*!< Source object */
                                       )
{
  if ((*a_u_dst) != NULL) {
    fprintf(stderr, "Error in SparseGridsSimulation::interpolate() - \n");
    fprintf(stderr, "  a_u_dst is not NULL!\n");
    exit(1);
  }

  /* get the source grid dimensions */
  GridDimensions dim_src;
  StdVecOps::copyFrom(dim_src, a_src->solver.dim_global, m_ndims);

  /* get number of vector components */
  int nvars = a_src->solver.nvars;

  /* gather the source on rank 0 */
  double *ug_src = NULL;
  if (!m_rank) allocateDataArrays(dim_src, nvars, &ug_src);
  MPIGatherArraynD(  m_ndims,
                     (void*) &(a_src->mpi),
                     ug_src,
                     a_src->solver.u,
                     a_src->solver.dim_global,
                     a_src->solver.dim_local,
                     a_src->solver.ghosts,
                     nvars);
  

  /* now do the interpolation, dimension-by-dimension */
  if (!m_rank) {

    GridDimensions dim_to(m_ndims,0);
    GridDimensions dim_from(m_ndims,0);

    double *u_from;
    double *u_to;

    dim_to = dim_src;
    u_to = ug_src;
    u_from = NULL;

    for (int dir = 0; dir < m_ndims; dir++) {

      dim_from = dim_to;
      dim_to[dir] = a_dim_dst[dir];

      if (dim_from[dir] == dim_to[dir]) continue;

      double fac = (dim_to[dir] > dim_from[dir] ? 
                        (double)dim_to[dir]/(double)dim_from[dir] 
                      : (double)dim_from[dir]/(double)dim_to[dir] );
      if (!isPowerOfTwo((int)fac)) {
        fprintf(stderr,"Error in SparseGridsSimulation::interpolate() - \n");
        fprintf(stderr,"  refinement/coarsening factor not a power of 2!\n");
        exit(1);
      }

      if (u_from != NULL) delete[] u_from;

      u_from = u_to;

      allocateDataArrays(dim_to, nvars, &u_to);
      if (dim_to[dir] < dim_from[dir]) {
        coarsen1D(dim_from, dim_to, u_from, u_to, dir, nvars);
      } else {
        refine1D(dim_from, dim_to, u_from, u_to, dir, nvars);
      }
  
    }

    /* dim_to should be equal to a_dim_dst now */
    for (int d = 0; d < m_ndims; d++) {
      if (dim_to[d] != a_dim_dst[d]) {
        fprintf(stderr,"Error in SparseGridsSimulation::interpolate() - \n");
        fprintf(stderr,"  dim_to[%d] != a_dim_dst[%d]!\n", d, d);
        exit(1);
      }
    }

    if (u_from != NULL) delete[] u_from;
    (*a_u_dst) = u_to;

  }

  return;
}

/*! Coarsen along a given dimension - the source and destination must have the same sizes
 *  along all other dimensions. 
 *
 *  Note that the solution must *not* have any ghost points!
*/
void SparseGridsSimulation::coarsen1D(  const GridDimensions& a_dim_src,  /*!< Grid size of source data */
                                        const GridDimensions& a_dim_dst,  /*!< Grid size of destination data */
                                        const double* const   a_u_src,    /*!< Source solution */
                                        double* const         a_u_dst,    /*!< Destination solution */
                                        const int             a_dir,      /*!< Dimension along which to coarsen */
                                        const int             a_nvars     /*!< Number of vector components of the solution */
                                     )
{
  for (int d = 0; d < m_ndims; d++) {
    if ((d != a_dir) && (a_dim_src[d] != a_dim_dst[d])) {
      fprintf(stderr, "Error in SparseGridsSimulation::coarsen1D() -\n");
      fprintf(stderr, " a_dim_src[%d] != a_dim_dst[%d]\n", d, d);
      exit(1);
    }
  }

  int n_src = a_dim_src[a_dir];
  int n_dst = a_dim_dst[a_dir];
  if (n_dst > n_src) {
    fprintf(stderr, "Error in SparseGridsSimulation::coarsen1D() -\n");
    fprintf(stderr, " destination grid is finer than source grid along a_dir!\n");
    exit(1);
  }

  double fac = ((double) n_src) / ((double) n_dst);
  int stride = (int) fac;
  if (std::abs(((double)stride)-fac) > _MACHINE_ZERO_) {
    fprintf(stderr, "Error in SparseGridsSimulation::coarsen1D() -\n");
    fprintf(stderr, "  non-integer coarsening factor!\n");
    exit(1);
  }

  /* create bounds for the transverse loop, i.e., to loop over 
   * all 1D lines along dimension "a_dir" */
  int bounds_transverse[m_ndims];
  _ArrayCopy1D_(a_dim_src, bounds_transverse, m_ndims); 
  bounds_transverse[a_dir] =  1;
  int n_transverse; 
  _ArrayProduct1D_(bounds_transverse, m_ndims, n_transverse);

  for (int ti = 0; ti < n_transverse; ti++) {

    int index_transverse[m_ndims];
    _ArrayIndexnD_(m_ndims, ti, bounds_transverse, index_transverse, 0);

    int index_dst[m_ndims], index_src[m_ndims];
    _ArrayCopy1D_(index_transverse,index_dst,m_ndims);
    _ArrayCopy1D_(index_transverse,index_src,m_ndims);

    for (index_dst[a_dir] = 0; index_dst[a_dir] < n_dst; index_dst[a_dir]++) {

      for (int v = 0; v < a_nvars; v++) {

        double avg = 0;
        for ( index_src[a_dir] = index_dst[a_dir]*stride;
              index_src[a_dir] < (index_dst[a_dir]+1)*stride;
              index_src[a_dir]++ ) {
          int p;
          _ArrayIndex1D_(m_ndims, a_dim_src, index_src, 0, p);
          avg += a_u_src[p*a_nvars+v];
        }
        avg /= (double) stride;

        {
          int p;
          _ArrayIndex1D_(m_ndims, a_dim_dst, index_dst, 0, p);
          a_u_dst[p*a_nvars+v] = avg;
        }


      }

    }

  }

  return;
}

/*! Refine along a given dimension - the source and destination must have the same sizes
 *  along all other dimensions. 
 *
 *  Note that the solution must *not* have any ghost points!
*/
void SparseGridsSimulation::refine1D( const GridDimensions& a_dim_src,  /*!< Grid size of source data */
                                      const GridDimensions& a_dim_dst,  /*!< Grid size of destination data */
                                      const double* const   a_u_src,    /*!< Source solution */
                                      double* const         a_u_dst,    /*!< Destination solution */
                                      const int             a_dir,      /*!< Dimension along which to coarsen */
                                      const int             a_nvars     /*!< Number of vector components of the solution */
                                    )
{
  for (int d = 0; d < m_ndims; d++) {
    if ((d != a_dir) && (a_dim_src[d] != a_dim_dst[d])) {
      fprintf(stderr, "Error in SparseGridsSimulation::coarsen1D() -\n");
      fprintf(stderr, " a_dim_src[%d] != a_dim_dst[%d]\n", d, d);
      exit(1);
    }
  }

  int n_src = a_dim_src[a_dir];
  int n_dst = a_dim_dst[a_dir];
  if (n_dst < n_src) {
    fprintf(stderr, "Error in SparseGridsSimulation::refine1D() -\n");
    fprintf(stderr, "  destination grid is coarser than source grid along a_dir!\n");
    exit(1);
  }

  double fac = ((double) n_dst) / ((double) n_src);
  int stride = (int) fac;
  if (std::abs(((double)stride)-fac) > _MACHINE_ZERO_) {
    fprintf(stderr, "Error in SparseGridsSimulation::refine1D() -\n");
    fprintf(stderr, "  non-integer refinement factor!\n");
    exit(1);
  }

  fprintf(stderr, "Error in SparseGridsSimulation::refine1D() -\n");
  fprintf(stderr, "  NOT YET IMPLEMENTED\n");
  exit(1);
  
  return;
}


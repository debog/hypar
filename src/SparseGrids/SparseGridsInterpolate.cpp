/*! @file SparseGridsInterpolate.cpp
    @brief Functions to interpolate from one grid to another
    @author Debojyoti Ghosh, John Loffeld, Lee Ricketson
*/

#include <cmath>
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
  /* get the destination grid dimensions */
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
  MPIPartitionArraynDwGhosts( m_ndims,
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

  /* get number of vector components and number of ghost points*/
  int nvars = a_src->solver.nvars;
  int ghosts = a_src->solver.ghosts;

  /* gather the source on rank 0 */
  double *ug_src = NULL;
  if (!m_rank) allocateDataArrays(dim_src, nvars, &ug_src, ghosts);
  MPIGatherArraynDwGhosts( m_ndims,
                           (void*) &(a_src->mpi),
                           ug_src,
                           a_src->solver.u,
                           a_src->solver.dim_global,
                           a_src->solver.dim_local,
                           ghosts,
                           nvars );

  interpolate(  a_dim_dst,
                a_u_dst,
                dim_src,
                ug_src,
                nvars,
                ghosts );
  return;
}

/*! Interpolate data from one grid to another of a desired resolution. 
    Note that along each dimension, the ratio of the number of grid points 
    in the source grid and that in the destination grid must be an integer 
    power of 2 (negative or positive).

    The source data *must* be global. It will get deallocated at the end
    of this function.

    The incoming pointer must be NULL. After this function is executed, it
    will point to a chunk of memory with the interpolated solution. This is
    the *global* solution.
*/
void SparseGridsSimulation::interpolate(const GridDimensions& a_dim_dst, /*!< grid dimensions to interpolate to */
                                        double** const a_u_dst, /*!< pointer to array containing interpolated data */
                                        const GridDimensions& a_dim_src, /*!< grid dimensions to interpolate from */
                                        double* const a_u_src, /*!< pointer to array containing data to interpolate from */
                                        const int a_nvars, /*!< Number of vector components of the solution */
                                        const int a_ghosts /*!< Number of ghost points */
                                       )
{
  if ((*a_u_dst) != NULL) {
    fprintf(stderr, "Error in SparseGridsSimulation::interpolate() - \n");
    fprintf(stderr, "  a_u_dst is not NULL!\n");
    exit(1);
  }

  /* now do the interpolation, dimension-by-dimension */
  if (!m_rank) {

    if (a_u_src == NULL) {
      fprintf(stderr, "Error in SparseGridsSimulation::interpolate() - \n");
      fprintf(stderr, "  a_u_src is NULL!\n");
      exit(1);
    }

    GridDimensions dim_to(m_ndims,0);
    GridDimensions dim_from(m_ndims,0);

    double *u_from;
    double *u_to;

    dim_to = a_dim_src;
    u_to = a_u_src;
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

      allocateDataArrays(dim_to, a_nvars, &u_to, a_ghosts);
      if (dim_to[dir] < dim_from[dir]) {
        coarsen1D(dim_from, dim_to, u_from, u_to, dir, a_nvars, a_ghosts);
      } else {
        refine1D(dim_from, dim_to, u_from, u_to, dir, a_nvars, a_ghosts);
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
                                        const int             a_nvars,    /*!< Number of vector components of the solution */
                                        const int             a_ngpt      /*!< Number of ghost points */
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

  /* compute dimensions with ghosts */
  GridDimensions dim_src_wg = a_dim_src;
  StdVecOps::add(dim_src_wg, 2*a_ngpt);
  GridDimensions dim_dst_wg = a_dim_dst;
  StdVecOps::add(dim_dst_wg, 2*a_ngpt);

  /* create bounds for the transverse loop, i.e., to loop over 
   * all 1D lines along dimension "a_dir" */
  int bounds_transverse[m_ndims];
  _ArrayCopy1D_(a_dim_src, bounds_transverse, m_ndims); 
  for (int d = 0; d < m_ndims; d++) bounds_transverse[d] += (2*a_ngpt);
  bounds_transverse[a_dir] =  1;

  int n_transverse; 
  _ArrayProduct1D_(bounds_transverse, m_ndims, n_transverse);

  for (int ti = 0; ti < n_transverse; ti++) {

    int index_transverse[m_ndims];
    _ArrayIndexnD_(m_ndims, ti, bounds_transverse, index_transverse, 0);

    int index_dst[m_ndims], index_src[m_ndims];
    _ArrayCopy1D_(index_transverse, index_dst, m_ndims);
    _ArrayCopy1D_(index_transverse, index_src, m_ndims);

    for (int i_dst = 0; i_dst < n_dst; i_dst++) {

      int i_m1 = i_dst*stride + (stride/2-1);
      int i_m2 = i_m1 - 1;
      int i_p1 = i_m1 + 1;
      int i_p2 = i_p1 + 1;

      int p;
      index_dst[a_dir] = i_dst+a_ngpt;
      _ArrayIndex1D_(m_ndims, dim_dst_wg, index_dst, 0, p);

      int p_m2;
      index_src[a_dir] = i_m2+a_ngpt;
      _ArrayIndex1D_(m_ndims, dim_src_wg, index_src, 0, p_m2);

      int p_m1;
      index_src[a_dir] = i_m1+a_ngpt;
      _ArrayIndex1D_(m_ndims, dim_src_wg, index_src, 0, p_m1);

      int p_p1;
      index_src[a_dir] = i_p1+a_ngpt;
      _ArrayIndex1D_(m_ndims, dim_src_wg, index_src, 0, p_p1);

      int p_p2;
      index_src[a_dir] = i_p2+a_ngpt;
      _ArrayIndex1D_(m_ndims, dim_src_wg, index_src, 0, p_p2);

      for (int v = 0; v < a_nvars; v++) {
        double val =  - ( 1.0/16.0) * a_u_src[p_m2*a_nvars+v]
                      + ( 9.0/16.0) * a_u_src[p_m1*a_nvars+v] 
                      + ( 9.0/16.0) * a_u_src[p_p1*a_nvars+v] 
                      - ( 1.0/16.0) * a_u_src[p_p2*a_nvars+v];
        a_u_dst[p*a_nvars+v] = val;
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
                                      const int             a_nvars,    /*!< Number of vector components of the solution */
                                      const int             a_ngpt      /*!< Number of ghost points */
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

  /* compute dimensions with ghosts */
  GridDimensions dim_src_wg = a_dim_src;
  StdVecOps::add(dim_src_wg, 2*a_ngpt);
  GridDimensions dim_dst_wg = a_dim_dst;
  StdVecOps::add(dim_dst_wg, 2*a_ngpt);

  /* create bounds for the transverse loop, i.e., to loop over 
   * all 1D lines along dimension "a_dir" */
  int bounds_transverse[m_ndims];
  _ArrayCopy1D_(a_dim_src, bounds_transverse, m_ndims); 
  for (int d = 0; d < m_ndims; d++) bounds_transverse[d] += (2*a_ngpt);
  bounds_transverse[a_dir] =  1;

  int n_transverse; 
  _ArrayProduct1D_(bounds_transverse, m_ndims, n_transverse);

  for (int ti = 0; ti < n_transverse; ti++) {

    int index_transverse[m_ndims];
    _ArrayIndexnD_(m_ndims, ti, bounds_transverse, index_transverse, 0);

    int index_dst [m_ndims], 
        index_src0[m_ndims], 
        index_src1[m_ndims], 
        index_src2[m_ndims],
        index_src3[m_ndims];
    _ArrayCopy1D_(index_transverse, index_dst  , m_ndims);
    _ArrayCopy1D_(index_transverse, index_src0, m_ndims);
    _ArrayCopy1D_(index_transverse, index_src1 , m_ndims);
    _ArrayCopy1D_(index_transverse, index_src2 , m_ndims);
    _ArrayCopy1D_(index_transverse, index_src3, m_ndims);

    for (int i_dst = 0; i_dst < n_dst; i_dst++) {

      double xi_dst = ((double) i_dst + 0.5) / ((double) stride) - 0.5;

      int i_src_1  = std::floor(xi_dst);
      int i_src_2  = std::ceil(xi_dst);
      int i_src_0 = i_src_1 - 1;
      int i_src_3 = i_src_2 + 1;

      double alpha = ((double)i_src_2 - xi_dst) / ((double)i_src_2 - (double)i_src_1);

      index_dst[a_dir] = i_dst + a_ngpt;
      int p; _ArrayIndex1D_(m_ndims, dim_dst_wg, index_dst, 0, p);

      index_src0[a_dir] = i_src_0 + a_ngpt;
      int q0; _ArrayIndex1D_(m_ndims, dim_src_wg, index_src0, 0, q0);

      index_src1[a_dir] = i_src_1 + a_ngpt;
      int q1; _ArrayIndex1D_(m_ndims, dim_src_wg, index_src1, 0, q1);

      index_src2[a_dir] = i_src_2 + a_ngpt;
      int q2; _ArrayIndex1D_(m_ndims, dim_src_wg, index_src2, 0, q2);

      index_src3[a_dir] = i_src_3 + a_ngpt;
      int q3; _ArrayIndex1D_(m_ndims, dim_src_wg, index_src3, 0, q3);

      double c0 = -((-2.0 + alpha)*(-1.0 + alpha)*alpha)/6.0;
      double c1 = ((-2.0 + alpha)*(-1.0 + alpha)*(1.0 + alpha))/2.0;
      double c2 = (alpha*(2.0 + alpha - alpha*alpha))/2.0;
      double c3 = (alpha*(-1.0 + alpha*alpha))/6.0;

      for (int v = 0; v < a_nvars; v++) {

        a_u_dst[p*a_nvars+v] =    c0 * a_u_src[q0*a_nvars+v]
                                + c1 * a_u_src[q1*a_nvars+v]
                                + c2 * a_u_src[q2*a_nvars+v]
                                + c3 * a_u_src[q3*a_nvars+v];

      }

    }

  }

  return;
}


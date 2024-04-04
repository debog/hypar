/*! @file InterpolateGlobalnDVar.c
    @brief Functions to interpolate from one grid to another
    @author Debojyoti Ghosh
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <arrayfunctions.h>
#include <mathfunctions.h>

/*! Is the input number an integer power of 2? */
static int isPowerOfTwo(int x)
{
  if (x == 0)  return 0;

  while (x > 1) {
    if (x%2 != 0) return 0;
    x /= 2;
  }
  return 1;
}

/*! Coarsen along a given dimension - the source and destination must have the same sizes
 *  along all other dimensions.
 *
 *  Note that the arrays must *not* have any ghost points!
 *
 *  Currently this function can only handle coarsening factors that are integer
 *  powers of 2.
*/
static int coarsen1D( const int* const    a_dim_src, /*!< Grid size of source data */
                      const int* const    a_dim_dst, /*!< Grid size of destination data */
                      const double* const a_u_src,   /*!< Source solution */
                      double* const       a_u_dst,   /*!< Destination solution */
                      const int           a_dir,     /*!< Dimension along which to coarsen */
                      const int           a_nvars,   /*!< Number of vector components of the solution */
                      const int           a_ngpt,    /*!< Number of ghost points */
                      const int           a_ndims    /*!< Number of spatial dimensions */
                   )
{
  for (int d = 0; d < a_ndims; d++) {
    if ((d != a_dir) && (a_dim_src[d] != a_dim_dst[d])) {
      fprintf(stderr, "Error in coarsen1D() -\n");
      fprintf(stderr, " a_dim_src[%d] != a_dim_dst[%d]\n", d, d);
      return 1;
    }
  }

  int n_src = a_dim_src[a_dir];
  int n_dst = a_dim_dst[a_dir];
  if (n_dst > n_src) {
    fprintf(stderr, "Error in coarsen1D() -\n");
    fprintf(stderr, " destination grid is finer than source grid along a_dir!\n");
    return 1;
  }

  double fac = ((double) n_src) / ((double) n_dst);
  int stride = (int) fac;
  if (absolute(((double)stride)-fac) > _MACHINE_ZERO_) {
    fprintf(stderr, "Error in coarsen1D() -\n");
    fprintf(stderr, "  non-integer coarsening factor!\n");
    return 1;
  }

  /* set interpolation coefficients depending on desired order */
  double c0, c1, c2, c3, c4, c5;
//  if (a_interp_order == 2) {
//    c0 = c5 = 0.0;
//    c1 = c4 = 0.0;
//    c2 = c3 = 0.5;
//  } else if (a_interp_order == 4) {
//    c0 = c5 = 0.0;
//    c1 = c4 = -1.0/16.0;
//    c2 = c3 = 9.0/16.0;
//  } else if (a_interp_order == 6) {
    c0 = c5 = 3.0/256.0;
    c1 = c4 = -25.0/256.0;
    c2 = c3 = 150.0/256.0;
//  } else {
//    fprintf(stderr,"Invalid value of interpolation order!\n");
//    return 1;
//  }

  /* create bounds for the transverse loop, i.e., to loop over
   * all 1D lines along dimension "a_dir" */
  int bounds_transverse[a_ndims];
  _ArrayCopy1D_(a_dim_src, bounds_transverse, a_ndims);
  bounds_transverse[a_dir] =  1;

  int index_transverse[a_ndims], done = 0;
  _ArraySetValue_(index_transverse, a_ndims, 0);
  while (!done) {

    int index_dst[a_ndims], index_src[a_ndims];
    _ArrayCopy1D_(index_transverse, index_dst, a_ndims);
    _ArrayCopy1D_(index_transverse, index_src, a_ndims);

    for (int i_dst = 0; i_dst < n_dst; i_dst++) {

      int i_m1 = i_dst*stride + (stride/2-1);
      int i_m3 = i_m1 - 2;
      int i_m2 = i_m1 - 1;
      int i_p1 = i_m1 + 1;
      int i_p2 = i_m1 + 2;
      int i_p3 = i_m1 + 3;

      int p;
      index_dst[a_dir] = i_dst;
      _ArrayIndex1D_(a_ndims, a_dim_dst, index_dst, a_ngpt, p);

      int p_m3;
      index_src[a_dir] = i_m3;
      _ArrayIndex1D_(a_ndims, a_dim_src, index_src, a_ngpt, p_m3);

      int p_m2;
      index_src[a_dir] = i_m2;
      _ArrayIndex1D_(a_ndims, a_dim_src, index_src, a_ngpt, p_m2);

      int p_m1;
      index_src[a_dir] = i_m1;
      _ArrayIndex1D_(a_ndims, a_dim_src, index_src, a_ngpt, p_m1);

      int p_p1;
      index_src[a_dir] = i_p1;
      _ArrayIndex1D_(a_ndims, a_dim_src, index_src, a_ngpt, p_p1);

      int p_p2;
      index_src[a_dir] = i_p2;
      _ArrayIndex1D_(a_ndims, a_dim_src, index_src, a_ngpt, p_p2);

      int p_p3;
      index_src[a_dir] = i_p3;
      _ArrayIndex1D_(a_ndims, a_dim_src, index_src, a_ngpt, p_p3);

      for (int v = 0; v < a_nvars; v++) {
        double val =    c0 * a_u_src[p_m3*a_nvars+v]
                      + c1 * a_u_src[p_m2*a_nvars+v]
                      + c2 * a_u_src[p_m1*a_nvars+v]
                      + c3 * a_u_src[p_p1*a_nvars+v]
                      + c4 * a_u_src[p_p2*a_nvars+v]
                      + c5 * a_u_src[p_p3*a_nvars+v];
        a_u_dst[p*a_nvars+v] = val;
      }

    }

    _ArrayIncrementIndex_(a_ndims, bounds_transverse, index_transverse, done);

  }

  return 0;
}

/*! Refine along a given dimension - the source and destination must have the same sizes
 *  along all other dimensions.
 *
 *  Note that the arrays must *not* have any ghost points!
 *
 *  Currently this function can only handle refinement factors that are integer
 *  powers of 2.
*/
static int refine1D(const int* const     a_dim_src, /*!< Grid size of source data */
                    const int* const     a_dim_dst, /*!< Grid size of destination data */
                    const double* const  a_u_src,   /*!< Source solution */
                    double* const        a_u_dst,   /*!< Destination solution */
                    const int            a_dir,     /*!< Dimension along which to coarsen */
                    const int            a_nvars,   /*!< Number of vector components of the solution */
                    const int            a_ngpt,    /*!< Number of ghost points */
                    const int            a_ndims    /*!< Number of spatial dimensions */
                  )
{
  for (int d = 0; d < a_ndims; d++) {
    if ((d != a_dir) && (a_dim_src[d] != a_dim_dst[d])) {
      fprintf(stderr, "Error in refine1D() -\n");
      fprintf(stderr, " a_dim_src[%d] != a_dim_dst[%d]\n", d, d);
      return 1;
    }
  }

  int n_src = a_dim_src[a_dir];
  int n_dst = a_dim_dst[a_dir];
  if (n_dst < n_src) {
    fprintf(stderr, "Error in refine1D() -\n");
    fprintf(stderr, "  destination grid is coarser than source grid along a_dir!\n");
    return 1;
  }

  double fac = ((double) n_dst) / ((double) n_src);
  int stride = (int) fac;
  if (absolute(((double)stride)-fac) > _MACHINE_ZERO_) {
    fprintf(stderr, "Error in refine1D() -\n");
    fprintf(stderr, "  non-integer refinement factor!\n");
    return 1;
  }

  /* create bounds for the transverse loop, i.e., to loop over
   * all 1D lines along dimension "a_dir" */
  int bounds_transverse[a_ndims];
  _ArrayCopy1D_(a_dim_src, bounds_transverse, a_ndims);
  bounds_transverse[a_dir] =  1;

  int index_transverse[a_ndims], done = 0;
  _ArraySetValue_(index_transverse, a_ndims, 0);
  while (!done) {

    int index_dst [a_ndims],
        index_src0[a_ndims],
        index_src1[a_ndims],
        index_src2[a_ndims],
        index_src3[a_ndims],
        index_src4[a_ndims],
        index_src5[a_ndims];
    _ArrayCopy1D_(index_transverse, index_dst , a_ndims);
    _ArrayCopy1D_(index_transverse, index_src0, a_ndims);
    _ArrayCopy1D_(index_transverse, index_src1, a_ndims);
    _ArrayCopy1D_(index_transverse, index_src2, a_ndims);
    _ArrayCopy1D_(index_transverse, index_src3, a_ndims);
    _ArrayCopy1D_(index_transverse, index_src4, a_ndims);
    _ArrayCopy1D_(index_transverse, index_src5, a_ndims);

    for (int i_dst = 0; i_dst < n_dst; i_dst++) {

      double xi_dst = ((double) i_dst + 0.5) / ((double) stride) - 0.5;

      int i_src_2  = floor(xi_dst);
      int i_src_3  = ceil(xi_dst);
      int i_src_0 = i_src_2 - 2;
      int i_src_1 = i_src_2 - 1;
      int i_src_4 = i_src_3 + 1;
      int i_src_5 = i_src_3 + 2;

      double alpha = (xi_dst - (double)i_src_2) / ((double)i_src_3 - (double)i_src_2);

      index_dst[a_dir] = i_dst;
      int p; _ArrayIndex1D_(a_ndims, a_dim_dst, index_dst, a_ngpt, p);

      index_src0[a_dir] = i_src_0;
      int q0; _ArrayIndex1D_(a_ndims, a_dim_src, index_src0, a_ngpt, q0);

      index_src1[a_dir] = i_src_1;
      int q1; _ArrayIndex1D_(a_ndims, a_dim_src, index_src1, a_ngpt, q1);

      index_src2[a_dir] = i_src_2;
      int q2; _ArrayIndex1D_(a_ndims, a_dim_src, index_src2, a_ngpt, q2);

      index_src3[a_dir] = i_src_3;
      int q3; _ArrayIndex1D_(a_ndims, a_dim_src, index_src3, a_ngpt, q3);

      index_src4[a_dir] = i_src_4;
      int q4; _ArrayIndex1D_(a_ndims, a_dim_src, index_src4, a_ngpt, q4);

      index_src5[a_dir] = i_src_5;
      int q5; _ArrayIndex1D_(a_ndims, a_dim_src, index_src5, a_ngpt, q5);

      /* set interpolation coefficients depending on desired order */
      double c0, c1, c2, c3, c4, c5;
//      if (a_interp_order == 2) {
//        c0 = 0.0;
//        c1 = 0.0;
//        c2 = (1.0-alpha);
//        c3 = alpha;
//        c4 = 0.0;
//        c5 = 0.0;
//      } else if (a_interp_order == 4) {
//        c0 = 0.0;
//        c1 = -((-2.0 + alpha)*(-1.0 + alpha)*alpha)/6.0;
//        c2 = ((-2.0 + alpha)*(-1.0 + alpha)*(1.0 + alpha))/2.0;
//        c3 = (alpha*(2.0 + alpha - alpha*alpha))/2.0;
//        c4 = (alpha*(-1.0 + alpha*alpha))/6.0;
//        c5 = 0.0;
//      } else if (a_interp_order == 6) {
        c0 = -((-3.0 + alpha)*(-2.0 + alpha)*(-1.0 + alpha)*alpha*(1.0 + alpha))/120.0;
        c1 = ((-3.0 + alpha)*(-2.0 + alpha)*(-1.0 + alpha)*alpha*(2.0 + alpha))/24.0;
        c2 = -((-3.0 + alpha)*(-2.0 + alpha)*(-1.0 + alpha)*(1.0 + alpha)*(2.0 + alpha))/12.0;
        c3 = ((-3.0 + alpha)*(-2.0 + alpha)*alpha*(1.0 + alpha)*(2.0 + alpha))/12.0;
        c4 = -((-3.0 + alpha)*(-1.0 + alpha)*alpha*(1.0 + alpha)*(2.0 + alpha))/24.0;
        c5 = (alpha*(4.0 - 5.0*alpha*alpha + alpha*alpha*alpha*alpha))/120.0;
//      } else {
//        fprintf(stderr,"Invalid value of interpolation order!\n");
//        return 1;
//      }

      for (int v = 0; v < a_nvars; v++) {

        a_u_dst[p*a_nvars+v] =    c0 * a_u_src[q0*a_nvars+v]
                                + c1 * a_u_src[q1*a_nvars+v]
                                + c2 * a_u_src[q2*a_nvars+v]
                                + c3 * a_u_src[q3*a_nvars+v]
                                + c4 * a_u_src[q4*a_nvars+v]
                                + c5 * a_u_src[q5*a_nvars+v];

      }

    }

    _ArrayIncrementIndex_(a_ndims, bounds_transverse, index_transverse, done);

  }

  return 0;
}

/*! Interpolate n-dimensional data from one grid to another of a desired resolution.
    Note that along each dimension, the ratio of the number of grid points
    in the source grid and that in the destination grid must be an integer
    power of 2 (negative or positive).

    The source data *must* be global. It will get deallocated at the end
    of this function. It *must* have the appropriate number of ghost points.

    The incoming pointer for destination data must be NULL. After this function
    is executed, it will point to a chunk of memory with the interpolated solution.
    This is the *global* solution. It will have the specified number of ghost
    points appropriately filled.
*/
int InterpolateGlobalnDVar( const int* const  a_dim_dst, /*!< grid dimensions to interpolate to */
                            double** const    a_u_dst,   /*!< pointer to array containing interpolated data */
                            const int* const  a_dim_src, /*!< grid dimensions to interpolate from */
                            double* const     a_u_src,   /*!< pointer to array containing data to interpolate from */
                            const int         a_nvars,   /*!< Number of vector components of the solution */
                            const int         a_ghosts,  /*!< Number of ghost points */
                            const int         a_ndims,   /*!< Number of spatial dimensions */
                            const int* const  a_periodic /*!< Is the domain periodic or not along each dimension */
                         )
{
  if ((*a_u_dst) != NULL) {
    fprintf(stderr, "Error in InterpolateGlobalnDVar() - \n");
    fprintf(stderr, "  a_u_dst is not NULL!\n");
    return 1;
  }

  if (a_u_src == NULL) {
    fprintf(stderr, "Error in InterpolateGlobalnDVar() - \n");
    fprintf(stderr, "  a_u_src is NULL!\n");
    return 1;
  }

  int dim_to[a_ndims], dim_from[a_ndims];

  double *u_from;
  double *u_to;

  _ArrayCopy1D_(a_dim_src, dim_to, a_ndims);
  u_to = a_u_src;
  u_from = NULL;

  for (int dir = 0; dir < a_ndims; dir++) {

    _ArrayCopy1D_(dim_to, dim_from, a_ndims);
    dim_to[dir] = a_dim_dst[dir];

    if (dim_from[dir] == dim_to[dir]) continue;

    double fac = (dim_to[dir] > dim_from[dir] ?
                      (double)dim_to[dir]/(double)dim_from[dir]
                    : (double)dim_from[dir]/(double)dim_to[dir] );

    if (!isPowerOfTwo((int)fac)) {
      fprintf(stderr,"Error in interpolate() - \n");
      fprintf(stderr,"  refinement/coarsening factor not a power of 2!\n");
      return 1;
    }

    if (u_from != NULL) free(u_from);

    u_from = u_to;

    {
      long size = (long) a_nvars;
      for (int d = 0; d < a_ndims; d++) {
        size *= (long) (dim_to[d] + 2*a_ghosts);
      }
      u_to = (double*) calloc (size, sizeof(double));
    }

    fillGhostCells( dim_from,
                    a_ghosts,
                    u_from,
                    a_nvars,
                    a_ndims,
                    a_periodic );

    if (dim_to[dir] < dim_from[dir]) {
      int retval = coarsen1D( dim_from,
                              dim_to,
                              u_from,
                              u_to,
                              dir,
                              a_nvars,
                              a_ghosts,
                              a_ndims );
      if (retval) return retval;
    } else {
      int retval = refine1D(  dim_from,
                              dim_to,
                              u_from,
                              u_to,
                              dir,
                              a_nvars,
                              a_ghosts,
                              a_ndims );
      if (retval) return retval;
    }

  }

  /* dim_to should be equal to a_dim_dst now */
  for (int d = 0; d < a_ndims; d++) {
    if (dim_to[d] != a_dim_dst[d]) {
      fprintf(stderr,"Error in InterpolateGlobalnDVar() - \n");
      fprintf(stderr,"  dim_to[%d] (%d) != a_dim_dst[%d] (%d)!\n",
              d, dim_to[d], d, a_dim_dst[d]);
      return 1;
    }
  }

  if (u_from != NULL) free(u_from);
  (*a_u_dst) = u_to;

  return 0;
}

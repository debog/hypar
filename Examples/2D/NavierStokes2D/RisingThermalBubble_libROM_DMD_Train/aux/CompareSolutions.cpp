/*
 * This code computes the difference between two sets of binary
 * solutions that are specified by their filename prefix, i.e.
 * "op_", as the L1, L2, and Linf norm of this difference.
 *
 * The solutions must all be in the binary format. The difference is
 * written out in the binary format as well.
 *
 * If the size of the two solutions differ, they must differ by an
 * integer power of 2 along each dimension.
 * 6th-order central interpolation is used to interpolate the
 * reference solution onto the grid of the solution.
 *
 * Please compile it using the C99 standard.
 *   gcc -std=c99 CompareSolutions.c -lm [-o <executable.name>]
*/

#include <string>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define _MACHINE_ZERO_ 1.0e-14
#define _MAX_STRING_SIZE_ 500

#define _ArraySetValue_(x,size,value)\
  {\
    int arraycounter;\
    for (arraycounter = 0; arraycounter < (size); arraycounter++)  x[arraycounter] = (value);\
  }

#define _ArrayCopy1D_(x,y,size) \
  { \
    int arraycounter; \
    for (arraycounter = 0; arraycounter < size; arraycounter++) y[arraycounter] = x[arraycounter];\
  }

#define _ArrayIndex1D_(N,imax,i,ghost,index)  \
  { \
    index = i[N-1]+(ghost); \
    int arraycounter; \
    for (arraycounter = (N)-2; arraycounter > -1; arraycounter--) { \
      index = ((index*(imax[arraycounter]+2*(ghost))) + (i[arraycounter]+(ghost))); \
    } \
  }

#define _ArrayIncrementIndex_(N,imax,i,done) \
  { \
    int arraycounter = 0; \
    while (arraycounter < (N)) { \
      if (i[arraycounter] == imax[arraycounter]-1) { \
        i[arraycounter] = 0; \
        arraycounter++; \
      } else { \
        i[arraycounter]++; \
        break; \
      } \
    } \
    if (arraycounter == (N)) done = 1; \
    else          done = 0; \
  }

#define _ArrayAXPY_(x,a,y,size)\
  { \
    int arraycounter; \
    for (arraycounter=0; arraycounter<size; arraycounter++) y[arraycounter] += a*x[arraycounter];\
  }

int ArrayCopynD(int                 ndims,
                const double* const x,
                double* const       y,
                const int* const    dim,
                int                 g1,
                int                 g2,
                int                 nvars )
{
  if (!y) {
    fprintf(stderr,"Error in ArrayCopynD(): array \"y\" not allocated.\n");
    return(1);
  }
  if (!x) {
    fprintf(stderr,"Error in ArrayCopynD(): array \"x\" not allocated.\n");
    return(1);
  }
  int done = 0;
  int index[ndims];
  _ArraySetValue_(index,ndims,0);
  while (!done) {
    int p1, p2;
    _ArrayIndex1D_(ndims,dim,index,g1,p1);
    _ArrayIndex1D_(ndims,dim,index,g2,p2);
    _ArrayCopy1D_((x+p1*nvars),(y+p2*nvars),nvars);
    _ArrayIncrementIndex_(ndims,dim,index,done);
  }
  return(0);
}

double ArrayMaxnD(int    nvars,
                  int    ndims,
                  int    *dim,
                  int    ghosts,
                  int    *index,
                  double *x )
{
  double sum = 0;
  int done = 0; _ArraySetValue_(index,ndims,0);
  while (!done) {
    int p; _ArrayIndex1D_(ndims,dim,index,ghosts,p);
    int v;
    for (v=0; v<nvars; v++) {
      double term = ( x[p*nvars+v]>0 ? x[p*nvars+v] : -x[p*nvars+v] );
      if (term > sum) sum = term;
    }
    _ArrayIncrementIndex_(ndims,dim,index,done);
  }
  return(sum);
}

double ArraySumAbsnD(int     nvars,
                     int     ndims,
                     int     *dim,
                     int     ghosts,
                     int     *index,
                     double  *x )
{
  double sum = 0;
  int done = 0; _ArraySetValue_(index,ndims,0);
  while (!done) {
    int p; _ArrayIndex1D_(ndims,dim,index,ghosts,p);
    int v; for (v=0; v<nvars; v++) sum += ( x[p*nvars+v]>0 ? x[p*nvars+v] : -x[p*nvars+v] );
    _ArrayIncrementIndex_(ndims,dim,index,done);
  }
  return(sum);
}

double ArraySumSquarenD(int    nvars,
                        int    ndims,
                        int    *dim,
                        int    ghosts,
                        int    *index,
                        double *x )
{
  double sum = 0;
  int done = 0; _ArraySetValue_(index,ndims,0);
  while (!done) {
    int p; _ArrayIndex1D_(ndims,dim,index,ghosts,p);
    int v; for (v=0; v<nvars; v++) sum += (x[p*nvars+v]*x[p*nvars+v]);
    _ArrayIncrementIndex_(ndims,dim,index,done);
  }
  return(sum);
}

void fillGhostCells( const int* const a_dim,
                     const int        a_ngpt,
                     double* const    a_u,
                     const int        a_nvars,
                     const int        a_ndims,
                     const int* const a_periodic )
{
  for (int d = 0; d < a_ndims; d++) {

    int bounds[a_ndims];
    _ArrayCopy1D_(a_dim, bounds, a_ndims);
    bounds[d] = a_ngpt;

    int index[a_ndims];
    _ArraySetValue_(index, a_ndims, 0);

    if (a_periodic[d]) {

      /* periodic case */

      int done = 0;
      while (!done) {

        {
          /* low end - face = 1 */

          int p_gpt = 0,
              p_int = 0;

          int index_gpt[a_ndims];
          _ArrayCopy1D_(index, index_gpt, a_ndims);
          index_gpt[d] -= a_ngpt;
          _ArrayIndex1D_(a_ndims, a_dim, index_gpt, a_ngpt, p_gpt);

          int index_int[a_ndims];
          _ArrayCopy1D_(index, index_int, a_ndims);
          index_int[d] += (a_dim[d]-a_ngpt);
          _ArrayIndex1D_(a_ndims, a_dim, index_int, a_ngpt, p_int);

          _ArrayCopy1D_((a_u+a_nvars*p_int), (a_u+a_nvars*p_gpt), a_nvars);
        }

        {
          /* high end - face = -1 */

          int p_gpt = 0,
              p_int = 0;

          int index_gpt[a_ndims];
          _ArrayCopy1D_(index, index_gpt, a_ndims);
          index_gpt[d] += a_dim[d];
          _ArrayIndex1D_(a_ndims, a_dim, index_gpt, a_ngpt, p_gpt);

          int index_int[a_ndims];
          _ArrayCopy1D_(index, index_int, a_ndims);
          _ArrayIndex1D_(a_ndims, a_dim, index_int, a_ngpt, p_int);

          _ArrayCopy1D_((a_u+a_nvars*p_int), (a_u+a_nvars*p_gpt), a_nvars);
        }

        _ArrayIncrementIndex_(a_ndims, bounds, index, done);

      }

    } else {

      /* not periodic - extrapolate */

      int done = 0;
      while (!done) {

        {
          /* low end - face = 1 */

          int p_gpt = 0,
              p_int_0 = 0,
              p_int_1 = 0,
              p_int_2 = 0,
              p_int_3 = 0;

          int index_gpt[a_ndims];
          _ArrayCopy1D_(index, index_gpt, a_ndims);
          index_gpt[d] -= a_ngpt;
          _ArrayIndex1D_(a_ndims, a_dim, index_gpt, a_ngpt, p_gpt);

          int index_int[a_ndims];
          _ArrayCopy1D_(index, index_int, a_ndims);

          index_int[d] = 0;
          _ArrayIndex1D_(a_ndims, a_dim, index_int, a_ngpt, p_int_0);
          index_int[d]++;
          _ArrayIndex1D_(a_ndims, a_dim, index_int, a_ngpt, p_int_1);
          index_int[d]++;
          _ArrayIndex1D_(a_ndims, a_dim, index_int, a_ngpt, p_int_2);
          index_int[d]++;
          _ArrayIndex1D_(a_ndims, a_dim, index_int, a_ngpt, p_int_3);

          double alpha = - (double) (a_ngpt - index[d]);
          double c0 = -((-2.0 + alpha)*(-1.0 + alpha)*alpha)/6.0;
          double c1 = ((-2.0 + alpha)*(-1.0 + alpha)*(1.0 + alpha))/2.0;
          double c2 = (alpha*(2.0 + alpha - alpha*alpha))/2.0;
          double c3 = (alpha*(-1.0 + alpha*alpha))/6.0;

          for (int v = 0; v < a_nvars; v++) {

            a_u[p_gpt*a_nvars+v] =    c0 * a_u[p_int_0*a_nvars+v]
                                    + c1 * a_u[p_int_1*a_nvars+v]
                                    + c2 * a_u[p_int_2*a_nvars+v]
                                    + c3 * a_u[p_int_3*a_nvars+v];

          }

        }

        {
          /* high end - face = -1 */

          int p_gpt = 0,
              p_int_0 = 0,
              p_int_1 = 0,
              p_int_2 = 0,
              p_int_3 = 0;

          int index_gpt[a_ndims];
          _ArrayCopy1D_(index, index_gpt, a_ndims);
          index_gpt[d] += a_dim[d];
          _ArrayIndex1D_(a_ndims, a_dim, index_gpt, a_ngpt, p_gpt);

          int index_int[a_ndims];
          _ArrayCopy1D_(index, index_int, a_ndims);

          index_int[d] = a_dim[d]-1;
          _ArrayIndex1D_(a_ndims, a_dim, index, a_ngpt, p_int_0);
          index_int[d]--;
          _ArrayIndex1D_(a_ndims, a_dim, index, a_ngpt, p_int_1);
          index_int[d]--;
          _ArrayIndex1D_(a_ndims, a_dim, index, a_ngpt, p_int_2);
          index_int[d]--;
          _ArrayIndex1D_(a_ndims, a_dim, index, a_ngpt, p_int_3);

          double alpha = - (double) (index[d]+1);
          double c0 = -((-2.0 + alpha)*(-1.0 + alpha)*alpha)/6.0;
          double c1 = ((-2.0 + alpha)*(-1.0 + alpha)*(1.0 + alpha))/2.0;
          double c2 = (alpha*(2.0 + alpha - alpha*alpha))/2.0;
          double c3 = (alpha*(-1.0 + alpha*alpha))/6.0;

          for (int v = 0; v < a_nvars; v++) {

            a_u[p_gpt*a_nvars+v] =    c0 * a_u[p_int_0*a_nvars+v]
                                    + c1 * a_u[p_int_1*a_nvars+v]
                                    + c2 * a_u[p_int_2*a_nvars+v]
                                    + c3 * a_u[p_int_3*a_nvars+v];

          }

        }

        _ArrayIncrementIndex_(a_ndims, bounds, index, done);

      }

    }

  }

  return;
}

int isPowerOfTwo(int x)
{
  if (x == 0)  return 0;

  while (x > 1) {
    if (x%2 != 0) return 0;
    x /= 2;
  }
  return 1;
}

int coarsen1D( const int* const    a_dim_src,
               const int* const    a_dim_dst,
               const double* const a_u_src,
               double* const       a_u_dst,
               const int           a_dir,
               const int           a_nvars,
               const int           a_ngpt,
               const int           a_ndims )
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
  if (abs(((double)stride)-fac) > _MACHINE_ZERO_) {
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

int refine1D(const int* const     a_dim_src,
             const int* const     a_dim_dst,
             const double* const  a_u_src,
             double* const        a_u_dst,
             const int            a_dir,
             const int            a_nvars,
             const int            a_ngpt,
             const int            a_ndims )
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
  if (abs(((double)stride)-fac) > _MACHINE_ZERO_) {
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

int InterpolateGlobalnDVar( const int* const    a_dim_dst,
                            double*  const      a_u_dst,
                            const int* const    a_dim_src,
                            const double* const a_u_src,
                            const int           a_nvars,
                            const int           a_ndims,
                            const int* const    a_periodic )
{
  static const int ghosts = 3;

  long size_src = a_nvars;
  for (int n = 0; n < a_ndims; n++) {
    size_src *= (a_dim_src[n] + 2*ghosts);
  }
  double* u_src = (double*) calloc (size_src, sizeof(double));
  ArrayCopynD(  a_ndims,
                a_u_src,
                u_src,
                a_dim_src,
                0,
                ghosts,
                a_nvars );

  int dim_to[a_ndims], dim_from[a_ndims];

  double *u_from;
  double *u_to;

  _ArrayCopy1D_(a_dim_src, dim_to, a_ndims);
  u_to = u_src;
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
        size *= (long) (dim_to[d] + 2*ghosts);
      }
      u_to = (double*) calloc (size, sizeof(double));
    }

    fillGhostCells( dim_from,
                    ghosts,
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
                              ghosts,
                              a_ndims );
      if (retval) return retval;
    } else {
      int retval = refine1D(  dim_from,
                              dim_to,
                              u_from,
                              u_to,
                              dir,
                              a_nvars,
                              ghosts,
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
  ArrayCopynD(  a_ndims,
                u_to,
                a_u_dst,
                a_dim_dst,
                ghosts,
                0,
                a_nvars );

  free(u_to);

  return 0;
}

typedef struct solution_object {
  int     ndims;
  int     nvars;
  int*    size;
  double* x;
  double* u;
  int*    is_periodic;
} SolObj;

void deleteSolObj(SolObj* const a_solution)
{
  if (a_solution->size) free(a_solution->size);
  if (a_solution->x) free(a_solution->x);
  if (a_solution->u) free(a_solution->u);
  if (a_solution->is_periodic) free(a_solution->is_periodic);
}

void createSolObj(const int         a_ndims,
                  const int         a_nvars,
                  const int* const  a_size,
                  SolObj* const     a_solution )
{
  a_solution->ndims = a_ndims;
  a_solution->nvars = a_nvars;

  a_solution->size = (int*) calloc(a_ndims, sizeof(int));
  _ArrayCopy1D_(a_size, a_solution->size, a_ndims);

  long sizex = 0;
  long sizeu = a_solution->nvars;
  for (int n=0; n<a_solution->ndims; n++) {
    sizex += a_solution->size[n];
    sizeu *= a_solution->size[n];
  }
  a_solution->x = (double*) calloc (sizex, sizeof(double));
  a_solution->u = (double*) calloc (sizeu, sizeof(double));

  a_solution->is_periodic = (int*) calloc (a_ndims, sizeof(int));
  _ArraySetValue_(a_solution->is_periodic, a_ndims, 1);
  return;
}

int readBinaryFile( const char* const a_filename,
                    SolObj* const     a_solution )
{
  FILE* in;
  in = fopen(a_filename,"rb");

  if (!in) {
    return 1;
  }

  printf("Reading binary file %s:\n", a_filename);

  fread(&a_solution->ndims, sizeof(int), 1, in);
  fread(&a_solution->nvars, sizeof(int), 1, in);
  printf("  ndims: %d\n", a_solution->ndims);
  printf("  nvars: %d\n", a_solution->nvars);

  a_solution->size = (int*) calloc(a_solution->ndims, sizeof(int));
  fread(a_solution->size, sizeof(int), a_solution->ndims, in);
  printf("  size: ");
  for (int n=0; n<a_solution->ndims; n++) {
    printf("%d ", a_solution->size[n]);
  }
  printf("\n");

  long sizex = 0;
  long sizeu = a_solution->nvars;
  for (int n=0; n<a_solution->ndims; n++) {
    sizex += a_solution->size[n];
    sizeu *= a_solution->size[n];
  }

  a_solution->x = (double*) calloc (sizex, sizeof(double));
  a_solution->u = (double*) calloc (sizeu, sizeof(double));
  fread(a_solution->x, sizeof(double), sizex, in);
  fread(a_solution->u, sizeof(double), sizeu, in);

  a_solution->is_periodic = (int*) calloc(a_solution->ndims,sizeof(int));
  _ArraySetValue_(a_solution->is_periodic, a_solution->ndims, 1);

  printf("  Done.\n");
  fclose(in);
  return 0;
}

int checkConsistency( const SolObj* const a_s1,
                      const SolObj* const a_s2 )
{
  if (a_s1->ndims != a_s2->ndims) {
    fprintf(stderr,"Error in checkConsistency(): ndims is not the same.\n");
    return 1;
  }
  if (a_s1->nvars != a_s2->nvars) {
    fprintf(stderr,"Error in checkConsistency(): nvars is not the same.\n");
    return 1;
  }
  return 0;
}

int ComputeDiffNorm(  SolObj* const a_s1,
                      SolObj* const a_s2,
                      double* const a_diff_norm )
{
  static const double tolerance = 1e-15;
  int ierr;
  int index[a_s2->ndims];

  ierr = checkConsistency(a_s1, a_s2);
  if (ierr) return ierr;

  long npoints = 1;
  long size = a_s2->nvars;
  for (int n=0; n<a_s2->ndims; n++) {
    npoints *= (a_s2->size[n]);
    size *= a_s2->size[n];
  }

  /* calculate solution norms (for relative error) */
  double solution_norm[3] = {0.0,0.0,0.0};

  /* L1 */
  {
    double sum = ArraySumAbsnD ( a_s2->nvars,
                          a_s2->ndims,
                          a_s2->size,
                          0,
                          index,
                          a_s2->u );
    solution_norm[0] = sum / ((double)npoints);
  }
  /* L2 */
  {
    double sum = ArraySumSquarenD( a_s2->nvars,
                            a_s2->ndims,
                            a_s2->size,
                            0,
                            index,
                            a_s2->u );
    solution_norm[1] = sqrt(sum/((double)npoints));
  }
  /* Linf */
  {
    double sum = ArrayMaxnD( a_s2->nvars,
                      a_s2->ndims,
                      a_s2->size,
                      0,
                      index,
                      a_s2->u );
    solution_norm[2] = sum;
  }

  _ArrayAXPY_(a_s1->u,-1.0,a_s2->u,size);

  /* L1 */
  {
    double sum = ArraySumAbsnD ( a_s2->nvars,
                          a_s2->ndims,
                          a_s2->size,
                          0,
                          index,
                          a_s2->u );
    a_diff_norm[0] = sum / ((double)npoints);
  }
  /* L2 */
  {
    double sum = ArraySumSquarenD( a_s2->nvars,
                            a_s2->ndims,
                            a_s2->size,
                            0,
                            index,
                            a_s2->u );
    a_diff_norm[1] = sqrt(sum/((double)npoints));
  }
  /* Linf */
  {
    double sum = ArrayMaxnD( a_s2->nvars,
                      a_s2->ndims,
                      a_s2->size,
                      0,
                      index,
                      a_s2->u );
    a_diff_norm[2] = sum;
  }

  if (    (solution_norm[0] > tolerance)
      &&  (solution_norm[1] > tolerance)
      &&  (solution_norm[2] > tolerance) ) {
    a_diff_norm[0] /= solution_norm[0];
    a_diff_norm[1] /= solution_norm[1];
    a_diff_norm[2] /= solution_norm[2];
  }

  return 0;
}

int writeBinary(  const int           ndims,
                  const int           nvars,
                  const int* const    dim,
                  const double* const x,
                  const double* const u,
                  const char* const   f )
{
  int size, d;
  size_t bytes;
  FILE *out;
  out = fopen(f,"wb");
  if (!out) {
    fprintf(stderr,"Error: could not open %s for writing.\n",f);
    return(1);
  }

  /* write ndims, nvars */
  bytes = fwrite(&ndims,sizeof(int),1,out);
  if ((int)bytes != 1) {
    fprintf(stderr,"Error in WriteBinary(): Unable to write ndims to output file.\n");
  }
  bytes = fwrite(&nvars,sizeof(int),1,out);
  if ((int)bytes != 1) {
    fprintf(stderr,"Error in WriteBinary(): Unable to write nvars to output file.\n");
  }

  /* write dimensions */
  bytes = fwrite(dim,sizeof(int),ndims,out);
  if ((int)bytes != ndims) {
    fprintf(stderr,"Error in WriteBinary(): Unable to write dimensions to output file.\n");
  }

  /* write grid */
  size = 0;
  for (d = 0; d < ndims; d++) size += dim[d];
  bytes = fwrite(x,sizeof(double),size,out);
  if ((int)bytes != size) {
    fprintf(stderr,"Error in WriteBinary(): Unable to write grid to output file.\n");
  }

  /* write solution */
  size = 1;
  for (d = 0; d < ndims; d++) size *= dim[d]; size *= nvars;
  bytes = fwrite(u,sizeof(double),size,out);
  if ((int)bytes != size) {
    fprintf(stderr,"Error in WriteBinary(): Unable to write solution to output file.\n");
  }

  fclose(out);
  return(0);
}

int ComputeDiff(const std::string&   a_fname1,
                const std::string&   a_fname2,
                const std::string&   a_diff_fname,
                std::vector<double>& a_diff_norm_arr )
{
  int ierr;
  const int ngpt = 3;

  SolObj solution_1;
  ierr = readBinaryFile(a_fname1.c_str(), &solution_1);
  if (ierr) return ierr;

  SolObj solution_2;
  ierr = readBinaryFile(a_fname2.c_str(), &solution_2);
  if (ierr) return ierr;

  SolObj solution_2_int;
  createSolObj( solution_1.ndims,
                solution_1.nvars,
                solution_1.size,
                &solution_2_int );

  ierr = InterpolateGlobalnDVar(solution_2_int.size,
                                solution_2_int.u,
                                solution_2.size,
                                solution_2.u,
                                solution_2_int.nvars,
                                solution_2_int.ndims,
                                solution_2_int.is_periodic );
  if (ierr) return ierr;
  deleteSolObj(&solution_2);

  a_diff_norm_arr.resize(3, -1);
  ierr = ComputeDiffNorm( &solution_1,
                          &solution_2_int,
                          a_diff_norm_arr.data() );
  if (ierr) return ierr;

  writeBinary(  solution_2_int.ndims,
                solution_2_int.nvars,
                solution_2_int.size,
                solution_1.x,
                solution_2_int.u,
                a_diff_fname.c_str() );

  deleteSolObj(&solution_1);
  deleteSolObj(&solution_2_int);

  printf("  Norms of the diff between solutions (L1, L2, Linf): %1.6e, %1.6e, %1.6e\n",
          a_diff_norm_arr[0],
          a_diff_norm_arr[1],
          a_diff_norm_arr[2] );

  return 0;
}

int main()
{
  FILE *inputs;
  char op_file_format[50], overwrite[50];

  std::string op_fname_extn = ".bin";
  char op_fname_root_1[50] = "lionocerangotangadinf";
  char op_fname_root_2[50] = "lionocerangotangadinf";
  char diff_fname_root[50] = "diff";

  printf("Reading compare.inp.\n");
  inputs = fopen("compare.inp","r");
  if (!inputs) {
    fprintf(stderr,"Error: File \"compare.inp\" not found.\n");
    return(1);
  } else {
    char word[100];
    fscanf(inputs,"%s",word);
    if (!strcmp(word, "begin")){
      while (strcmp(word, "end")){
        fscanf(inputs,"%s",word);
         if (!strcmp(word, "op_fname_root_1")) {
          fscanf(inputs,"%s" ,op_fname_root_1);
         } else if (!strcmp(word, "op_fname_root_2")) {
          fscanf(inputs,"%s" ,op_fname_root_2);
         } else if (!strcmp(word, "diff_fname_root")) {
          fscanf(inputs,"%s" ,diff_fname_root);
        }
      }
    }
    fclose(inputs);
  }
  printf("Prefix 1: %s\n", op_fname_root_1);
  printf("Prefix 2: %s\n", op_fname_root_2);
  printf("Diff Prefix: %s\n", diff_fname_root);


  int n_iter, file_op_iter = 1;
  double dt;

  printf("Reading solver.inp.\n");
  inputs = fopen("solver.inp","r");
  if (!inputs) {
    fprintf(stderr,"Error: File \"solver.inp\" not found.\n");
    return(1);
  } else {
    char word[100];
    fscanf(inputs,"%s",word);
    if (!strcmp(word, "begin")){
      while (strcmp(word, "end")){
        fscanf(inputs,"%s",word);
         if (!strcmp(word, "op_file_format")) {
          fscanf(inputs,"%s" ,op_file_format);
        } else if (!strcmp(word, "op_overwrite")) {
          fscanf(inputs,"%s" ,overwrite);
        }  else if (!strcmp(word, "n_iter")) {
          fscanf(inputs,"%d" ,&n_iter);
        }  else if (!strcmp(word, "file_op_iter")) {
          fscanf(inputs,"%d" ,&file_op_iter);
        }  else if (!strcmp(word, "dt")) {
          fscanf(inputs,"%f" ,&dt);
        }
      }
    }
    fclose(inputs);
  }

  if (strcmp(op_file_format,"binary") && strcmp(op_file_format,"bin")) {
    printf("Error: solution output needs to be in binary files.\n");
    return(0);
  }

  double t_final = dt * n_iter;
  double dt_file_op = file_op_iter * dt;

  printf("dt = %lf\n", dt);
  printf("n_iter = %d\n", n_iter);
  printf("file_op_iter = %d\n", file_op_iter);
  printf("t_final = %lf\n", t_final);
  printf("dt_file_op = %lf\n", dt_file_op);

  FILE* out;
  out = fopen("compare.out","w");
  if (!strcmp(overwrite,"no")) {
    int counter = 0;
    while(1) {
      /* set filename */
      char counter_str[6];
      sprintf(counter_str, "%05d", counter);
      std::string filename1 = std::string(op_fname_root_1)
                              + "_"
                              + std::string(counter_str)
                              + op_fname_extn;
      std::string filename2 = std::string(op_fname_root_2)
                              + "_"
                              + std::string(counter_str)
                              + op_fname_extn;
      std::string diff_name = std::string(diff_fname_root)
                              + "_"
                              + std::string(counter_str)
                              + op_fname_extn;
      std::vector<double> diff_norm_arr;
      int err = ComputeDiff(filename1, filename2, diff_name, diff_norm_arr);
      if (err) {
        printf("No more files with prefix %s,%s found.\n",
                op_fname_root_1,
                op_fname_root_2);
        break;
      }
      fprintf(out, "%lf  %1.16e  %1.16e  %1.16e\n",
              counter*dt_file_op,
              diff_norm_arr[0],
              diff_norm_arr[1],
              diff_norm_arr[2] );
      counter++;
    }
  } else if (!strcmp(overwrite,"yes")) {
    /* set filename */
    std::string filename1 = std::string(op_fname_root_1)
                            + op_fname_extn;
    std::string filename2 = std::string(op_fname_root_2)
                            + op_fname_extn;
    std::string diff_name = std::string(diff_fname_root)
                            + op_fname_extn;
    std::vector<double> diff_norm_arr;
    int err = ComputeDiff(filename1, filename2, diff_name, diff_norm_arr);
    if (err) {
      printf("Unable to read one or both of specified files (prefix: %s, %s)\n",
              op_fname_root_1,
              op_fname_root_2);
    }
    fprintf(out, "%lf  %1.16e  %1.16e  %1.16e\n",
            t_final,
            diff_norm_arr[0],
            diff_norm_arr[1],
            diff_norm_arr[2] );
  }
  fclose(out);

  return(0);
}

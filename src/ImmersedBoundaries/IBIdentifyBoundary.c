/*! @file IBIdentifyBoundary.c
    @brief Identiy the boundary nodes of the immersed body
    @author Debojyoti Ghosh
*/

#include <stdlib.h>
#include <arrayfunctions.h>
#include <mpivars.h>
#include <immersedboundaries.h>

/*! Count the number of immersed boundary points: boundary points are those
    grid points inside the immersed body that are within stencil-width-distance of
    a grid point outside the body.
*/
static int CountBoundaryPoints(
                                int     a_imax,   /*!< Number of grid points in x */
                                int     a_jmax,   /*!< Number of grid points in y */
                                int     a_kmax,   /*!< Number of grid points in z */
                                int     a_ghosts, /*!< Number of ghost points */
                                double *a_blank /*!< blanking array where entries are zero
                                                     for grid points inside, and one for
                                                     grid points outside. */
                              )
{
  static int dim[_IB_NDIMS_], indexC[_IB_NDIMS_], indexN[_IB_NDIMS_],
             i, j, k, p, q, count;
  dim[0] = a_imax;
  dim[1] = a_jmax;
  dim[2] = a_kmax;

  count = 0;
  for (i = 0; i < a_imax; i++) {
    for (j = 0; j < a_jmax; j++) {
      for (k = 0; k < a_kmax; k++) {
        indexC[0] = i;
        indexC[1] = j;
        indexC[2] = k;
        _ArrayIndex1D_(_IB_NDIMS_,dim,indexC,a_ghosts,p);
        /* if this point is inside the body (0), find out if any */
        /* of the neighboring points are outside (1)              */
        if (!a_blank[p]){

          int g, flag = 0;
          for (g = 1; g <= a_ghosts; g++){

            _ArrayCopy1D_(indexC,indexN,_IB_NDIMS_); indexN[0] += g;
            _ArrayIndex1D_(_IB_NDIMS_,dim,indexN,a_ghosts,q);
            if (a_blank[q])  flag = 1;

            _ArrayCopy1D_(indexC,indexN,_IB_NDIMS_); indexN[0] -= g;
            _ArrayIndex1D_(_IB_NDIMS_,dim,indexN,a_ghosts,q);
            if (a_blank[q])  flag = 1;

            _ArrayCopy1D_(indexC,indexN,_IB_NDIMS_); indexN[1] += g;
            _ArrayIndex1D_(_IB_NDIMS_,dim,indexN,a_ghosts,q);
            if (a_blank[q])  flag = 1;

            _ArrayCopy1D_(indexC,indexN,_IB_NDIMS_); indexN[1] -= g;
            _ArrayIndex1D_(_IB_NDIMS_,dim,indexN,a_ghosts,q);
            if (a_blank[q])  flag = 1;

            _ArrayCopy1D_(indexC,indexN,_IB_NDIMS_); indexN[2] += g;
            _ArrayIndex1D_(_IB_NDIMS_,dim,indexN,a_ghosts,q);
            if (a_blank[q])  flag = 1;

            _ArrayCopy1D_(indexC,indexN,_IB_NDIMS_); indexN[2] -= g;
            _ArrayIndex1D_(_IB_NDIMS_,dim,indexN,a_ghosts,q);
            if (a_blank[q])  flag = 1;

          }
          if (flag)   count++;
        }
      }
    }
  }
  return(count);
}

/*! Set the indices of the immersed boundary points.*/
static int SetBoundaryPoints(
                                int     a_imax,     /*!< Number of grid points in x */
                                int     a_jmax,     /*!< Number of grid points in y */
                                int     a_kmax,     /*!< Number of grid points in z */
                                int     a_ghosts,   /*!< Number of ghost points */
                                double  *a_blank,   /*!< blanking array where entries are zero
                                                     for grid points inside, and one for
                                                     grid points outside. */
                                void    *a_b        /*!< Array of immersed boundary points of type #IBNode */
                            )
{
  IBNode *boundary = (IBNode*) a_b;
  static int dim[_IB_NDIMS_], indexC[_IB_NDIMS_], indexN[_IB_NDIMS_],
             i, j, k, p, q, count;
  dim[0] = a_imax;
  dim[1] = a_jmax;
  dim[2] = a_kmax;

  count = 0;
  for (i = 0; i < a_imax; i++) {
    for (j = 0; j < a_jmax; j++) {
      for (k = 0; k < a_kmax; k++) {
        indexC[0] = i;
        indexC[1] = j;
        indexC[2] = k;
        _ArrayIndex1D_(_IB_NDIMS_,dim,indexC,a_ghosts,p);
        /* if this point is inside the body (0), find out if any */
        /* of the neighboring points are outside (1)              */
        if (!a_blank[p]){

          int g, flag = 0;
          for (g = 1; g <= a_ghosts; g++){

            _ArrayCopy1D_(indexC,indexN,_IB_NDIMS_); indexN[0] += g;
            _ArrayIndex1D_(_IB_NDIMS_,dim,indexN,a_ghosts,q);
            if (a_blank[q])  flag = 1;

            _ArrayCopy1D_(indexC,indexN,_IB_NDIMS_); indexN[0] -= g;
            _ArrayIndex1D_(_IB_NDIMS_,dim,indexN,a_ghosts,q);
            if (a_blank[q])  flag = 1;

            _ArrayCopy1D_(indexC,indexN,_IB_NDIMS_); indexN[1] += g;
            _ArrayIndex1D_(_IB_NDIMS_,dim,indexN,a_ghosts,q);
            if (a_blank[q])  flag = 1;

            _ArrayCopy1D_(indexC,indexN,_IB_NDIMS_); indexN[1] -= g;
            _ArrayIndex1D_(_IB_NDIMS_,dim,indexN,a_ghosts,q);
            if (a_blank[q])  flag = 1;

            _ArrayCopy1D_(indexC,indexN,_IB_NDIMS_); indexN[2] += g;
            _ArrayIndex1D_(_IB_NDIMS_,dim,indexN,a_ghosts,q);
            if (a_blank[q])  flag = 1;

            _ArrayCopy1D_(indexC,indexN,_IB_NDIMS_); indexN[2] -= g;
            _ArrayIndex1D_(_IB_NDIMS_,dim,indexN,a_ghosts,q);
            if (a_blank[q])  flag = 1;

          }
          if (flag) {
            boundary[count].m_i = i;
            boundary[count].m_j = j;
            boundary[count].m_k = k;
            boundary[count].m_p = p;
            count++;
          }
        }
      }
    }
  }
  return(count);
}

/*! Identify the immersed boundary points: an immersed boundary point is any grid point
    inside the immersed body that is within stencil-width-distance of a grid point outside
    the immersed body. This function does the following:
    + count the number of immersed boundary points.
    + allocate the array of immersed boundary points and set their indices.
*/
int IBIdentifyBoundary(
                        void   *a_ib,     /*!< Immersed boundary object of type #ImmersedBoundary */
                        void   *a_m,      /*!< MPI object of type #MPIVariables */
                        int    *a_dim_l,  /*!< local dimensions */
                        int    a_ghosts,  /*!< number of ghost points */
                        double *a_blank /*!< Blanking array: for grid points within the
                                             immersed body, this value will be set to 0 */
                      )
{
  ImmersedBoundary  *IB     = (ImmersedBoundary*) a_ib;
  MPIVariables      *mpi    = (MPIVariables*) a_m;
  Body3D            *body   = IB->m_body;

  int imax = a_dim_l[0],
      jmax = a_dim_l[1],
      kmax = a_dim_l[2];

  int n_boundary_nodes = CountBoundaryPoints(imax,jmax,kmax,a_ghosts,a_blank);
  IB->m_n_boundary_nodes = n_boundary_nodes;
  if (n_boundary_nodes == 0) IB->m_boundary = NULL;
  else {
    IB->m_boundary = (IBNode*) calloc (n_boundary_nodes, sizeof(IBNode));
    int check = SetBoundaryPoints(imax,jmax,kmax,a_ghosts,a_blank,IB->m_boundary);
    if (check != n_boundary_nodes) {
      fprintf(stderr,"Error in IBIdentifyBoundary(): Inconsistency encountered when setting boundary indices. ");
      fprintf(stderr,"on rank %d.\n",mpi->m_rank);
      fprintf(stderr,"SetBoundaryPoints() returned %d, while n_boundary_nodes is %d.\n",check,n_boundary_nodes);
    }
  }

  return(n_boundary_nodes);
}

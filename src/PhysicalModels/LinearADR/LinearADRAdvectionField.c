/*! @file LinearADRAdvectionField.c
    @author Debojyoti Ghosh
    @brief Contains the function to read in a non-constant advection field
*/

#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <basic.h>
#include <common.h>
#include <arrayfunctions.h>
#include <io.h>
#include <physicalmodels/linearadr.h>
#include <hypar.h>
#include <mpivars.h>

/*! Set the advection field over the domain - reads the spatially-varying
 * advection data from a file, if available. The array to store the field
 * \b must already be allocated.
 *
 * For a single simulation, it reads in the data from a file. The data
 * must have the same grid dimensions as the solution.
 *
 * For an ensemble simulation, it reads in the advection field from the
 * file corresponding to the \b idx (index of this simulation). The
 * data in this file must have the same grid dimensions as the domain
 * it is being read in for.
 *
 * For a sparse grids simulation, it reads in the advection field from
 * the file with data that has the grid dimensions as the full grid. The
 * field on the current grid is obtained by interpolation.
*/
int LinearADRAdvectionField(void *s,      /*!< Solver object of type #HyPar */
                            void *m,      /*!< MPI object of type #MPIVariables*/
                            int  idx,     /*!< Index of this simulation */
                            int  nsims,   /*!< Total number of simulations */
                            int *dim_adv  /*!< Grid dimensions of the advection field stored in file */
                           )
{
  HyPar         *solver = (HyPar*) s;
  MPIVariables  *mpi    = (MPIVariables*) m;
  LinearADR     *param  = (LinearADR*) solver->m_physics;
  double        *adv    = param->m_a;

  /* sanity checks */
  if (param->m_constant_advection == 1) {
    fprintf(stderr,"Error in LinearADRAdvectionField(): param->m_constant_advection == 1!\n");
    return(1);
  }
  if (!strcmp(param->m_adv_filename,"none")) {
    fprintf(stderr,"Error in LinearADRAdvectionField(): param->m_adv_filename is \"none\"!\n");
    return(1);
  }

  int d, done, flag;
  int *dim = solver->m_dim_local;
  int ghosts = solver->m_ghosts;

  /* number of advection field components at a grid point is
   * number of spatial dimensions times number of solution
   * components */
  int adv_nvar = solver->m_ndims * solver->m_nvars;
  if (param->m_adv_arr_size != adv_nvar*solver->m_npoints_local_wghosts) {
    fprintf(stderr,"Error in LinearADRAdvectionField(): Incorrect adv_arr_size\n");
    return(1);
  }

  char fname_root[_MAX_STRING_SIZE_];
  strcpy(fname_root, param->m_adv_filename);
  if (idx >= 0) {
    if (nsims > 1) {
      char index[_MAX_STRING_SIZE_];
      GetStringFromInteger(idx, index, (int)log10(nsims)+1);
      strcat(fname_root, "_");
      strcat(fname_root, index);
    }
  }

  /* read spatially-varying advection field from provided file */
  if (dim_adv == NULL) {
    int ierr = ReadArray( solver->m_ndims,
                          adv_nvar,
                          solver->m_dim_global,
                          solver->m_dim_local,
                          solver->m_ghosts,
                          solver,
                          mpi,
                          NULL,
                          adv,
                          fname_root,
                          &flag);
    if (ierr) {
      fprintf(stderr,"Error in LinearADRAdvectionField()\n");
      fprintf(stderr,"  ReadArray() returned error!\n");
      return ierr;
    }
  } else {
    int ierr = ReadArraywInterp(solver->m_ndims,
                                adv_nvar,
                                solver->m_dim_global,
                                solver->m_dim_local,
                                dim_adv,
                                solver->m_ghosts,
                                solver,
                                mpi,
                                NULL,
                                adv,
                                fname_root,
                                &flag);
    if (ierr) {
      fprintf(stderr,"Error in LinearADRAdvectionField()\n");
      fprintf(stderr,"  ReadArraywInterp() returned error!\n");
      return ierr;
    }
  }

  if (!flag) {
    /* if advection file not available, set it to zero */
    _ArraySetValue_(adv,solver->m_npoints_local_wghosts*adv_nvar,0.0);
  }

  /* if parallel, exchange MPI-boundary ghost point data */
  IERR MPIExchangeBoundariesnD( solver->m_ndims,
                                adv_nvar,
                                solver->m_dim_local,
                                solver->m_ghosts,
                                mpi,
                                adv );

  /* Along each dimension:
     - If boundaries are periodic, then set ghost points accordingly if
       number of MPI ranks along that dimension is 1 (if more than 1,
       MPIExchangeBoundariesnD() called above has taken care of it).
     - Else, extrapolate at the boundaries. */
  int indexb[solver->m_ndims],
      indexi[solver->m_ndims],
      bounds[solver->m_ndims],
      offset[solver->m_ndims];
  for (d = 0; d < solver->m_ndims; d++) {
    if (solver->m_is_periodic[d] && (mpi->m_iproc[d] == 1)) {
      _ArrayCopy1D_(dim,bounds,solver->m_ndims); bounds[d] = ghosts;
      /* left boundary */
      done = 0; _ArraySetValue_(indexb,solver->m_ndims,0);
      _ArraySetValue_(offset,solver->m_ndims,0); offset[d] = -ghosts;
      while (!done) {
        _ArrayCopy1D_(indexb,indexi,solver->m_ndims); indexi[d] = indexb[d] + dim[d] - ghosts;
        int p1; _ArrayIndex1DWO_(solver->m_ndims,dim,indexb,offset,ghosts,p1);
        int p2; _ArrayIndex1D_  (solver->m_ndims,dim,indexi,ghosts,p2);
        _ArrayCopy1D_((adv+p2*adv_nvar), (adv+p1*adv_nvar), adv_nvar);
        _ArrayIncrementIndex_(solver->m_ndims,bounds,indexb,done);
      }
      /* right boundary */
      done = 0; _ArraySetValue_(indexb,solver->m_ndims,0);
      _ArraySetValue_(offset,solver->m_ndims,0); offset[d] = dim[d];
      while (!done) {
        _ArrayCopy1D_(indexb,indexi,solver->m_ndims);
        int p1; _ArrayIndex1DWO_(solver->m_ndims,dim,indexb,offset,ghosts,p1);
        int p2; _ArrayIndex1D_  (solver->m_ndims,dim,indexi,ghosts,p2);
        _ArrayCopy1D_((adv+p2*adv_nvar), (adv+p1*adv_nvar), adv_nvar);
        _ArrayIncrementIndex_(solver->m_ndims,bounds,indexb,done);
      }
    } else {
      /* left boundary */
      if (!mpi->m_ip[d]) {
        _ArrayCopy1D_(dim,bounds,solver->m_ndims); bounds[d] = ghosts;
        _ArraySetValue_(offset,solver->m_ndims,0); offset[d] = -ghosts;
        done = 0; _ArraySetValue_(indexb,solver->m_ndims,0);
        while (!done) {
          _ArrayCopy1D_(indexb,indexi,solver->m_ndims); indexi[d] = ghosts-1-indexb[d];
          int p1; _ArrayIndex1DWO_(solver->m_ndims,dim,indexb,offset,ghosts,p1);
          int p2; _ArrayIndex1D_  (solver->m_ndims,dim,indexi,ghosts,p2);
          _ArrayCopy1D_((adv+p2*adv_nvar), (adv+p1*adv_nvar), adv_nvar);
          _ArrayIncrementIndex_(solver->m_ndims,bounds,indexb,done);
        }
      }
      /* right boundary */
      if (mpi->m_ip[d] == mpi->m_iproc[d]-1) {
        _ArrayCopy1D_(dim,bounds,solver->m_ndims); bounds[d] = ghosts;
        _ArraySetValue_(offset,solver->m_ndims,0); offset[d] = dim[d];
        done = 0; _ArraySetValue_(indexb,solver->m_ndims,0);
        while (!done) {
          _ArrayCopy1D_(indexb,indexi,solver->m_ndims); indexi[d] = dim[d]-1-indexb[d];
          int p1; _ArrayIndex1DWO_(solver->m_ndims,dim,indexb,offset,ghosts,p1);
          int p2; _ArrayIndex1D_  (solver->m_ndims,dim,indexi,ghosts,p2);
          _ArrayCopy1D_((adv+p2*adv_nvar), (adv+p1*adv_nvar), adv_nvar);
          _ArrayIncrementIndex_(solver->m_ndims,bounds,indexb,done);
        }
      }
    }
  }

  return(0);
}

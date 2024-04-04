/*! @file ShallowWater2DTopography.c
    @author Debojyoti Ghosh
    @brief Contains the function to set the bottom topography
*/

#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <basic.h>
#include <common.h>
#include <arrayfunctions.h>
#include <io.h>
#include <physicalmodels/shallowwater2d.h>
#include <hypar.h>
#include <mpivars.h>

/*! Set the bottom topography over the domain - reads the topography
    data from a file, if available, else sets it to a constant
*/
int ShallowWater2DTopography(
                              void *s,       /*!< Solver object of type #HyPar */
                              void *m,       /*!< MPI object of type #MPIVariables*/
                              int  idx,      /*!< Index of this simulation */
                              int  nsims,    /*!< Total number of simulations */
                              int  *dim_topo /*!< Grid dimensions of the advection field stored in file */
                            )
{
  HyPar          *solver = (HyPar*) s;
  MPIVariables   *mpi    = (MPIVariables*) m;
  ShallowWater2D *param  = (ShallowWater2D*) solver->physics;
  double         *S      = param->b;
  int            d, done, *dim = solver->dim_local,
                 ghosts = solver->ghosts;
  _DECLARE_IERR_;

  char fname_root[_MAX_STRING_SIZE_] = "topography";
  if (idx >= 0) {
    if (nsims > 1) {
      char index[_MAX_STRING_SIZE_];
      GetStringFromInteger(idx, index, (int)log10(nsims)+1);
      strcat(fname_root, "_");
      strcat(fname_root, index);
    }
  }

  /* read topography from provided file, if available */
  if (dim_topo == NULL) {
    int ierr = ReadArray( solver->ndims,
                          1,
                          solver->dim_global,
                          solver->dim_local,
                          solver->ghosts,
                          solver,
                          mpi,
                          NULL,
                          S,
                          fname_root,
                          &param->topo_flag);
    if (ierr) {
      fprintf(stderr,"Error in ShallowWater2DTopography()\n");
      fprintf(stderr,"  ReadArray() returned error!\n");
      return ierr;
    }
  } else {
    int ierr = ReadArraywInterp(  solver->ndims,
                                  1,
                                  solver->dim_global,
                                  solver->dim_local,
                                  dim_topo,
                                  solver->ghosts,
                                  solver,
                                  mpi,
                                  NULL,
                                  S,
                                  fname_root,
                                  &param->topo_flag);
    if (ierr) {
      fprintf(stderr,"Error in ShallowWater2DTopography()\n");
      fprintf(stderr,"  ReadArraywInterp() returned error!\n");
      return ierr;
    }
  }

  if (!param->topo_flag) {
    /* if topography file not available, set it to zero */
    _ArraySetValue_(S,solver->npoints_local_wghosts,0.0);
  }

  /* if parallel, exchange MPI-boundary ghost point data */
  IERR MPIExchangeBoundariesnD(_MODEL_NDIMS_,1,solver->dim_local,
                                solver->ghosts,mpi,S); CHECKERR(ierr);


  if (param->bt_type) {
    /* if topography is periodic, then the overall problem must also be periodic
       (i.e. boundary conditions will be specified as periodic). Hence,
       MPIExchangeBoundariesnD() will take care of setting the ghosts points
       for multi-processor simulations. For single processor, set the ghost
       points accordingly. */
    int indexb[_MODEL_NDIMS_], indexi[_MODEL_NDIMS_], bounds[_MODEL_NDIMS_],
        offset[_MODEL_NDIMS_];
    for (d = 0; d < _MODEL_NDIMS_; d++) {
      if (mpi->iproc[d] == 1) {
        _ArrayCopy1D_(dim,bounds,_MODEL_NDIMS_); bounds[d] = ghosts;
        /* left boundary */
        done = 0; _ArraySetValue_(indexb,_MODEL_NDIMS_,0);
        _ArraySetValue_(offset,_MODEL_NDIMS_,0); offset[d] = -ghosts;
        while (!done) {
          _ArrayCopy1D_(indexb,indexi,_MODEL_NDIMS_); indexi[d] = indexb[d] + dim[d] - ghosts;
          int p1; _ArrayIndex1DWO_(_MODEL_NDIMS_,dim,indexb,offset,ghosts,p1);
          int p2; _ArrayIndex1D_  (_MODEL_NDIMS_,dim,indexi,ghosts,p2);
          S[p1] = S[p2];
          _ArrayIncrementIndex_(_MODEL_NDIMS_,bounds,indexb,done);
        }
        /* right boundary */
        done = 0; _ArraySetValue_(indexb,_MODEL_NDIMS_,0);
        _ArraySetValue_(offset,_MODEL_NDIMS_,0); offset[d] = dim[d];
        while (!done) {
          _ArrayCopy1D_(indexb,indexi,_MODEL_NDIMS_);
          int p1; _ArrayIndex1DWO_(_MODEL_NDIMS_,dim,indexb,offset,ghosts,p1);
          int p2; _ArrayIndex1D_  (_MODEL_NDIMS_,dim,indexi,ghosts,p2);
          S[p1] = S[p2];
          _ArrayIncrementIndex_(_MODEL_NDIMS_,bounds,indexb,done);
        }
      }
    }
  } else {
    /* if topography is not periodic, extrapolate it at the boundaries */
    int indexb[_MODEL_NDIMS_], indexi[_MODEL_NDIMS_], bounds[_MODEL_NDIMS_],
        offset[_MODEL_NDIMS_];
    for (d = 0; d < _MODEL_NDIMS_; d++) {
      /* left boundary */
      if (!mpi->ip[d]) {
        _ArrayCopy1D_(dim,bounds,_MODEL_NDIMS_); bounds[d] = ghosts;
        _ArraySetValue_(offset,_MODEL_NDIMS_,0); offset[d] = -ghosts;
        done = 0; _ArraySetValue_(indexb,_MODEL_NDIMS_,0);
        while (!done) {
          _ArrayCopy1D_(indexb,indexi,_MODEL_NDIMS_); indexi[d] = ghosts-1-indexb[d];
          int p1; _ArrayIndex1DWO_(_MODEL_NDIMS_,dim,indexb,offset,ghosts,p1);
          int p2; _ArrayIndex1D_  (_MODEL_NDIMS_,dim,indexi,ghosts,p2);
          S[p1] = S[p2];
          _ArrayIncrementIndex_(_MODEL_NDIMS_,bounds,indexb,done);
        }
      }
      /* right boundary */
      if (mpi->ip[d] == mpi->iproc[d]-1) {
        _ArrayCopy1D_(dim,bounds,_MODEL_NDIMS_); bounds[d] = ghosts;
        _ArraySetValue_(offset,_MODEL_NDIMS_,0); offset[d] = dim[d];
        done = 0; _ArraySetValue_(indexb,_MODEL_NDIMS_,0);
        while (!done) {
          _ArrayCopy1D_(indexb,indexi,_MODEL_NDIMS_); indexi[d] = dim[d]-1-indexb[d];
          int p1; _ArrayIndex1DWO_(_MODEL_NDIMS_,dim,indexb,offset,ghosts,p1);
          int p2; _ArrayIndex1D_  (_MODEL_NDIMS_,dim,indexi,ghosts,p2);
          S[p1] = S[p2];
          _ArrayIncrementIndex_(_MODEL_NDIMS_,bounds,indexb,done);
        }
      }
    }
  }

  return(0);
}

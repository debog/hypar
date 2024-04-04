/*! @file PetscComputePreconMatImpl.cpp
    @brief Contains the function to assemble the preconditioning matrix
    @author Debojyoti Ghosh
*/

#ifdef with_petsc

#include <stdio.h>
#include <arrayfunctions.h>
#include <simulation_object.h>
#include <mpivars_cpp.h>
#include <petscinterface.h>

#undef __FUNCT__
#define __FUNCT__ "PetscComputePreconMatImpl"
/*!
  Compute and assemble the preconditioning matrix for the implicit time integration
  of the governing equations: The ODE, obtained after discretizing the governing PDE in space,
  is expressed as follows:
  \f{equation}{
    \frac {d{\bf U}}{dt} = {\bf L}\left({\bf U}\right)
    \Rightarrow \frac {d{\bf U}}{dt} - {\bf L}\left({\bf U}\right) = 0,
  \f}
  where \f${\bf L}\f$ is the spatially discretized right-hand-side, and \f${\bf U}\f$
  represents the entire solution vector.

  The Jacobian is thus given by:
  \f{equation}{
    {\bf J} = \left[\alpha{\bf I} - \frac {\partial {\bf L}} {\partial {\bf U}} \right]
  \f}
  where \f$\alpha\f$ is the shift coefficient (#PETScContext::shift) of the time integration method.

  Let \f$\mathcal{D}_{\rm hyp}\f$ and \f$\mathcal{D}_{\rm par}\f$ represent the spatial discretization methods for
  the hyperbolic and parabolic terms. Thus, the Jacobian can be written as
  follows:
  \f{equation}{
    {\bf J} = \left[\alpha{\bf I} - \left(\mathcal{D}_{\rm hyp}\left\{\frac {\partial {\bf F}_{\rm hyp}} {\partial {\bf u}}\right\} + \mathcal{D}_{\rm par}\left\{\frac {\partial {\bf F}_{\rm par}} {\partial {\bf u}}\right\}\right)\right]
  \f}
  The preconditioning matrix is usually a close approximation of the actual Jacobian matrix, where the actual
  Jacobian may be too expensive to evaluate and assemble. In this function, the preconditioning matrix is
  the following approximation of the actual Jacobian:
  \f{equation}{
    {\bf J}_p = \left[\alpha{\bf I} - \left(\mathcal{D}^{\left(l\right)}_{\rm hyp}\left\{\frac {\partial {\bf F}_{\rm hyp}} {\partial {\bf u}}\right\} + \mathcal{D}^{\left(l\right)}_{\rm par}\left\{\frac {\partial {\bf F}_{\rm par}} {\partial {\bf u}}\right\}\right) \right] \approx {\bf J},
  \f}
  where \f$\mathcal{D}^{\left(l\right)}_{\rm hyp,par}\f$ represent lower order discretizations of the hyperbolic and parabolic terms. The matrix \f${\bf J}_p\f$
  is provided to the preconditioner.

  Note that the specific physics model provides the following functions:
  + #HyPar::JFunction computes \f$\partial {\bf F}_{\rm hyp}/ \partial {\bf u}\f$ at a grid point.
  + #HyPar::KFunction computes \f$\partial {\bf F}_{\rm par}/ \partial {\bf u}\f$ at a grid point.

  Currently, this function doesn't include the source term.

  + See https://petsc.org/release/docs/manualpages/PC/index.html for more information on PETSc preconditioners.
  + All functions and variables whose names start with Vec, Mat, PC, KSP, SNES, and TS are defined by PETSc. Refer to
    the PETSc documentation (https://petsc.org/release/docs/). Usually, googling with the function
    or variable name yields the specific doc page dealing with that function/variable.
*/
int PetscComputePreconMatImpl(  Mat Pmat,   /*!< Preconditioning matrix to construct */
                                Vec Y,      /*!< Solution vector */
                                void *ctxt  /*!< Application context */ )
{
  PetscErrorCode ierr;
  PETScContext* context = (PETScContext*) ctxt;
  SimulationObject* sim = (SimulationObject*) context->simobj;
  int nsims = context->nsims;

  PetscFunctionBegin;
  /* initialize preconditioning matrix to zero */
  MatZeroEntries(Pmat);

  /* copy solution from PETSc vector */
  for (int ns = 0; ns < nsims; ns++) {

    TransferVecFromPETSc( sim[ns].solver.u,
                          Y,
                          context,
                          ns,
                          context->offsets[ns]);

    HyPar* solver( &(sim[ns].solver) );
    MPIVariables* mpi( &(sim[ns].mpi) );

    int ndims = solver->ndims,
        nvars = solver->nvars,
        npoints = solver->npoints_local,
        ghosts = solver->ghosts,
        *dim = solver->dim_local,
        *isPeriodic = solver->isPeriodic,
        *points = context->points[ns],
        index[ndims],indexL[ndims],indexR[ndims],
        rows[nvars],cols[nvars];

    double *u = solver->u,
           *iblank = solver->iblank,
           dxinv, values[nvars*nvars];

    /* apply boundary conditions and exchange data over MPI interfaces */
    solver->ApplyBoundaryConditions(solver,mpi,u,NULL,context->waqt);
    MPIExchangeBoundariesnD(ndims,nvars,dim,ghosts,mpi,u);

    /* loop through all computational points */
    for (int n = 0; n < npoints; n++) {
      int *this_point = points + n*(ndims+1);
      int p = this_point[ndims];
      int index[ndims]; _ArrayCopy1D_(this_point,index,ndims);

      double iblank = solver->iblank[p];

      /* compute the contributions from the hyperbolic flux derivatives along each dimension */
      if (solver->JFunction) {
        for (int dir = 0; dir < ndims; dir++) {

          /* compute indices and global 1D indices for left and right neighbors */
          _ArrayCopy1D_(index,indexL,ndims); indexL[dir]--;
          int pL;  _ArrayIndex1D_(ndims,dim,indexL,ghosts,pL);

          _ArrayCopy1D_(index,indexR,ndims); indexR[dir]++;
          int pR;  _ArrayIndex1D_(ndims,dim,indexR,ghosts,pR);

          int pg, pgL, pgR;
          pg  = (int) context->globalDOF[ns][p];
          pgL = (int) context->globalDOF[ns][pL];
          pgR = (int) context->globalDOF[ns][pR];

          /* Retrieve 1/delta-x at this grid point */
          _GetCoordinate_(dir,index[dir],dim,ghosts,solver->dxinv,dxinv);

          /* diagonal element */
          for (int v=0; v<nvars; v++) { rows[v] = nvars*pg + v; cols[v] = nvars*pg + v; }
          solver->JFunction(values,(u+nvars*p),solver->physics,dir,nvars,0);
          _ArrayScale1D_(values,(dxinv*iblank),(nvars*nvars));
          MatSetValues(Pmat,nvars,rows,nvars,cols,values,ADD_VALUES);

          /* left neighbor */
          if (pgL >= 0) {
            for (int v=0; v<nvars; v++) { rows[v] = nvars*pg + v; cols[v] = nvars*pgL + v; }
            solver->JFunction(values,(u+nvars*pL),solver->physics,dir,nvars,1);
            _ArrayScale1D_(values,(-dxinv*iblank),(nvars*nvars));
            MatSetValues(Pmat,nvars,rows,nvars,cols,values,ADD_VALUES);
          }

          /* right neighbor */
          if (pgR >= 0) {
            for (int v=0; v<nvars; v++) { rows[v] = nvars*pg + v; cols[v] = nvars*pgR + v; }
            solver->JFunction(values,(u+nvars*pR),solver->physics,dir,nvars,-1);
            _ArrayScale1D_(values,(-dxinv*iblank),(nvars*nvars));
            MatSetValues(Pmat,nvars,rows,nvars,cols,values,ADD_VALUES);
          }
        }
      }

      /* compute the contributions from the parabolic term derivatives along each dimension */
      if (solver->KFunction) {
        for (int dir = 0; dir < ndims; dir++) {

          /* compute indices and global 1D indices for left and right neighbors */
          _ArrayCopy1D_(index,indexL,ndims); indexL[dir]--;
          int pL;  _ArrayIndex1D_(ndims,dim,indexL,ghosts,pL);

          _ArrayCopy1D_(index,indexR,ndims); indexR[dir]++;
          int pR;  _ArrayIndex1D_(ndims,dim,indexR,ghosts,pR);

          int pg, pgL, pgR;
          pg  = (int) context->globalDOF[ns][p];
          pgL = (int) context->globalDOF[ns][pL];
          pgR = (int) context->globalDOF[ns][pR];

          /* Retrieve 1/delta-x at this grid point */
          _GetCoordinate_(dir,index[dir],dim,ghosts,solver->dxinv,dxinv);

          /* diagonal element */
          for (int v=0; v<nvars; v++) { rows[v] = nvars*pg + v; cols[v] = nvars*pg + v; }
          solver->KFunction(values,(u+nvars*p),solver->physics,dir,nvars);
          _ArrayScale1D_(values,(-2*dxinv*dxinv*iblank),(nvars*nvars));
          MatSetValues(Pmat,nvars,rows,nvars,cols,values,ADD_VALUES);

          /* left neighbor */
          if (pgL >= 0) {
            for (int v=0; v<nvars; v++) { rows[v] = nvars*pg + v; cols[v] = nvars*pgL + v; }
            solver->KFunction(values,(u+nvars*pL),solver->physics,dir,nvars);
            _ArrayScale1D_(values,(dxinv*dxinv*iblank),(nvars*nvars));
            MatSetValues(Pmat,nvars,rows,nvars,cols,values,ADD_VALUES);
          }

          /* right neighbor */
          if (pgR >= 0) {
            for (int v=0; v<nvars; v++) { rows[v] = nvars*pg + v; cols[v] = nvars*pgR + v; }
            solver->KFunction(values,(u+nvars*pR),solver->physics,dir,nvars);
            _ArrayScale1D_(values,(dxinv*dxinv*iblank),(nvars*nvars));
            MatSetValues(Pmat,nvars,rows,nvars,cols,values,ADD_VALUES);
          }
        }
      }
    }
  }

  MatAssemblyBegin(Pmat,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd  (Pmat,MAT_FINAL_ASSEMBLY);

  MatShift(Pmat,context->shift);
  PetscFunctionReturn(0);
}

#endif

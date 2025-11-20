/*! @file PetscJacobianMatNonzeroEntriesImpl.cpp
    @brief Contains the function to set nonzero entries of the Jacobian matrix
    @author Debojyoti Ghosh
*/

#ifdef with_petsc

#include <stdio.h>
#include <vector>
#include <arrayfunctions.h>
#include <simulation_object.h>
#include <mpivars_cpp.h>
#include <petscinterface.h>

#undef __FUNCT__
#define __FUNCT__ "PetscJacobianMatNonzeroEntriesImpl"
/*!
  Set non-zero entries of the Jacobian matrix for the implicit time integration
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

  All functions and variables whose names start with Vec, Mat, PC, KSP, SNES, and TS are defined by PETSc. Refer to
  the PETSc documentation (https://petsc.org/release/docs/). Usually, googling with the function
  or variable name yields the specific doc page dealing with that function/variable.
*/
int PetscJacobianMatNonzeroEntriesImpl( Mat a_Amat,   /*!< Matrix */
                                        int a_width,  /*!< Stencil width */
                                        void *a_ctxt  /*!< Application context */ )
{
  PETScContext* context = (PETScContext*) a_ctxt;
  SimulationObject* sim = (SimulationObject*) context->m_simobj;

  PetscFunctionBegin;
  int nsims = context->m_nsims;
  /* initialize matrix to zero */
  MatZeroEntries(a_Amat);

  for (int ns = 0; ns < nsims; ns++) {

    HyPar* solver( &(sim[ns].solver) );
    MPIVariables* mpi( &(sim[ns].mpi) );

    int ndims = solver->m_ndims,
        nvars = solver->m_nvars,
        npoints = solver->m_npoints_local,
        ghosts = solver->m_ghosts,
        *dim = solver->m_dim_local,
        *points = context->m_points[ns],
        index[ndims],indexL[ndims],indexR[ndims],
        rows[nvars],cols[nvars];

    std::vector<double> values(nvars*nvars, 1.0);

    /* loop through all computational points */
    for (int n = 0; n < npoints; n++) {

      int *this_point = points + n*(ndims+1);
      int p = this_point[ndims];
      int index[ndims]; _ArrayCopy1D_(this_point,index,ndims);

      for (int dir = 0; dir < ndims; dir++) {

        int pg = (int) context->m_globalDOF[ns][p];
        /* diagonal element */
        for (int v=0; v<nvars; v++) { rows[v] = nvars*pg + v; cols[v] = nvars*pg + v; }
        MatSetValues(a_Amat,nvars,rows,nvars,cols,values.data(),ADD_VALUES);

        for (int d = 1; d <= a_width; d++) {

          /* left neighbor */
          _ArrayCopy1D_(index,indexL,ndims);
          indexL[dir] -= d;
          int pL;  _ArrayIndex1D_(ndims,dim,indexL,ghosts,pL);
          int pgL = (int) context->m_globalDOF[ns][pL];
          if (pgL >= 0) {
            for (int v=0; v<nvars; v++) { rows[v] = nvars*pg + v; cols[v] = nvars*pgL + v; }
            MatSetValues(a_Amat,nvars,rows,nvars,cols,values.data(),ADD_VALUES);
          }

          _ArrayCopy1D_(index,indexR,ndims);
          indexR[dir] += d;
          int pR;  _ArrayIndex1D_(ndims,dim,indexR,ghosts,pR);
          int pgR = (int) context->m_globalDOF[ns][pR];
          /* right neighbor */
          if (pgR >= 0) {
            for (int v=0; v<nvars; v++) { rows[v] = nvars*pg + v; cols[v] = nvars*pgR + v; }
            MatSetValues(a_Amat,nvars,rows,nvars,cols,values.data(),ADD_VALUES);
          }

        }

      }
    }
  }

  MatAssemblyBegin(a_Amat,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd  (a_Amat,MAT_FINAL_ASSEMBLY);

  PetscFunctionReturn(0);
}

#endif

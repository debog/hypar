/*! @file InitializePhysics.c
    @author Debojyoti Ghosh
    @brief Initialize the physical model
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <basic.h>
#include <bandedmatrix.h>
#include <interpolation.h>
#include <mpivars.h>
#include <simulation_object.h>

/* include header files for each physical model */
#include <physicalmodels/linearadr.h>
#include <physicalmodels/burgers.h>
#include <physicalmodels/fpdoublewell.h>
#include <physicalmodels/fppowersystem.h>
#include <physicalmodels/fppowersystem1bus.h>
#include <physicalmodels/fppowersystem3bus.h>
#include <physicalmodels/euler1d.h>
#include <physicalmodels/euler2d.h>
#include <physicalmodels/navierstokes2d.h>
#include <physicalmodels/navierstokes3d.h>
#include <physicalmodels/numa2d.h>
#include <physicalmodels/numa3d.h>
#include <physicalmodels/shallowwater1d.h>
#include <physicalmodels/shallowwater2d.h>
#include <physicalmodels/vlasov.h>

/*! Initialize the physical model for a simulation: Depending on the
    physical model specified, this function calls the initialization
    function for that physical model. The latter is responsible for
    setting all the physics-specific functions that are required
    by the model.
*/
int InitializePhysics(  void  *a_s,   /*!< Array of simulation objects of type #SimulationObject */
                        int   a_nsims /*!< number of simulation objects */
                     )
{
  SimulationObject *sim = (SimulationObject*) a_s;
  int ns;
  _DECLARE_IERR_;

  if (a_nsims == 0) return 0;

  if (!sim[0].mpi.m_rank) {
    printf("Initializing physics. Model = \"%a_s\"\n",sim[0].solver.m_model);
  }

  for (ns = 0; ns < a_nsims; ns++) {

    HyPar        *solver   = &(sim[ns].solver);
    MPIVariables *mpi      = &(sim[ns].mpi);

    /* Initialize physics-specific functions to NULL */
    solver->ComputeCFL            = NULL;
    solver->ComputeDiffNumber     = NULL;
    solver->FFunction             = NULL;
    solver->dFFunction            = NULL;
    solver->FdFFunction           = NULL;
    solver->GFunction             = NULL;
    solver->HFunction             = NULL;
    solver->SFunction             = NULL;
    solver->UFunction             = NULL;
    solver->JFunction             = NULL;
    solver->KFunction             = NULL;
    solver->Upwind                = NULL;
    solver->UpwinddF              = NULL;
    solver->UpwindFdF             = NULL;
    solver->PreStage              = NULL;
    solver->PostStage             = NULL;
    solver->PreStep               = NULL;
    solver->PostStep              = NULL;
    solver->PrintStep             = NULL;
    solver->PhysicsOutput         = NULL;
    solver->PhysicsInput          = NULL;
    solver->AveragingFunction     = NULL;
    solver->GetLeftEigenvectors   = NULL;
    solver->GetRightEigenvectors  = NULL;
    solver->IBFunction            = NULL;

    if (!strcmp(solver->m_model,_LINEAR_ADVECTION_DIFFUSION_REACTION_)) {

      solver->m_physics = (LinearADR*) calloc (1,sizeof(LinearADR));
      IERR LinearADRInitialize(solver,mpi); CHECKERR(ierr);

    } else if (!strcmp(solver->m_model,_BURGERS_)) {

    solver->m_physics = (Burgers*) calloc (1,sizeof(Burgers));
    IERR BurgersInitialize(solver,mpi); CHECKERR(ierr);

    } else if (!strcmp(solver->m_model,_FP_DOUBLE_WELL_)) {

      solver->m_physics = (FPDoubleWell*) calloc (1,sizeof(FPDoubleWell));
      IERR FPDoubleWellInitialize(solver,mpi); CHECKERR(ierr);

    } else if (!strcmp(solver->m_model,_FP_POWER_SYSTEM_)) {

      solver->m_physics = (FPPowerSystem*) calloc (1,sizeof(FPPowerSystem));
      IERR FPPowerSystemInitialize(solver,mpi); CHECKERR(ierr);

    } else if (!strcmp(solver->m_model,_FP_POWER_SYSTEM_1BUS_)) {

      solver->m_physics = (FPPowerSystem1Bus*) calloc (1,sizeof(FPPowerSystem1Bus));
      IERR FPPowerSystem1BusInitialize(solver,mpi); CHECKERR(ierr);

    } else if (!strcmp(solver->m_model,_FP_POWER_SYSTEM_3BUS_)) {

      solver->m_physics = (FPPowerSystem3Bus*) calloc (1,sizeof(FPPowerSystem3Bus));
      IERR FPPowerSystem3BusInitialize(solver,mpi); CHECKERR(ierr);

    } else if (!strcmp(solver->m_model,_EULER_1D_)) {

      solver->m_physics = (Euler1D*) calloc (1,sizeof(Euler1D));
      IERR Euler1DInitialize(solver,mpi); CHECKERR(ierr);

    } else if (!strcmp(solver->m_model,_EULER_2D_)) {

      solver->m_physics = (Euler2D*) calloc (1,sizeof(Euler2D));
      IERR Euler2DInitialize(solver,mpi); CHECKERR(ierr);

    } else if (!strcmp(solver->m_model,_NAVIER_STOKES_2D_)) {

      solver->m_physics = (NavierStokes2D*) calloc (1,sizeof(NavierStokes2D));
      IERR NavierStokes2DInitialize(solver,mpi); CHECKERR(ierr);

    } else if (!strcmp(solver->m_model,_NAVIER_STOKES_3D_)) {

      solver->m_physics = (NavierStokes3D*) calloc (1,sizeof(NavierStokes3D));
      IERR NavierStokes3DInitialize(solver,mpi); CHECKERR(ierr);

    } else if (!strcmp(solver->m_model,_NUMA2D_)) {

      solver->m_physics = (Numa2D*) calloc (1,sizeof(Numa2D));
      IERR Numa2DInitialize(solver,mpi); CHECKERR(ierr);

    } else if (!strcmp(solver->m_model,_NUMA3D_)) {

      solver->m_physics = (Numa3D*) calloc (1,sizeof(Numa3D));
      IERR Numa3DInitialize(solver,mpi); CHECKERR(ierr);

    } else if (!strcmp(solver->m_model,_SHALLOW_WATER_1D_)) {

      solver->m_physics = (ShallowWater1D*) calloc (1,sizeof(ShallowWater1D));
      IERR ShallowWater1DInitialize(solver,mpi); CHECKERR(ierr);

    } else if (!strcmp(solver->m_model,_SHALLOW_WATER_2D_)) {

      solver->m_physics = (ShallowWater2D*) calloc (1,sizeof(ShallowWater2D));
      IERR ShallowWater2DInitialize(solver,mpi); CHECKERR(ierr);

    } else if (!strcmp(solver->m_model,_VLASOV_)) {

      solver->m_physics = (Vlasov*) calloc (1,sizeof(Vlasov));
      IERR VlasovInitialize(solver,mpi); CHECKERR(ierr);

    }else {

      fprintf(stderr,"Error (domain %d): %a_s is not a supported physical model.\n",
              ns, solver->m_model);
      return(1);

    }

    /* some checks */
    if ( ( (solver->GetLeftEigenvectors == NULL) || (solver->GetRightEigenvectors == NULL) )
        && (!strcmp(solver->m_interp_type,_CHARACTERISTIC_)) && (solver->m_nvars > 1) ) {
      if (!mpi->m_rank) {
        fprintf(stderr,"Error (domain %d): Interpolation type is defined as characteristic ", ns);
        fprintf(stderr,"but physics initializations returned NULL pointers for ");
        fprintf(stderr,"Get(Left,Right)Eigenvectors needed for characteristic-based ");
        fprintf(stderr,"reconstruction.\n");
      }
      return(1);
    }

    if (!strcmp(solver->m_split_hyperbolic_flux,"yes")) {
      if ((!solver->dFFunction) || (!solver->UpwinddF)) {
        if (!mpi->m_rank) {
          fprintf(stderr,"Error (domain %d): Splitting of hyperbolic flux requires a dFFunction ", ns);
          fprintf(stderr,"and its upwinding function UpwinddF.\n");
          fprintf(stderr,"Error: f(u) = [f(u) - df(u)] + df(u).\n");
          fprintf(stderr,"Error: dFFunction or UpwinddF (or both) is (are) NULL.\n");
        }
        return(1);
      }
      if (solver->FdFFunction && solver->UpwindFdF) solver->m_flag_fdf_specified = 1;
      else                                          solver->m_flag_fdf_specified = 0;
    }

    if ((solver->IBFunction == NULL) && (solver->m_flag_ib)) {
      if (!mpi->m_rank) {
        fprintf(stderr,"Error in InitializePhysics() (domain %d): Physical model %a_s does not yet have an immersed boundary treatment.\n",
                ns, solver->m_model);
      }
      return(1);
    }

  }

  return(0);
}

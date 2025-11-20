/*! @file petscinterface_struct.h
    @brief Contains structure that defines the interface for time integration with PETSc (https://petsc.org/release/)
    @author Debojyoti Ghosh
 */

#ifdef with_petsc

#ifndef _PETSC_INTERFACE_STRUCT_H_
#define _PETSC_INTERFACE_STRUCT_H_

#include <sys/time.h>
#include <vector>
#include <string>

/* include PETSc header files */
#include <petscsys.h>
#include <petscvec.h>
#include <petscmat.h>
#include <petscts.h>
#include <petscdmshell.h>

/* some definitions */

/*! Maximum size of character strings */
#define _MAX_STRING_SIZE_ 500
/*! Explicit time integration */
#define _EXPLICIT_  0
/*! Implicit time integration */
#define _IMPLICIT_  1

/*! \def PETScContext
    \brief Structure containing the variables for time-integration with PETSc
 * This structure contains all the variables needed to integration in time using
 * PETSc's TS module (https://petsc.org/release/docs/manualpages/TS/index.html).
*/

/*! \brief  Structure containing the variables for time-integration with PETSc
 *
 * This structure contains all the variables needed to integration in time using
 * PETSc's TS module (https://petsc.org/release/docs/manualpages/TS/index.html).
*/
typedef struct _petsccontext_ {

  /*! MPI rank */
  int m_rank;

  /*! number of MPI ranks */
  int m_nproc;

  /*! array of simulation objects of type #SimulationObject */
  void *m_simobj;

  /*! number of simulation objects */
  int m_nsims;

  /*! The shift variable in implicit time-integration */
  PetscReal m_shift;

  /*! Current time */
  double m_waqt;

  /*! Simulation time at start of step */
  double m_t_start;

  /*! Current dt */
  double m_dt;

  /*! Stage times */
  std::vector<double> m_stage_times;

  /*! Current stage */
  int m_stage_index;

  /*! A counter variable */
  int m_tic;

  /*! Number of computational points (local), i.e., total number of grid points
      not counting the ghost points */
  int m_npoints;

  /*! Number of computational DOFs (local), i.e., total number of grid points times
   *  the number of vector components at each grid point. */
  int m_ndofs;

  /*! Vector of index lists of local points in each simulation domain whose size is
     (#HyPar::m_ndims+1)*#PETScContext::m_npoints.
      For each point, it stores its ndims-dimensional index and its 1D index
      of its location in the array #HyPar::m_u. */
  std::vector<int*> m_points;

  /*! Vector offsets for each simulation domain */
  int* m_offsets;

  /* flags for implicit treatment */
  /*! Flag to indicate if hyperbolic term is treated implicitly or explicitly    */
  int m_flag_hyperbolic;
  /*! Flag to indicate if the split hyperbolic term (\a f - \a df) is treated implicitly or explicitly  */
  int m_flag_hyperbolic_f;
  /*! Flag to indicate if the split hyperbolic term \a df is treated implicitly or explicitly  */
  int m_flag_hyperbolic_df;
  /*! Flag to indicate if the parabolic term is treated implicitly or explicitly      */
  int m_flag_parabolic;
  /*! Flag to indicate if the source term is treated implicitly or explicitly         */
  int m_flag_source;

  /*! Flag to turn on preconditioning (off by default) */
  int m_flag_use_precon;

  /*! Construct or use provided PC matrix */
  std::string m_precon_matrix_type;

  /*! Flag to indicate if the system being solved for implicit time-integration is linear/nonlinear. */
  int m_flag_is_linear;

  /*! \f$\epsilon\f$ parameter for the Jacobian-free Newton-Krylov directional derivative computation */
  double m_jfnk_eps;

  /*! A vector of essentialy integer arrays of the same size and layout at the solutions
      (with ghost points) containing the global DOF index for each grid point in each
      simulation domain.
      It is declared as a \a double type so that its calculation can use functions
      defined for double-type variables.
      \sa PetscGlobalDOF() */
  std::vector<double*> m_globalDOF;

  /*! iteration start time */
  struct timeval m_iter_start_time;
  /*! iteration end time */
  struct timeval m_iter_end_time;
  /*! iteration wallclock time (in seconds) */
  double m_iter_wctime;
  /*! total time integration wall time (in seconds) */
  double m_ti_runtime;

#ifdef with_librom
  /*! libROM interface */
  void* m_rom_interface;
  /*! ROM mode \sa #libROMInterface::m_mode */
  std::string m_rom_mode;
  /*! Array of simulation times to write solution output at */
  std::vector<double> m_op_times_arr;
#endif
} PETScContext;

#endif

#endif


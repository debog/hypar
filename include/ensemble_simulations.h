/*! @file ensemble_simulations.h
    @brief Ensemble simulation class
    @author Debojyoti Ghosh
*/

#ifndef _ENSEMBLE_SIM_H_
#define _ENSEMBLE_SIM_H_

#include <stdio.h>
#include <vector>
#include <simulation.h>

#define _ENSEMBLE_SIM_INP_FNAME_ "simulation.inp"

/*! \class EnsembleSimulation
    \brief Class describing ensemble simulations
 *  This class contains all data and functions needed
 *  to run ensemble simulations, i.e., multiple simulations
 *  of the same physics with the same numerical methods but
 *  on multiple grids.
*/

/*! \brief Class describing ensemble simulations
 *
 * This class contains all data and functions needed
 * to run ensemble simulations, i.e., multiple simulations
 * of the same physics with the same numerical methods but
 * on multiple grids.
*/
class EnsembleSimulation : public Simulation
{
  public:

    /*! Constructor */
    EnsembleSimulation()
    {
      m_is_defined = false;
      m_sims.clear();
      m_nsims = 0;
      m_nproc = -1;
      m_rank = -1;
    }

    /*! Destructor */
    ~EnsembleSimulation()
    {
      int err = Cleanup((void*) m_sims.data(), m_nsims);
      if (err) {
        printf( "Error: CleanUp() returned with status %d on process %d.\n",
                err, m_rank );
      }
      m_sims.clear();
    }

    /*! Define this object */
    int define(int, int);

    /*! Read solver inputs */
    inline int ReadInputs()
    {
      int retval = ::ReadInputs(  (void*) m_sims.data(),
                                  m_nsims,
                                  m_rank );
      ::WriteInputs( (void*) m_sims.data(),
                     m_nsims,
                     m_rank );
      return retval;
    }

    /*! Initialize the simulations */
    inline int Initialize()
    {
      int retval = ::Initialize(  (void*) m_sims.data(),
                                  m_nsims );
      return retval;
    }

    /*! Read initial solution for each simulation */
    inline int InitialSolution()
    {
      int retval = ::InitialSolution( (void*) m_sims.data(),
                                      m_nsims );
      return retval;
    }

    /*! Initialize the boundary conditions of the simulations */
    inline int InitializeBoundaries()
    {
      int retval = ::InitializeBoundaries(  (void*) m_sims.data(),
                                            m_nsims );
      return retval;
    }

    /*! Initialize the immersed boundary conditions of the simulations */
    inline int InitializeImmersedBoundaries()
    {
      int retval = ::InitializeImmersedBoundaries(  (void*) m_sims.data(),
                                                    m_nsims );
      return retval;
    }

    /*! Initialize the physics of the simulations */
    inline int InitializePhysics()
    {
      int retval = ::InitializePhysics( (void*) m_sims.data(),
                                        m_nsims );
      return retval;
    }

    /*! Initialize the physics data of the simulations */
    inline int InitializePhysicsData()
    {
      for (int ns = 0; ns < m_nsims; ns++) {
        int retval = ::InitializePhysicsData( (void*) &(m_sims[ns]),
                                              ns, m_nsims, NULL );
        if (retval) {
          fprintf(stderr, "Error in EnsembleSimulations::InitializePhysicsData()\n");
          fprintf(stderr, "  InitializePhysicsData returned with error code %d on rank %d.\n",
                  retval, m_sims[ns].mpi.rank);
          return retval;
        }
      }
      return 0;
    }

    /*! Initialize the numerical solvers of the simulations */
    inline int InitializeSolvers()
    {
      int retval = ::InitializeSolvers( (void*) m_sims.data(),
                                        m_nsims );
      return retval;
    }

    /*! Run the simulation using native time integrators */
    inline int Solve()
    {
      int retval = ::Solve( (void*) m_sims.data(),
                            m_nsims,
                            m_rank,
                            m_nproc );
      return retval;
    }

    /*! Write simulation errors and wall times to file */
    inline void WriteErrors(double a_wt_solver,
                            double a_wt_total )
    {
      ::SimWriteErrors( (void*) m_sims.data(),
                        m_nsims,
                        m_rank,
                        a_wt_solver,
                        a_wt_total );
      return;
    }

    /*! Return whether object is defined or not */
    inline bool isDefined() const { return m_is_defined; }

#ifndef serial
    /*! Create duplicate MPI communicators */
    inline int mpiCommDup()
    {
      for (int n = 0; n < m_nsims; n++) {
        MPI_Comm_dup(MPI_COMM_WORLD, &(m_sims[n].mpi.world));
      }
      return 0;
    }
#endif

#ifdef with_petsc
    /*! Set flag whether to use PETSc time integration */
    inline void usePetscTS(PetscBool a_flag)
    {
      for (int n = 0; n < m_nsims; n++) {
        m_sims[n].solver.use_petscTS  = a_flag;
      }
    }

    /*! Run the simulation using PETSc time integrators */
    inline int SolvePETSc()
    {
      int retval = ::SolvePETSc(  (void*) m_sims.data(),
                                  m_nsims,
                                  m_rank,
                                  m_nproc );
      return retval;
    }
#endif

  protected:

    bool  m_is_defined;     /*!< Boolean to show if this object is defined */
    int   m_nsims;          /*!< Number of ensemble simulations */
    int   m_rank,           /*!< MPI rank of this process */
          m_nproc;          /*!< Total number of MPI ranks */

    std::vector<SimulationObject> m_sims; /*!< vector of simulation objects */

  private:

};

#endif

/*! @file single_simulation.h
    @brief Single simulation class
    @author Debojyoti Ghosh
*/

#ifndef _SINGLE_SIM_H_
#define _SINGLE_SIM_H_

#include <stdio.h>
#include <vector>
#include <simulation.h>

#define _NSIMS_ 1

/*! \class SingleSimulation
    \brief Class describing a single simulation
 *  This class contains all data and functions needed
 *  to run a single simulation.
*/

/*! \brief Class describing a single simulation
 *
 * This class contains all data and functions needed
 * to run a single simulation.
*/
class SingleSimulation : public Simulation
{
  public:

    /*! Constructor */
    SingleSimulation()
    {
      m_is_defined = false;
      m_sim = nullptr;
      m_nproc = -1;
      m_rank = -1;
    }

    /*! Destructor */
    ~SingleSimulation()
    {
      int err = Cleanup((void*) m_sim, _NSIMS_);
      if (err) {
        printf( "Error: CleanUp() returned with status %d on process %d.\n",
                err, m_rank );
      }
      delete m_sim;
    }

    /*! Define this object */
    int define(int, int);

    /*! Read solver inputs */
    inline int ReadInputs()
    {
      int retval = ::ReadInputs(  (void*) m_sim,
                                  _NSIMS_,
                                  m_rank );
      ::WriteInputs( (void*) m_sim,
                     _NSIMS_,
                     m_rank );
      return retval;
    }

    /*! Initialize the simulations */
    inline int Initialize()
    {
      int retval = ::Initialize(  (void*) m_sim,
                                  _NSIMS_ );
      return retval;
    }

    /*! Read initial solution for each simulation */
    inline int InitialSolution()
    {
      int retval = ::InitialSolution( (void*) m_sim,
                                      _NSIMS_ );
      return retval;
    }

    /*! Initialize the boundary conditions of the simulations */
    inline int InitializeBoundaries()
    {
      int retval = ::InitializeBoundaries(  (void*) m_sim,
                                            _NSIMS_ );
      return retval;
    }

    /*! Initialize the immersed boundary conditions of the simulations */
    inline int InitializeImmersedBoundaries()
    {
      int retval = ::InitializeImmersedBoundaries(  (void*) m_sim,
                                                    _NSIMS_ );
      return retval;
    }

    /*! Initialize the physics of the simulations */
    inline int InitializePhysics()
    {
      int retval = ::InitializePhysics( (void*) m_sim,
                                        _NSIMS_ );
      return retval;
    }

    /*! Initialize the physics data of the simulations */
    inline int InitializePhysicsData()
    {
      int retval = ::InitializePhysicsData( (void*) m_sim,
                                            0, _NSIMS_, NULL );
      return retval;
    }

    /*! Initialize the numerical solvers of the simulations */
    inline int InitializeSolvers()
    {
      int retval = ::InitializeSolvers( (void*) m_sim,
                                        _NSIMS_ );
      return retval;
    }

    /*! Run the simulation using native time integrators */
    inline int Solve()
    {
      int retval = ::Solve( (void*) m_sim,
                            _NSIMS_,
                            m_rank,
                            m_nproc );
      return retval;
    }

    /*! Write simulation errors and wall times to file */
    inline void WriteErrors(double a_wt_solver,
                            double a_wt_total )
    {
      ::SimWriteErrors( (void*) m_sim,
                        _NSIMS_,
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
      MPI_Comm_dup(MPI_COMM_WORLD, &(m_sim->mpi.world));
      return 0;
    }
#endif

#ifdef with_petsc
    /*! Set flag whether to use PETSc time integration */
    inline void usePetscTS(PetscBool a_flag)
    {
      m_sim->solver.use_petscTS  = a_flag;
    }

    /*! Run the simulation using PETSc time integrators */
    inline int SolvePETSc()
    {
      int retval = ::SolvePETSc(  (void*) m_sim,
                                  _NSIMS_,
                                  m_rank,
                                  m_nproc );
      return retval;
    }
#endif

  protected:

    bool  m_is_defined;       /*!< Boolean to show if this object is defined */
    int   m_rank,             /*!< MPI rank of this process */
          m_nproc;            /*!< Total number of MPI ranks */

    SimulationObject* m_sim;  /*< Simulation object */

  private:

};

#endif

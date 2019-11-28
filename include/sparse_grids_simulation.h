/*! @file sparse_grids_simulations.h
    @brief Sparse grids simulation class
    @author Debojyoti Ghosh, John Loffeld, Lee Ricketson
*/

#ifndef _SPARSE_GRIDS_SIM_H_
#define _SPARSE_GRIDS_SIM_H_

#include <vector>
#include <utility>
#include <simulation.h>

#define _SPARSEGRIDS_SIM_INP_FNAME_ "sparse_grids.inp"

typedef std::vector<int> GridDimensions;
typedef std::vector<int> ProcDistribution;
typedef std::pair<double, GridDimensions> SGCTElem;

#define _coeff_ first
#define _dim_ second

/*! \class SparseGridsSimulation
    \brief Class describing sparse grids simulations
 *  This class contains all data and functions needed
 *  to run a sparse grids simulation.
*/

/*! \brief Class describing sparse grids simulations
 *
 * This class contains all data and functions needed
 * to run sparse grids simulations.
*/
class SparseGridsSimulation : public Simulation
{
  public:

    /*! Constructor */
    SparseGridsSimulation()
    {
      m_is_defined = false;

      m_nsims_sg = -1;
      m_ndims = -1;
      m_rank = -1;
      m_nproc = -1;

      m_sim_fg = NULL;
      m_sims_sg.clear();

      m_n_fg = -1;
      m_imin = -1;

      m_combination.clear();
    }

    /*! Destructor */
    ~SparseGridsSimulation()
    {
      int err;
      /* clean up full grid object */
      err = CleanupBarebones(m_sim_fg);
      delete m_sim_fg;
      /* clean up and delete sparse grids objects */
      err = Cleanup((void*) m_sims_sg.data(), m_nsims_sg);
      if (err) {
        printf( "Error: CleanUp() returned with status %d on process %d.\n",
                err, m_rank );
      }
      m_sims_sg.clear();
      /* done */
      return;
    }

    /*! Define this object */
    int define(int, int);

    /*! Read solver inputs for the full grid simulation object */
    inline int ReadInputs()
    {
      int retval = ::ReadInputs(  (void*) m_sim_fg,
                                  1, 
                                  m_rank );
      return retval;
    }

    /*! Initialize the simulations */
    int Initialize();

    /*! Read initial solution for each simulation */
    inline int InitialSolution()
    {
      int retval = ::InitialSolution( (void*) m_sim_fg, 1 );
      return retval;
    }

    /*! Initialize the boundary conditions of the simulations */
    inline int InitializeBoundaries()
    {
      int retval = ::InitializeBoundaries(  (void*) m_sims_sg.data(),
                                            m_nsims_sg );
      return retval;
    }

    /*! Initialize the immersed boundary conditions of the simulations */
    inline int InitializeImmersedBoundaries()
    {
      int retval = ::InitializeImmersedBoundaries(  (void*) m_sims_sg.data(),
                                                    m_nsims_sg );
      return retval;
    }

    /*! Initialize the physics of the simulations */
    inline int InitializePhysics()
    {
      int retval = ::InitializePhysics( (void*) m_sims_sg.data(),
                                        m_nsims_sg );
      return retval;
    }

    /*! Initialize the numerical solvers of the simulations */
    inline int InitializeSolvers()
    {
      int retval = ::InitializeSolvers( (void*) m_sims_sg.data(),
                                        m_nsims_sg );

      /* disable output from individual sparse grids */
      for (int i = 0; i < m_nsims_sg; i++) {
        m_sims_sg[i].solver.WriteOutput = NULL;
      }

      return retval;
    }

    /*! Run the simulation using native time integrators */
    int Solve();

    /*! Write simulation errors and wall times to file */
    inline void WriteErrors(double a_wt_solver,
                            double a_wt_total )
    {
      ::SimWriteErrors( (void*) m_sim_fg,
                        1,
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
      if (!m_sim_fg) {
        fprintf(stderr, "Error in SparseGridsSimulation::mpiCommDup()\n");
        fprintf(stderr, "  m_sim_fg is not allocated on rank %d!\n", m_rank);
        return 1;
      }
      MPI_Comm_dup(MPI_COMM_WORLD, &(m_sim_fg->mpi.world));
      return 0;
    }
#endif

#ifdef with_petsc
    /*! Set flag whether to use PETSc time integration */
    inline void usePetscTS(PetscBool a_flag)
    {
      m_sim_fg->solver.use_petscTS = a_flag;
    }

    /*! Run the simulation using PETSc time integrators */
    inline int SolvePETSc()
    {
      int retval = ::SolvePETSc(  (void*) m_sims_sg.data(),
                                  m_nsims_sg,
                                  m_rank,
                                  m_nproc );
      return retval;
    }
#endif

  protected:

    bool  m_is_defined;   /*!< Boolean to show if this object is defined */

    int   m_nsims_sg;     /*!< Number of sparse grids simulations */
    int   m_ndims;        /*!< Number of spatial dimensions */
    int   m_rank,         /*!< MPI rank of this process */
          m_nproc;        /*!< Total number of MPI ranks */

    SimulationObject* m_sim_fg;               /*!< full grid simulation object */
    std::vector<SimulationObject> m_sims_sg;  /*!< vector of sparse grids simulation objects */

    int m_n_fg; /*!< Base2 log of the number of grid points along a dimension of the full grid */
    int m_imin; /*!< Base2 log of the minimum number of grid points along any dimension */

    std::vector<SGCTElem> m_combination;  /*!< Coefficients and grid dimensions for the combination technique */
    std::vector<ProcDistribution> m_iprocs; /*!< MPI ranks along each dimension for the grids in the combination technique */

    /*! Some basic sanity checks */
    int SanityChecks();

    /*! Compute the number, coefficients, and dimensions of all the sparse
     *  grids that go into the combination technique. */
    int ComputeSGDimsAndCoeffs();

    /*! Compute the load-balanced MPI ranks distributions */
    int ComputeProcessorDistribution( ProcDistribution&,
                                      const GridDimensions&);

    /*! Compute all grids such that the base2 logs of the size in each dimension
     *  adds up to the input argument */
    void GetCTGridSizes(const int, std::vector<GridDimensions>&);

    /*! Initialize a "barebones" simulation object */
    int InitializeBarebones(SimulationObject*);

    /*! Cleanup (deallocate) a "barebones" simulation object */
    int CleanupBarebones(SimulationObject*);

    /*! Set solver parameters for a simulation object */
    int SetSolverParameters(  SimulationObject&,
                              const GridDimensions&,
                              const ProcDistribution&,
                              const SimulationObject&,
                              const int,
                              const int );

    /*! Checks if an integer is a power of 2 */
    inline bool isPowerOfTwo(int x) 
    {
      if (x == 0)  return false;

      while (x > 1) {
        if (x%2 != 0) return false;
        x /= 2;
      }
      return true;
    }

    /*! The combination/choose function from probability stuff */
    inline double cfunc(int a, int b)
    {
      return ( ( (double)factorial(a) ) / ( (double) (factorial(b)*factorial(a-b)) ) );
    }

    /*! Factorial function */
    inline int factorial(int a)
    {
      int retval = 1.0;
      for (int i=1; i<=a; i++) retval *= a;
      return retval;
    }

  private:

};

#endif

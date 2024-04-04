/*! @file sparse_grids_simulations.h
    @brief Sparse grids simulation class
    @author Debojyoti Ghosh, John Loffeld, Lee Ricketson
*/

#ifndef _SPARSE_GRIDS_SIM_H_
#define _SPARSE_GRIDS_SIM_H_

#include <stdio.h>
#include <cstdlib>
#include <vector>
#include <string.h>
#include <string>
#include <utility>
#include <std_vec_ops.h>
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

      m_interp_order = 6;

      m_is_periodic.clear();
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

    /*! Initialize the simulation */
    int Initialize();

    /*! Read initial solution */
    int InitialSolution();

    /*! Initialize the boundary conditions of the simulation */
    inline int InitializeBoundaries()
    {
      int retval = ::InitializeBoundaries(  (void*) m_sims_sg.data(),
                                            m_nsims_sg );
      if (m_nsims_sg > 0) {
        m_is_periodic.resize(m_ndims);
        for (int d=0; d<m_ndims; d++) {
          m_is_periodic[d] = ( m_sims_sg[0].solver.isPeriodic[d] == 1 ? true : false);
          m_sim_fg->solver.isPeriodic[d] = m_sims_sg[0].solver.isPeriodic[d];
        }
      }
      return retval;
    }

    /*! Initialize the immersed boundary conditions of the simulation */
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

    /*! Initialize the physics of the simulations */
    inline int InitializePhysicsData()
    {
      for (int ns = 0; ns < m_sims_sg.size(); ns++) {
        int retval = ::InitializePhysicsData( (void*) &(m_sims_sg[ns]),
                                              -1,
                                              m_nsims_sg,
                                              m_sim_fg->solver.dim_global);
        if (retval) {
          fprintf(stderr, "Error in SparseGridsSimulation::InitializePhysicsData()\n");
          fprintf(stderr, "  InitializePhysicsData returned with error code %d on rank %d.\n",
                  retval, m_sims_sg[ns].mpi.rank);
          return retval;
        }
      }
      return 0;
    }

    /*! Initialize the numerical solvers of the simulations */
    inline int InitializeSolvers()
    {
      int retval = ::InitializeSolvers( (void*) m_sims_sg.data(),
                                        m_nsims_sg );
      InitializeSolversBarebones(m_sim_fg);

      /* some modifications to output filename roots */
      for (int i=0; i<m_nsims_sg; i++) {
        strcpy(m_sims_sg[i].solver.op_fname_root, "op_sg");
        strcpy(m_sims_sg[i].solver.aux_op_fname_root, "ts0_sg");
      }
      strcpy(m_sim_fg->solver.op_fname_root, "op_fg");
      strcpy(m_sim_fg->solver.aux_op_fname_root, "ts0_fg");
      return retval;
    }

    /*! Wrap up initializations */
    int InitializationWrapup();

    /*! Run the simulation using native time integrators */
    int Solve();

    /*! Write simulation errors and wall times to file */
    void WriteErrors(double, double);

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

    /*! Order of interpolation between grids (input - \b sparse_grids.inp ) */
    int   m_interp_order;

    std::vector<bool> m_is_periodic; /*!< Periodicity along each dimension */

    /*! Write out the sparse grid solutions to file? (input - \b sparse_grids.inp ) */
    int m_write_sg_solutions;

    /*! Print and write out the sparse grid errors? (input - \b sparse_grids.inp ) */
    int m_print_sg_errors;

    SimulationObject* m_sim_fg;               /*!< full grid simulation object */
    std::vector<SimulationObject> m_sims_sg;  /*!< vector of sparse grids simulation objects */

    int m_n_fg; /*!< Base2 log of the number of grid points along a dimension of the full grid */
    int m_imin; /*!< Base2 log of the minimum number of grid points along any dimension (input - \b sparse_grids.inp ) */

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

    /*! Initialize some function pointers for a "barebones" simulation object */
    int InitializeSolversBarebones(SimulationObject*);

    /*! Cleanup (deallocate) a "barebones" simulation object */
    int CleanupBarebones(SimulationObject*);

    /*! Set solver parameters for a simulation object */
    int SetSolverParameters(  SimulationObject&,
                              const GridDimensions&,
                              const ProcDistribution&,
                              const SimulationObject&,
                              const int,
                              const int );

    /*! Output solutions to file */
    void OutputSolution(double);

    /*! Combination technique */
    void CombinationTechnique(SimulationObject* const);

    /*! Interpolate data from one simulation object to another */
    void interpolate( SimulationObject* const,
                      const SimulationObject* const);

    /*! Interpolate data from a simulation to a global C-array */
    void interpolate( const GridDimensions&,
                      double** const,
                      const SimulationObject* const);

    /*! Interpolate data from one simulation object to another */
    void interpolateGrid( SimulationObject* const,
                          const SimulationObject* const);

    /*! Coarsen along a given dimension */
    void coarsenGrid1D( const GridDimensions&,
                        const GridDimensions&,
                        const double* const,
                        double* const,
                        int );

    /*! Refine along a given dimension */
    void refineGrid1D( const GridDimensions&,
                       const GridDimensions&,
                       const double* const,
                       double* const,
                       int );

    /*! Fill ghost cells for interpolation */
    void fillGhostCells(  const GridDimensions&,
                          const int,
                          double* const,
                          const int );

    /*! Calculate errors */
    void CalculateError();

    /*! Compute error for a simulation object */
    void computeError( SimulationObject&, double* const);

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
      for (int i=1; i<=a; i++) retval *= i;
      return retval;
    }

    /*! Allocate grid array given grid dimension */
    inline void allocateGridArrays( const GridDimensions& a_dim,
                                    double** const        a_x,
                                    const int             a_ngpt = 0)
    {
      GridDimensions dim_wghosts = a_dim;
      for (int i=0; i<dim_wghosts.size(); i++) {
        dim_wghosts[i] += (2*a_ngpt);
      }
      long size_x = StdVecOps::sum(dim_wghosts);
      (*a_x) = (double*) calloc (size_x, sizeof(double));
      for (int i=0; i<size_x; i++) (*a_x)[i] = 0.0;
      return;
    }

    /*! Allocate data arrays given grid dimension */
    inline void allocateDataArrays( const GridDimensions& a_dim,
                                    const int             a_nvars,
                                    double** const        a_u,
                                    const int             a_ngpt = 0)
    {
      GridDimensions dim_wghosts = a_dim;
      for (int i=0; i<dim_wghosts.size(); i++) {
        dim_wghosts[i] += (2*a_ngpt);
      }
      long size_u = a_nvars*StdVecOps::product(dim_wghosts);
      (*a_u) = (double*) calloc (size_u, sizeof(double));
      for (int i=0; i<size_u; i++) (*a_u)[i] = 0.0;
      return;
    }

  private:

};

#endif

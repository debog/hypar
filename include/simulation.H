#include <simulation.h>

extern "C" int ReadInputs(void*,int,int);/*!< Read the input parameters */
extern "C" int Initialize(void*,int);/*!< Initialize the solver */
extern "C" int InitialSolution(void*,int);/*!< Read the initial solution */
extern "C" int InitializeBoundaries(void*,int);/*!< Initialize the boundary conditions */
extern "C" int InitializeImmersedBoundaries(void*,int);/*!< Initialize the immersed boundary conditions */
extern "C" int InitializePhysics(void*,int);/*!< Initialize the physics */
extern "C" int InitializeSolvers(void*,int);/*!< Initialize the solvers */
extern "C" int Solve(void*,int, int, int);/*!< Solve the PDE - time-integration */
#ifdef with_petsc
extern "C" int SolvePETSc(void*,int, int, int);  /*!< Solve the PDE using PETSc TS */
#endif
extern "C" int Cleanup(void*,int);/*!< Clean up: deallocate all arrays and objects */

/*! Write errors for each simulation */
extern "C" void SimWriteErrors(void*, int, int, double, double);



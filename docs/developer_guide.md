# Developer Guide

This guide provides information for developers who want to extend HyPar or add new physical models.

## Code Structure

HyPar is organized into several main directories:

```
hypar/
├── include/               # Header files
│   ├── basic.h           # Basic definitions
│   ├── hypar.h           # Main solver structure
│   ├── mpivars.h         # MPI variables
│   ├── boundaryconditions.h  # Boundary condition types
│   ├── physicalmodels/   # Physical model headers
│   ├── arrayfunctions.h  # Array manipulation macros
│   └── ...
├── src/
│   ├── HyParFunctions/   # Core solver functions
│   ├── PhysicalModels/   # Physical model implementations
│   ├── BoundaryConditions/  # Boundary condition implementations
│   ├── TimeIntegration/  # Time stepping schemes
│   ├── InterpolationFunctions/  # Spatial discretization
│   ├── MathFunctions/    # Mathematical utilities
│   └── Simulation/       # Main simulation drivers
├── Examples/             # Example test cases
└── docs/                 # Documentation
```

## Adding a New Physical Model

### Overview

Physical models in HyPar define the governing equations, including:
- **Hyperbolic flux** functions (advection/convection)
- **Parabolic terms** (diffusion) - optional
- **Source terms** - optional
- **Upwinding schemes** for flux computation
- **Eigenstructure** for characteristic-based methods - optional

### Prerequisites

Before implementing a new model:

1. **Understand your PDE**: Know the mathematical form of all terms
2. **Choose a reference model**: Use an existing similar model as a template:
   - Simple scalar: `Burgers` (`src/PhysicalModels/Burgers/`)
   - Linear ADR: `LinearADR` (`src/PhysicalModels/LinearADR/`)
   - Hyperbolic system: `Euler1D` (`src/PhysicalModels/Euler1D/`)
   - Hyperbolic-parabolic: `NavierStokes2D` (`src/PhysicalModels/NavierStokes2D/`)

3. **Understand HyPar conventions**:
   - Solution stored as 1D array: `u[nvars*p+v]` where `p` = grid point, `v` = variable
   - Use array index macros: `_ArrayIndex1D_`, `_ArrayIncrementIndex_`, etc.
   - Include ghost points in calculations

### Implementation Steps

#### 1. Create Header File

Create `include/physicalmodels/yourmodel.h`:

```c
#ifndef _YOURMODEL_H_
#define _YOURMODEL_H_

#define _YOUR_MODEL_ "yourmodel"

#include <basic.h>

/* Model-specific parameters */
typedef struct yourmodel_parameters {
  double param1;     // Physical parameter
  double param2;     // Another parameter
  // Add more parameters as needed
} YourModel;

/* Function declarations */
int YourModelInitialize  (void*,void*);
int YourModelCleanup     (void*);

#endif
```

#### 2. Create Source Directory

```bash
mkdir src/PhysicalModels/YourModel
```

#### 3. Implement Required Functions

Every model must implement these functions:

##### Initialize (`YourModelInitialize.c`)

```c
#include <physicalmodels/yourmodel.h>
#include <hypar.h>
#include <mpivars.h>

/* Forward declarations */
double YourModelComputeCFL(void*,void*,double,double);
int    YourModelAdvection(double*,double*,int,void*,double);
int    YourModelUpwind(double*,double*,double*,double*,
                       double*,double*,int,void*,double);

int YourModelInitialize(void *s, void *m)
{
  HyPar        *solver  = (HyPar*)        s;
  MPIVariables *mpi     = (MPIVariables*) m;
  YourModel    *physics = (YourModel*)    solver->physics;

  /* Set default parameters */
  physics->param1 = 1.0;
  physics->param2 = 0.0;

  /* Read from physics.inp (on rank 0) */
  if (!mpi->rank) {
    FILE *in = fopen("physics.inp","r");
    if (in) {
      char word[100];
      fscanf(in,"%s",word);
      if (!strcmp(word, "begin")) {
        while (strcmp(word, "end")) {
          fscanf(in,"%s",word);
          if (!strcmp(word, "param1")) {
            fscanf(in,"%lf",&physics->param1);
          } else if (!strcmp(word, "param2")) {
            fscanf(in,"%lf",&physics->param2);
          }
        }
      }
      fclose(in);
    }
  }

  /* Broadcast parameters to all MPI ranks */
#ifndef serial
  MPIBroadcast_double(&physics->param1,1,0,&mpi->world);
  MPIBroadcast_double(&physics->param2,1,0,&mpi->world);
#endif

  /* Print parameters */
  if (!mpi->rank) {
    printf("YourModel: param1 = %lf\n", physics->param1);
    printf("YourModel: param2 = %lf\n", physics->param2);
  }

  /* Register function pointers */
  solver->ComputeCFL = YourModelComputeCFL;
  solver->FFunction  = YourModelAdvection;
  solver->Upwind     = YourModelUpwind;

  return(0);
}
```

##### Compute CFL (`YourModelComputeCFL.c`)

```c
#include <physicalmodels/yourmodel.h>
#include <hypar.h>

double YourModelComputeCFL(void *s, void *m, double dt, double t)
{
  HyPar     *solver  = (HyPar*)     s;
  YourModel *physics = (YourModel*) solver->physics;

  int ndims  = solver->ndims;
  int *dim   = solver->dim_local;
  int ghosts = solver->ghosts;

  double max_cfl = 0;
  int done = 0;
  int index[ndims];

  _ArraySetValue_(index,ndims,0);
  while (!done) {
    int p;
    _ArrayIndex1D_(ndims,dim,index,ghosts,p);

    for (int dir=0; dir<ndims; dir++) {
      double dxinv;
      _GetCoordinate_(dir,index[dir],dim,ghosts,solver->dxinv,dxinv);

      /* Compute maximum wave speed */
      double wave_speed = /* model-dependent */;
      double local_cfl = fabs(wave_speed) * dt * dxinv;

      if (local_cfl > max_cfl) max_cfl = local_cfl;
    }
    _ArrayIncrementIndex_(ndims,dim,index,done);
  }

  return(max_cfl);
}
```

##### Hyperbolic Flux (`YourModelAdvection.c`)

```c
#include <physicalmodels/yourmodel.h>
#include <hypar.h>

int YourModelAdvection(double *f, double *u, int dir, 
                       void *s, double t)
{
  HyPar     *solver  = (HyPar*)     s;
  YourModel *physics = (YourModel*) solver->physics;

  int ndims  = solver->ndims;
  int nvars  = solver->nvars;
  int *dim   = solver->dim_local;
  int ghosts = solver->ghosts;

  /* Loop over all grid points and compute flux */
  int done = 0;
  int index[ndims];
  _ArraySetValue_(index,ndims,0);

  while (!done) {
    int p;
    _ArrayIndex1D_(ndims,dim,index,ghosts,p);

    for (int v=0; v<nvars; v++) {
      /* Compute flux: f[nvars*p+v] = ... */
      f[nvars*p+v] = /* model-specific flux */;
    }

    _ArrayIncrementIndex_(ndims,dim,index,done);
  }

  return(0);
}
```

##### Upwinding (`YourModelUpwind.c`)

```c
#include <physicalmodels/yourmodel.h>
#include <hypar.h>

int YourModelUpwind(double *fI, double *fL, double *fR,
                    double *uL, double *uR, double *u,
                    int dir, void *s, double t)
{
  HyPar     *solver  = (HyPar*)     s;
  YourModel *physics = (YourModel*) solver->physics;

  int ndims  = solver->ndims;
  int nvars  = solver->nvars;
  int *dim   = solver->dim_local;

  /* Loop over interfaces and compute upwind flux */
  /* Implementation depends on upwinding method:
     - Simple upwind (check wave direction)
     - Roe scheme
     - Local Lax-Friedrichs
     - HLLC, etc.
  */

  return(0);
}
```

##### Cleanup (`YourModelCleanup.c`)

```c
#include <physicalmodels/yourmodel.h>

int YourModelCleanup(void *s)
{
  YourModel *physics = (YourModel*) s;

  /* Free any allocated arrays */
  /* if (physics->some_array) free(physics->some_array); */

  return(0);
}
```

#### 4. Register Model with HyPar

Edit `src/Simulation/InitializePhysics.c`:

Add include:
```c
#include <physicalmodels/yourmodel.h>
```

Add registration block:
```c
} else if (!strcmp(solver->model,_YOUR_MODEL_)) {
  solver->physics = (YourModel*) calloc (1,sizeof(YourModel));
  YourModelInitialize(solver,mpi);
```

#### 5. Update Build System

Create `src/PhysicalModels/YourModel/Makefile.am`:

```makefile
noinst_LIBRARIES = libyourmodel.a
libyourmodel_a_SOURCES = \
  YourModelInitialize.c \
  YourModelCleanup.c \
  YourModelAdvection.c \
  YourModelUpwind.c \
  YourModelComputeCFL.c
```

Add to `src/PhysicalModels/Makefile.am`:
```makefile
SUBDIRS = ... YourModel
```

#### 6. Rebuild

```bash
autoreconf -i
./configure [options]
make
```

### Optional Functions

Implement as needed:

- **Parabolic terms**: `YourModelDiffusion` → `solver->GFunction`
- **Source terms**: `YourModelSource` → `solver->SFunction`
- **Eigenstructure**: `YourModelLeftEigenvectors`, `YourModelRightEigenvectors`
- **Jacobian** (for implicit methods): `YourModelJacobian` → `solver->JFunction`

### Testing Your Model

1. Create test case in `Examples/YourModel/`
2. Create input files: `solver.inp`, `physics.inp`, `boundary.inp`, `initial.inp`
3. Run: `../../bin/HyPar`
4. Verify: Check CFL, conservation, stability
5. Validate: Compare with exact/reference solutions

## Code Conventions

### Array Indexing

Solution array is 1D: `u[nvars*grid_point + variable]`

Use provided macros:
```c
_ArrayIndex1D_(ndims, dim, index, ghosts, p);  // Get 1D index p
_ArrayIncrementIndex_(ndims, dim, index, done); // Increment multi-D index
_ArraySetValue_(array, size, value);            // Set all elements
_ArrayCopy1D_(source, dest, size);              // Copy array
```

### Ghost Points

Grid layout includes ghost points:
- `dim_local[]` - local grid size (without ghosts)
- `ghosts` - number of ghost points on each side
- Physical domain: points with `ghosts` ≤ index < `dim_local[d] + ghosts`

### MPI Communication

For parallel runs:
```c
#ifndef serial
  // MPI-specific code
  MPIBroadcast_double(data, count, root, comm);
  // etc.
#endif
```

Always broadcast parameters read from files!

## Debugging Tips

1. **Enable verbose output**: Set `screen_op_iter 1` in `solver.inp`
2. **Check array bounds**: Use small grids initially
3. **Serial first**: Test with `--enable-serial` before MPI
4. **Valgrind**: Check for memory leaks: `valgrind --leak-check=full ./bin/HyPar`
5. **GDB**: Debug crashes: `gdb ./bin/HyPar`, then `run`

## Documentation

Use Doxygen style comments:

```c
/*! @file YourModelAdvection.c
    @brief Computes the hyperbolic flux
    @author Your Name
*/

/*!
  Calculate the hyperbolic flux for the your model equations.
  
  @param f     Computed flux (output)
  @param u     Solution array (input)
  @param dir   Spatial direction
  @param s     Solver object
  @param t     Current time
  @return      Error code (0 = success)
*/
int YourModelAdvection(double *f, double *u, int dir, void *s, double t);
```

## Contributing

1. Fork the repository on GitHub
2. Create a feature branch: `git checkout -b feature/your-model`
3. Commit your changes with clear messages
4. Test thoroughly (serial and parallel)
5. Submit a pull request

## Additional Resources

- **Full developer guide**: See `doc/Adding_Physical_Models.md` in the source
- **Existing models**: Browse `src/PhysicalModels/` for examples
- **Doxygen documentation**: http://hypar.github.io/
- **HyPar paper**: Ghosh & Constantinescu (2018), https://doi.org/10.1016/j.jcp.2018.02.009

## Getting Help

- **GitHub Issues**: https://github.com/debog/hypar/issues
- **Email**: debojyoti.ghosh@gmail.com
- **Check examples**: Most questions answered by existing model implementations

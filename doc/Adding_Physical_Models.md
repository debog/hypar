# Developer Guide: Adding New Physical Models to HyPar

This guide provides step-by-step instructions for implementing a new physical model in HyPar.

## Table of Contents
1. [Overview](#overview)
2. [Prerequisites](#prerequisites)
3. [Step-by-Step Implementation](#step-by-step-implementation)
4. [Required Functions](#required-functions)
5. [Optional Functions](#optional-functions)
6. [Testing and Validation](#testing-and-validation)
7. [Example: Implementing a Simple Model](#example-implementing-a-simple-model)

---

## Overview

Physical models in HyPar define the mathematical equations to be solved, including:
- Hyperbolic flux functions
- Parabolic (diffusion) terms
- Source terms
- Upwinding schemes
- Eigenstructure (for characteristic-based methods)

Each physical model is implemented as a self-contained module with its own:
- Header file (in `include/physicalmodels/`)
- Source directory (in `src/PhysicalModels/`)
- Physics-specific data structure
- Function pointers registered with the solver

---

## Prerequisites

Before implementing a new physical model, you should:

1. **Understand your PDE**: Know the mathematical form of:
   - Hyperbolic flux: \f$\frac{\partial {\bf f}({\bf u})}{\partial x}\f$
   - Parabolic terms: \f$\frac{\partial^2 {\bf u}}{\partial x^2}\f$ (if applicable)
   - Source terms: \f${\bf S}({\bf u})\f$ (if applicable)

2. **Choose a reference model**: Pick an existing similar model as a template:
   - Simple hyperbolic: `Burgers` (nonlinear scalar)
   - With diffusion: `LinearADR` (linear advection-diffusion-reaction)
   - Complex system: `Euler1D` (hyperbolic system with eigenstructure)
   - With parabolic terms: `NavierStokes2D` (hyperbolic-parabolic system)

3. **Understand HyPar conventions**:
   - Solution array layout: 1D array with ghost points
   - Index calculations using `_ArrayIndex1D_` macros
   - Multi-dimensional array traversal patterns

---

## Step-by-Step Implementation

### Step 1: Create Header File

Create `include/physicalmodels/yourmodel.h`:

```c
/*! @file yourmodel.h
    @brief Your model description
    @author Your Name

    Mathematical description of your model:
    \f{equation}{
      \frac{\partial u}{\partial t} + \frac{\partial f(u)}{\partial x} = S(u)
    \f}
*/

#ifndef _YOURMODEL_H_
#define _YOURMODEL_H_

#define _YOUR_MODEL_ "yourmodel"

#include <basic.h>

/*! Structure containing model-specific parameters */
typedef struct yourmodel_parameters {

  /* Physical parameters */
  double param1;        /*!< Description of parameter 1 */
  double param2;        /*!< Description of parameter 2 */

  /* Model-specific arrays (if needed) */
  double *field;        /*!< Description of field */

  /* Flags and options */
  int option_flag;      /*!< Description of option */

} YourModel;

/* Function declarations */
int YourModelInitialize  (void*,void*);
int YourModelCleanup     (void*);

#endif
```

### Step 2: Create Source Directory

```bash
mkdir src/PhysicalModels/YourModel
```

### Step 3: Implement Initialize Function

Create `src/PhysicalModels/YourModel/YourModelInitialize.c`:

```c
/*! @file YourModelInitialize.c
    @brief Initialize the YourModel module
    @author Your Name
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <physicalmodels/yourmodel.h>
#include <mpivars.h>
#include <hypar.h>

/* Forward declarations of required functions */
double YourModelComputeCFL    (void*,void*,double,double);
int    YourModelAdvection     (double*,double*,int,void*,double);
int    YourModelUpwind        (double*,double*,double*,double*,
                               double*,double*,int,void*,double);

/* Optional function declarations */
int    YourModelDiffusion     (double*,double*,int,void*,double);
int    YourModelSource        (double*,double*,void*,void*,double);
double YourModelComputeDiffNumber (void*,void*,double,double);

int YourModelInitialize(void *s,  /*!< Solver object of type #HyPar */
                        void *m   /*!< MPI object of type #MPIVariables */
                       )
{
  HyPar         *solver  = (HyPar*)        s;
  MPIVariables  *mpi     = (MPIVariables*) m;
  YourModel     *physics = (YourModel*)    solver->physics;
  int           ferr;

  static int count = 0;

  /* Default parameter values */
  physics->param1 = 1.0;
  physics->param2 = 0.0;
  physics->field = NULL;

  /* Read physical model parameters from physics.inp */
  if (!mpi->rank) {
    FILE *in;
    if (!count) printf("Reading physical model inputs from file \"physics.inp\".\n");
    in = fopen("physics.inp","r");
    if (in) {
      char word[_MAX_STRING_SIZE_];
      ferr = fscanf(in,"%s",word); if (ferr != 1) return(1);
      if (!strcmp(word, "begin")) {
        while (strcmp(word, "end")) {
          ferr = fscanf(in,"%s",word); if (ferr != 1) return(1);
          if (!strcmp(word, "param1")) {
            ferr = fscanf(in,"%lf",&physics->param1);
            if (ferr != 1) return(1);
          } else if (!strcmp(word, "param2")) {
            ferr = fscanf(in,"%lf",&physics->param2);
            if (ferr != 1) return(1);
          } else if (strcmp(word,"end")) {
            char useless[_MAX_STRING_SIZE_];
            ferr = fscanf(in,"%s",useless);
            if (ferr != 1) return(ferr);
            printf("Warning: keyword %s in file \"physics.inp\" with value %s not ",
                   word, useless);
            printf("recognized or extraneous. Ignoring.\n");
          }
        }
      } else {
        fprintf(stderr,"Error: Illegal format in file \"physics.inp\".\n");
        return(1);
      }
      fclose(in);
    } else {
      if (!count) {
        fprintf(stderr,"Warning: File \"physics.inp\" not found. ");
        fprintf(stderr,"Using default values.\n");
      }
    }
  }

  /* Broadcast parameters to all processes */
#ifndef serial
  IERR MPIBroadcast_double(&physics->param1,1,0,&mpi->world); CHECKERR(ierr);
  IERR MPIBroadcast_double(&physics->param2,1,0,&mpi->world); CHECKERR(ierr);
#endif

  /* Check for incompatible options */
  if (!strcmp(solver->SplitHyperbolicFlux,"yes")) {
    if (!mpi->rank) {
      fprintf(stderr,"Error: YourModel does not have flux splitting defined.\n");
    }
    return(1);
  }

  /* Print parameters */
  if (!mpi->rank) {
    printf("YourModel: Parameter 1 = %lf\n", physics->param1);
    printf("YourModel: Parameter 2 = %lf\n", physics->param2);
  }

  /* Allocate model-specific arrays if needed */
  /* Example: physics->field = (double*) calloc(...); */

  /* Register function pointers with solver */
  solver->ComputeCFL         = YourModelComputeCFL;
  solver->FFunction          = YourModelAdvection;
  solver->Upwind             = YourModelUpwind;

  /* Optional functions - only set if implemented */
  /* solver->GFunction          = YourModelDiffusion; */
  /* solver->SFunction          = YourModelSource; */
  /* solver->ComputeDiffNumber  = YourModelComputeDiffNumber; */

  count++;
  return(0);
}
```

### Step 4: Implement Required Functions

#### 4.1 Compute CFL (`YourModelComputeCFL.c`)

```c
/*! @file YourModelComputeCFL.c
    @brief Compute maximum CFL number
*/

#include <basic.h>
#include <arrayfunctions.h>
#include <physicalmodels/yourmodel.h>
#include <mpivars.h>
#include <hypar.h>

double YourModelComputeCFL(void   *s,   /*!< Solver object */
                           void   *m,   /*!< MPI object */
                           double dt,   /*!< Time step */
                           double t     /*!< Current time */
                          )
{
  HyPar     *solver  = (HyPar*)     s;
  YourModel *physics = (YourModel*) solver->physics;

  int ndims  = solver->ndims;
  int nvars  = solver->nvars;
  int ghosts = solver->ghosts;
  int *dim   = solver->dim_local;
  double *u  = solver->u;

  double max_cfl = 0;
  int done = 0;
  int index[ndims];

  _ArraySetValue_(index,ndims,0);
  while (!done) {
    int p;
    _ArrayIndex1D_(ndims,dim,index,ghosts,p);

    for (int v=0; v<nvars; v++) {
      for (int dir=0; dir<ndims; dir++) {
        double dxinv;
        _GetCoordinate_(dir,index[dir],dim,ghosts,solver->dxinv,dxinv);

        /* Compute local wave speed */
        double wave_speed = /* physics-dependent calculation */;
        double local_cfl = fabs(wave_speed) * dt * dxinv;

        if (local_cfl > max_cfl) max_cfl = local_cfl;
      }
    }
    _ArrayIncrementIndex_(ndims,dim,index,done);
  }

  return(max_cfl);
}
```

#### 4.2 Hyperbolic Flux (`YourModelAdvection.c`)

```c
/*! @file YourModelAdvection.c
    @brief Compute hyperbolic flux
*/

#include <stdlib.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <physicalmodels/yourmodel.h>
#include <hypar.h>

int YourModelAdvection(double *f,    /*!< Computed flux */
                       double *u,    /*!< Solution */
                       int    dir,   /*!< Spatial dimension */
                       void   *s,    /*!< Solver object */
                       double t      /*!< Current time */
                      )
{
  HyPar     *solver  = (HyPar*)     s;
  YourModel *physics = (YourModel*) solver->physics;

  int *dim   = solver->dim_local;
  int ghosts = solver->ghosts;
  int ndims  = solver->ndims;
  int nvars  = solver->nvars;

  int index[ndims], bounds[ndims], offset[ndims];

  /* Set bounds to include ghost points */
  _ArrayCopy1D_(dim,bounds,ndims);
  for (int i=0; i<ndims; i++) bounds[i] += 2*ghosts;

  /* Set offset for ghost point arrangement */
  _ArraySetValue_(offset,ndims,-ghosts);

  /* Loop over all points */
  int done = 0;
  _ArraySetValue_(index,ndims,0);
  while (!done) {
    int p;
    _ArrayIndex1DWO_(ndims,dim,index,offset,ghosts,p);

    /* Compute flux for each variable */
    for (int v=0; v<nvars; v++) {
      /* Your flux calculation here */
      f[nvars*p+v] = /* physics-dependent flux function */;
    }

    _ArrayIncrementIndex_(ndims,bounds,index,done);
  }

  return(0);
}
```

#### 4.3 Upwinding Scheme (`YourModelUpwind.c`)

```c
/*! @file YourModelUpwind.c
    @brief Upwinding scheme
*/

#include <stdlib.h>
#include <math.h>
#include <basic.h>
#include <arrayfunctions.h>
#include <mathfunctions.h>
#include <physicalmodels/yourmodel.h>
#include <hypar.h>

int YourModelUpwind(double *fI,   /*!< Interface flux */
                    double *fL,   /*!< Left flux */
                    double *fR,   /*!< Right flux */
                    double *uL,   /*!< Left solution */
                    double *uR,   /*!< Right solution */
                    double *u,    /*!< Solution */
                    int    dir,   /*!< Spatial dimension */
                    void   *s,    /*!< Solver object */
                    double t      /*!< Current time */
                   )
{
  HyPar     *solver  = (HyPar*)     s;
  YourModel *physics = (YourModel*) solver->physics;

  int ndims  = solver->ndims;
  int nvars  = solver->nvars;
  int *dim   = solver->dim_local;
  int ghosts = solver->ghosts;

  int index_outer[ndims], index_inter[ndims];
  int bounds_outer[ndims], bounds_inter[ndims];

  _ArrayCopy1D_(dim,bounds_outer,ndims);
  bounds_outer[dir] = 1;
  _ArrayCopy1D_(dim,bounds_inter,ndims);
  bounds_inter[dir] += 1;

  int done = 0;
  _ArraySetValue_(index_outer,ndims,0);

  while (!done) {
    _ArrayCopy1D_(index_outer,index_inter,ndims);

    for (index_inter[dir]=0; index_inter[dir]<bounds_inter[dir]; index_inter[dir]++) {
      int p;
      _ArrayIndex1D_(ndims,bounds_inter,index_inter,0,p);

      int indexL[ndims], indexR[ndims];
      _ArrayCopy1D_(index_inter,indexL,ndims);
      indexL[dir]--;
      _ArrayCopy1D_(index_inter,indexR,ndims);

      int pL, pR;
      _ArrayIndex1D_(ndims,dim,indexL,ghosts,pL);
      _ArrayIndex1D_(ndims,dim,indexR,ghosts,pR);

      /* Compute upwind flux for each variable */
      for (int v=0; v<nvars; v++) {

        /* Calculate wave speeds (eigenvalues) */
        double lambdaL = /* left wave speed */;
        double lambdaR = /* right wave speed */;

        /* Simple upwinding: Roe averaging or local Lax-Friedrichs */
        if ((lambdaL > 0) && (lambdaR > 0)) {
          /* Both waves moving right - use left flux */
          fI[nvars*p+v] = fL[nvars*p+v];
        } else if ((lambdaL < 0) && (lambdaR < 0)) {
          /* Both waves moving left - use right flux */
          fI[nvars*p+v] = fR[nvars*p+v];
        } else {
          /* Mixed or near-zero - use dissipative flux */
          double alpha = max(fabs(lambdaL), fabs(lambdaR));
          fI[nvars*p+v] = 0.5 * (fL[nvars*p+v] + fR[nvars*p+v]
                               - alpha * (uR[nvars*p+v] - uL[nvars*p+v]));
        }
      }
    }

    _ArrayIncrementIndex_(ndims,bounds_outer,index_outer,done);
  }

  return(0);
}
```

#### 4.4 Cleanup Function (`YourModelCleanup.c`)

```c
/*! @file YourModelCleanup.c
    @brief Cleanup
*/

#include <stdlib.h>
#include <physicalmodels/yourmodel.h>

int YourModelCleanup(void *s /*!< Solver object */)
{
  YourModel *physics = (YourModel*) s;

  /* Free any allocated arrays */
  if (physics->field) free(physics->field);

  return(0);
}
```

### Step 5: Register Model with HyPar

Edit `src/Simulation/InitializePhysics.c`:

1. Add include at the top:
```c
#include <physicalmodels/yourmodel.h>
```

2. Add initialization block in `InitializePhysics()` function:
```c
} else if (!strcmp(solver->model,_YOUR_MODEL_)) {

  solver->physics = (YourModel*) calloc (1,sizeof(YourModel));
  IERR YourModelInitialize(solver,mpi); CHECKERR(ierr);

} else {
```

### Step 6: Update Build System

Edit `src/PhysicalModels/YourModel/Makefile.am`:

```makefile
noinst_LIBRARIES = libyourmodel.a
libyourmodel_a_SOURCES = \
  YourModelInitialize.c \
  YourModelCleanup.c \
  YourModelAdvection.c \
  YourModelUpwind.c \
  YourModelComputeCFL.c
```

Then add to parent `Makefile.am` in `src/PhysicalModels/`:
```makefile
SUBDIRS = ... YourModel
```

### Step 7: Reconfigure and Rebuild

```bash
autoreconf -i
./configure
make
```

---

## Required Functions

Every physical model **must** implement:

| Function | Purpose | Registered As |
|----------|---------|---------------|
| `YourModelInitialize` | Allocate and initialize | Called from `InitializePhysics` |
| `YourModelCleanup` | Free memory | Called from physics cleanup |
| `YourModelComputeCFL` | Calculate CFL number | `solver->ComputeCFL` |
| `YourModelAdvection` | Compute hyperbolic flux | `solver->FFunction` |
| `YourModelUpwind` | Compute interface flux | `solver->Upwind` |

---

## Optional Functions

Implement as needed for your model:

### Parabolic Terms (Diffusion)

```c
int YourModelDiffusion(double *f,    /*!< Parabolic flux */
                       double *u,    /*!< Solution */
                       int    dir,   /*!< Direction */
                       void   *s,    /*!< Solver */
                       double t      /*!< Time */
                      );
```
Register as: `solver->GFunction = YourModelDiffusion;`

### Source Terms

```c
int YourModelSource(double *source,  /*!< Source term */
                    double *u,       /*!< Solution */
                    void   *s,       /*!< Solver */
                    void   *m,       /*!< MPI */
                    double t         /*!< Time */
                   );
```
Register as: `solver->SFunction = YourModelSource;`

### Diffusion Number

```c
double YourModelComputeDiffNumber(void   *s,
                                  void   *m,
                                  double dt,
                                  double t);
```
Register as: `solver->ComputeDiffNumber = YourModelComputeDiffNumber;`

### Eigenstructure (for characteristic-based methods)

```c
int YourModelLeftEigenvectors(double *u, double *L, void *s, int dir);
int YourModelRightEigenvectors(double *u, double *R, void *s, int dir);
```
Register as: `solver->GetLeftEigenvectors` and `solver->GetRightEigenvectors`

### Pre/Post Processing

- `PreStep`: Called before each time step
- `PostStep`: Called after each time step
- `PreStage`: Called before each RK stage
- `PostStage`: Called after each RK stage
- `PrintStep`: Custom output during simulation

---

## Testing and Validation

### 1. Create Test Case

In `Examples/YourModel/`, create:
- `solver.inp`: Grid, time-stepping, spatial discretization
- `physics.inp`: Model-specific parameters
- `initial.inp`: Initial conditions
- `boundary.inp`: Boundary conditions (if not periodic)

Example `solver.inp`:
```
begin
  ndims             1
  nvars             1
  size              1 100
  iproc             1 1
  ghost             3
  n_iter            1000
  dt                0.001
  screen_op_iter    100
  file_op_iter      100
  input_mode        serial
  ip_file_type      ascii
  output_mode       serial
  op_file_format    text
  op_overwrite      yes
  model             yourmodel
end
```

### 2. Verify Basic Functionality

```bash
cd Examples/YourModel
../../bin/HyPar
```

Check:
- No segmentation faults
- CFL number is reasonable
- Solution files are created
- Conservation/stability properties

### 3. Test Edge Cases

- Very small/large time steps
- Different grid resolutions
- Different spatial schemes (WENO, MUSCL, etc.)
- MPI parallel runs: `mpiexec -n 4 ../../bin/HyPar`

### 4. Validate Against Known Solutions

- Compare with exact solutions (if available)
- Order of accuracy tests (refine grid, check convergence rate)
- Benchmark against published results

---

## Example: Implementing a Simple Model

Let's implement the **linear advection equation**:
\f$\frac{\partial u}{\partial t} + a \frac{\partial u}{\partial x} = 0\f$

### Header File (`linearadvection.h`)

```c
#define _LINEAR_ADVECTION_ "linearadvection"

typedef struct linearadvection_parameters {
  double a;  /* advection speed */
} LinearAdvection;

int LinearAdvectionInitialize(void*,void*);
int LinearAdvectionCleanup(void*);
```

### Key Implementation Details

**Flux function**: \f$f(u) = a \cdot u\f$
```c
f[nvars*p+v] = physics->a * u[nvars*p+v];
```

**Upwinding**: Simple upwind based on sign of `a`
```c
if (physics->a > 0) {
  fI[nvars*p+v] = fL[nvars*p+v];
} else {
  fI[nvars*p+v] = fR[nvars*p+v];
}
```

**CFL**: \f$CFL = |a| \Delta t / \Delta x\f$
```c
double wave_speed = physics->a;
double local_cfl = fabs(wave_speed) * dt * dxinv;
```

---

## Common Pitfalls

1. **Index calculations**: Always use HyPar macros (`_ArrayIndex1D_`, etc.)
2. **Ghost points**: Remember to include/exclude appropriately
3. **MPI broadcasting**: Broadcast parameters read from file
4. **Memory management**: Free all allocated arrays in cleanup
5. **Function registration**: Verify all function pointers are set
6. **Array layouts**: Solution array is 1D: `u[nvars*p+v]` where `p` is point index, `v` is variable
7. **Flux array sizes**: Interface arrays have one extra point in the flux direction

---

## Additional Resources

- HyPar documentation: http://hypar.github.io/
- Existing models: `src/PhysicalModels/` directory
- Array functions: `include/arrayfunctions.h`
- MPI functions: `include/mpivars.h`

---

## Checklist

Before submitting your model:

- [ ] Header file created with proper documentation
- [ ] All required functions implemented
- [ ] Model registered in `InitializePhysics.c`
- [ ] Makefile.am updated
- [ ] Code compiles without warnings
- [ ] Example test case runs successfully
- [ ] CFL/stability verified
- [ ] Validated against known solution (if available)
- [ ] Code documented with Doxygen comments
- [ ] MPI parallel execution tested

---

**Questions?** Check existing models in `src/PhysicalModels/` or contact the HyPar development team.

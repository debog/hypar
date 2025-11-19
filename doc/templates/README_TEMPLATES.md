# Template Files for Adding New Physical Models

This directory contains template files to help you implement a new physical model in HyPar.

## Quick Start

1. **Copy templates to your model directory:**
   ```bash
   mkdir src/PhysicalModels/YourModel
   cp doc/templates/TemplateModel*.c src/PhysicalModels/YourModel/
   cp doc/templates/Makefile.am src/PhysicalModels/YourModel/
   ```

2. **Copy and customize header file:**
   ```bash
   cp doc/templates/template_model.h include/physicalmodels/yourmodel.h
   ```

3. **Global search and replace:**
   - Replace `TemplateModel` with `YourModel`
   - Replace `template_model` with `yourmodel`
   - Replace `TEMPLATE_MODEL` with `YOUR_MODEL`
   - Replace `[REPLACE]` tags with your actual code/documentation
   - Replace `[YOUR NAME]` with your name

4. **Implement physics-specific code:**
   - Look for `[REPLACE]` comments in each file
   - Implement flux functions, upwinding, CFL calculation
   - Add parameter reading in Initialize function
   - Allocate/free arrays in Initialize/Cleanup

5. **Register with HyPar:**
   - Add `#include <physicalmodels/yourmodel.h>` to `src/Simulation/InitializePhysics.c`
   - Add initialization block in `InitializePhysics()` function
   - Update `src/PhysicalModels/Makefile.am` to include your subdirectory

6. **Build and test:**
   ```bash
   autoreconf -i
   ./configure
   make
   ```

## Template Files

### Required Files

| File | Purpose |
|------|---------|
| `template_model.h` | Header with model struct and function declarations |
| `TemplateModelInitialize.c` | Read parameters, allocate memory, register functions |
| `TemplateModelAdvection.c` | Compute hyperbolic flux f(u) |
| `TemplateModelUpwind.c` | Compute interface flux using upwinding |
| `TemplateModelComputeCFL.c` | Compute CFL number for stability |
| `TemplateModelCleanup.c` | Free allocated memory |
| `Makefile.am` | Build system configuration |

### Optional Files (implement as needed)

- `TemplateModelDiffusion.c` - Parabolic (diffusion) terms
- `TemplateModelSource.c` - Source terms
- `TemplateModelComputeDiffNumber.c` - Diffusion stability number
- `TemplateModelEigen.c` - Eigenvalues and eigenvectors (for characteristic methods)
- `TemplateModelPreStep.c` - Called before each time step
- `TemplateModelPostStep.c` - Called after each time step
- `TemplateModelFunctions.c` - Helper/utility functions

## Key Sections to Customize

### 1. Header File (`template_model.h`)
- Define unique model identifier: `#define _YOUR_MODEL_ "yourmodel"`
- Add physical parameters to struct
- Declare all function prototypes

### 2. Initialize Function
- Set default parameter values
- Read parameters from `physics.inp` file
- Broadcast parameters to all MPI processes
- Validate parameters
- Allocate arrays
- Register function pointers with solver

### 3. Advection Function (Flux)
- Implement flux function f(u) for your PDE
- Loop over all grid points including ghosts
- Handle multi-dimensional and multi-variable cases

### 4. Upwind Function
- Calculate wave speeds (eigenvalues)
- Choose upwinding scheme (simple upwind, Lax-Friedrichs, Roe, HLLC, etc.)
- Combine left and right biased fluxes

### 5. CFL Computation
- Calculate maximum wave speed at each point
- Compute CFL = |wave_speed| * dt / dx
- Return maximum over local domain

### 6. Cleanup Function
- Free all allocated arrays
- Clean up any temporary storage

## Important Notes

1. **Array Indexing**: HyPar uses 1D arrays with macro-based indexing
   - Use `_ArrayIndex1D_` to convert multi-dimensional indices to 1D
   - Solution array format: `u[nvars*p + v]` where p=point, v=variable

2. **Ghost Points**:
   - Interior arrays include ghosts, interface arrays do not
   - Use `_ArrayIndex1D_` for interior points with ghosts
   - Use `_ArrayIndex1D_` with offset=0 for interface arrays

3. **MPI Communication**:
   - Read parameters only on rank 0
   - Broadcast to all processes before use
   - Return local values in CFL computation (global reduction done by solver)

4. **Function Registration**:
   - Set function pointers in Initialize: `solver->FFunction = YourModelAdvection;`
   - NULL pointers indicate feature not implemented

5. **Testing**:
   - Create test case in `Examples/YourModel/`
   - Start with simple 1D case
   - Verify conservation, stability, convergence

## References

- See `doc/Adding_Physical_Models.md` for detailed developer guide
- Study existing models in `src/PhysicalModels/` for examples
- HyPar documentation:
  - https://hypar.readthedocs.io/en/latest/
  - http://hypar.github.io/

## Example Models to Study

| Model | Complexity | Features |
|-------|-----------|----------|
| `Burgers` | Simple | Scalar nonlinear hyperbolic |
| `LinearADR` | Medium | Advection-diffusion-reaction |
| `Euler1D` | Medium | System with eigenstructure |
| `NavierStokes2D` | Complex | Hyperbolic-parabolic system |
| `ShallowWater1D` | Medium | Source terms, modified solution |

## Common Patterns

### Reading a double parameter:
```c
} else if (!strcmp(word, "parameter_name")) {
  ferr = fscanf(in,"%lf",&physics->parameter_name);
  if (ferr != 1) return(1);
```

### Reading an integer flag:
```c
} else if (!strcmp(word, "flag_name")) {
  ferr = fscanf(in,"%d",&physics->flag_name);
  if (ferr != 1) return(1);
```

### Reading a string option:
```c
} else if (!strcmp(word, "option_name")) {
  ferr = fscanf(in,"%s",physics->option_string);
  if (ferr != 1) return(1);
```

### Broadcasting parameters:
```c
#ifndef serial
IERR MPIBroadcast_double(&physics->param,1,0,&mpi->world); CHECKERR(ierr);
IERR MPIBroadcast_integer(&physics->flag,1,0,&mpi->world); CHECKERR(ierr);
#endif
```

### Allocating arrays:
```c
int size = solver->npoints_local_wghosts * solver->nvars;
physics->array = (double*) calloc(size, sizeof(double));
```

## Troubleshooting

**Problem**: Segmentation fault on startup
- Check that all parameters are initialized before use
- Verify array allocations are correct size
- Ensure MPI broadcast matches variable types

**Problem**: Unstable simulation
- Verify CFL calculation is correct
- Check upwinding scheme provides dissipation
- Ensure flux function is implemented correctly

**Problem**: Wrong results
- Validate flux function against analytical expression
- Check array indexing (especially ghost points)
- Verify boundary conditions are appropriate

**Problem**: Compile errors
- Check all includes are present
- Verify function signatures match declarations
- Ensure model is registered in InitializePhysics.c

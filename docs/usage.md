# Usage Guide

This guide explains how to set up and run simulations with HyPar.

## Overview

HyPar requires several input files to run a simulation. These files specify:
- Simulation parameters (grid size, time stepping, numerical methods)
- Boundary conditions
- Initial conditions
- Physical model parameters

**Note:** It is best to start with the provided examples in the `Examples/` directory, understand their input files based on this guide, and then modify them for new cases.

## Running HyPar

### Basic Execution

For a serial build:

```bash
cd path/to/run/directory
/path/to/hypar/bin/HyPar
```

For a parallel (MPI) build:

```bash
cd path/to/run/directory
mpiexec -n 4 /path/to/hypar/bin/HyPar
```

where `4` is the number of MPI processes.

### With PETSc Time Integration

If using PETSc time integrators, create a `.petscrc` file in the run directory:

```bash
cd path/to/run/directory
/path/to/hypar/bin/HyPar
```

or specify PETSc options on the command line:

```bash
/path/to/hypar/bin/HyPar -use-petscts -ts_type rk -ts_rk_type 3
```

## Required Input Files

### 1. solver.inp (Mandatory)

Specifies main simulation parameters.

**Format:**
```
begin
    <keyword>   <value>
    <keyword>   <value>
    ...
end
```

**Key Parameters:**

| Keyword | Type | Description | Default |
|---------|------|-------------|---------|
| `ndims` | int | Number of spatial dimensions | 1 |
| `nvars` | int | Number of solution components | 1 |
| `size` | int[ndims] | Global grid size in each dimension | **must specify** |
| `iproc` | int[ndims] | MPI processes in each dimension | **must specify** |
| `ghost` | int | Number of ghost points | 1 |
| `n_iter` | int | Number of time iterations | 0 |
| `restart_iter` | int | Iteration to restart from | 0 |
| `time_scheme` | string | Time integration scheme | `euler` |
| `time_scheme_type` | string | Time scheme type | `none` |
| `hyp_space_scheme` | string | Spatial scheme for hyperbolic term | `1` |
| `hyp_flux_split` | string | Split hyperbolic flux? (`yes`/`no`) | `no` |
| `hyp_interp_type` | string | Interpolation type | `characteristic` |
| `par_space_type` | string | Parabolic discretization type | `nonconservative-1stage` |
| `par_space_scheme` | string | Spatial scheme for parabolic term | `2` |
| `dt` | double | Time step size | 0.0 |
| `conservation_check` | string | Check conservation? (`yes`/`no`) | `no` |
| `screen_op_iter` | int | Iterations between screen output | 1 |
| `file_op_iter` | int | Iterations between file output | 1000 |
| `op_file_format` | string | Output format (`text`/`binary`/`tecplot`) | `text` |
| `ip_file_type` | string | Input file type | `ascii` |
| `input_mode` | string | Input mode (`serial`/`parallel`/`mpi-io`) | `serial` |
| `output_mode` | string | Output mode (`serial`/`parallel`/`mpi-io`) | `serial` |
| `op_overwrite` | string | Overwrite output files? (`yes`/`no`) | `no` |
| `model` | string | Physical model name | **must specify** |
| `immersed_body` | string | Immersed body STL filename | `none` |
| `use_gpu` | string | Use GPU? (`yes`/`no`) | `no` |
| `gpu_device_no` | int | GPU device number | -1 |

**Example:**
```
begin
    ndims              2
    nvars              3
    size               128 128
    iproc              2 2
    ghost              3
    n_iter             1000
    time_scheme        rk
    time_scheme_type   ssprk3
    hyp_space_scheme   weno5
    dt                 0.001
    conservation_check no
    screen_op_iter     10
    file_op_iter       100
    op_file_format     text
    op_overwrite       no
    model              navierstokes2d
end
```

**Notes:**
- `ndims` **must** be specified **before** `size` and `iproc`
- For parallel I/O modes, specify the number of I/O ranks: `input_mode parallel 4`

### 2. boundary.inp (Mandatory)

Specifies boundary conditions.

**Format:**
```
nb
boundary_type   dimension   face   [extent_coords]
boundary_type   dimension   face   [extent_coords]
...
```

Where:
- `nb`: Number of boundary conditions
- `boundary_type`: Type (e.g., `periodic`, `extrapolate`, `dirichlet`, `slip-wall`, etc.)
- `dimension`: Spatial dimension (0, 1, ..., ndims-1)
- `face`: Face index (1 = left/lower, -1 = right/upper)
- `[extent_coords]`: Optional spatial extent: `xmin_0 xmax_0 xmin_1 xmax_1 ...`

**Available Boundary Types:**
- `periodic` - Periodic boundaries
- `extrapolate` - Extrapolation boundary
- `dirichlet` - Dirichlet (fixed value) boundary
- `slip-wall` - Slip wall (inviscid)
- `noslip-wall` - No-slip wall (viscous)
- `subsonic-inflow` - Subsonic inflow (specify density, velocity)
- `subsonic-outflow` - Subsonic outflow (specify pressure)
- `supersonic-inflow` - Supersonic inflow (specify all variables)
- `supersonic-outflow` - Supersonic outflow (extrapolation)

**Example (1D periodic):**
```
2
periodic    0    1
periodic    0   -1
```

**Example (2D with walls):**
```
4
slip-wall    1    1    0.0 1.0  0.0 0.0
0.0 0.0
slip-wall    1   -1    0.0 1.0  1.0 1.0
0.0 0.0
periodic     0    1
periodic     0   -1
```

**Special Boundary Requirements:**

Some boundaries require additional data on the next line:

- **dirichlet**: Specify boundary values: `u[0] u[1] ... u[nvars-1]`
- **slip-wall** / **noslip-wall**: Specify wall velocity: `u[0] u[1] ... u[ndims-1]`
- **subsonic-inflow**: Specify: `rho u[0] u[1] ... u[ndims-1]`
- **subsonic-outflow**: Specify: `p`
- **supersonic-inflow**: Specify: `rho u[0] u[1] ... u[ndims-1] p`

### 3. initial.inp (Mandatory)

Contains the initial solution. This file is typically generated programmatically, not created by hand.

**Format:** Depends on `input_mode` and `ip_file_type` specified in `solver.inp`

**Generation:** See the `aux/init.c` files in the example directories for code to generate this file. The pattern is:
1. Define the initial solution analytically
2. Evaluate it on the grid
3. Write to file in the appropriate format

### 4. physics.inp (Model-dependent)

Contains physics-specific parameters. Requirements depend on the physical model.

**Format:**
```
begin
    <keyword>   <value>
    <keyword>   <value>
    ...
end
```

**Example (Euler equations):**
```
begin
    gamma    1.4
    upwinding    roe
end
```

**Example (Navier-Stokes):**
```
begin
    gamma      1.4
    upwinding  roe
    Pr         0.72
    Re         1000.0
    Minf       0.3
end
```

Refer to the physical model documentation for specific keywords for each model.

## Optional Input Files

### 5. weno.inp (Optional)

Parameters for WENO schemes (relevant only if using WENO spatial discretization).

**Format:**
```
begin
    <keyword>   <value>
    ...
end
```

| Keyword | Type | Description | Default |
|---------|------|-------------|---------|
| `mapped` | int | Use mapped WENO? (0/1) | 0 |
| `borges` | int | Use Borges-style mapping? (0/1) | 0 |
| `yc` | int | Use Yamaleev-Carpenter mapping? (0/1) | 0 |
| `no_limiting` | int | Disable limiting? (0/1) | 0 |
| `epsilon` | double | Small parameter to avoid division by zero | 1e-6 |
| `rc` | double | CRWENO parameter | 0.3 |
| `xi` | double | CRWENO parameter | 0.001 |
| `tol` | double | Tolerance | 1e-16 |

**Example:**
```
begin
    mapped     1
    borges     0
    epsilon    1e-6
end
```

### 6. lusolver.inp (Optional)

Parameters for LU solvers for tridiagonal systems (needed for some CRWENO schemes).

**Format:**
```
begin
    <keyword>   <value>
    ...
end
```

| Keyword | Type | Description | Default |
|---------|------|-------------|---------|
| `evaluate_norm` | int | Evaluate norm? (0/1) | 1 |
| `maxiter` | int | Maximum iterations | 10 |
| `atol` | double | Absolute tolerance | 1e-12 |
| `rtol` | double | Relative tolerance | 1e-10 |
| `verbose` | int | Verbosity level | 0 |

### 7. .petscrc (Optional, PETSc required)

PETSc time integration parameters. If this file exists (or flags are specified on command line), PETSc time integration is used.

**Format:** Standard PETSc options format

**HyPar-specific flags:**
- `-use-petscts` - Enable PETSc time integration
- `-jfnk_epsilon <value>` - Epsilon for Jacobian approximation (default: 1e-6)
- `-with_pc` - Use preconditioning

**IMEX-specific flags** (for implicit-explicit methods):
- `-hyperbolic_explicit` / `-hyperbolic_implicit`
- `-parabolic_explicit` / `-parabolic_implicit` (default: implicit)
- `-source_explicit` / `-source_implicit` (default: implicit)

**Example:**
```
-use-petscts
-ts_type arkimex
-ts_arkimex_type 3
-ts_adapt_type none
-ts_max_time 1.0
-ts_dt 0.001
-parabolic_implicit
-hyperbolic_explicit
```

### 8. Immersed Body STL File (Optional)

Required only if using immersed boundaries. Filename specified by `immersed_body` in `solver.inp`.

**Format:** ASCII STL format

**Requirements:**
- Normals must point **outward** from the body
- Geometry must be closed
- Sample files available in `Examples/STLGeometries/`

### 9. simulation.inp (Optional)

Required for ensemble/multidomain simulations.

**Format:**
```
begin
    nsims    <number_of_simulations>
end
```

### 10. sparse_grids.inp (Optional)

Required for sparse grids method.

**Format:**
```
begin
    log2_imin           2
    interp_order        6
    write_sg_solution   no
    write_sg_errors     no
end
```

### 11. librom.inp (Optional)

Required for reduced-order modeling with libROM.

**Key Parameters:**

| Keyword | Type | Description | Default |
|---------|------|-------------|---------|
| `rdim` | int | ROM dimension | **must specify** |
| `sampling_frequency` | int | Sampling frequency | 1 |
| `mode` | string | `train` or `predict` | `train` |
| `type` | string | ROM type (e.g., `DMD`) | `DMD` |
| `save_to_file` | string | Save ROM to file? | `true` |

**Example:**
```
begin
    rdim                 10
    sampling_frequency   1
    mode                 train
    type                 DMD
    save_to_file         true
end
```

## Output Files

After running, HyPar generates:

- **op_XXXXX.dat** - Solution files (numbered by output iteration)
- **errors.dat** - Numerical errors (if exact solution provided)
- **conservation.dat** - Conservation error (if `conservation_check` is `yes`)
- **function_counts.dat** - Function evaluation counts
- Model-specific output files (varies by physical model)

## Quick Start Example

Here's a minimal setup for a 1D linear advection problem:

**solver.inp:**
```
begin
    ndims        1
    nvars        1
    size         100
    iproc        1
    ghost        3
    n_iter       100
    time_scheme  rk
    time_scheme_type  ssprk3
    hyp_space_scheme  weno5
    dt           0.01
    model        linearadr
end
```

**boundary.inp:**
```
2
periodic  0   1
periodic  0  -1
```

**physics.inp:**
```
begin
    advection  1.0
end
```

Then generate `initial.inp` using a custom initialization program and run:

```bash
/path/to/hypar/bin/HyPar
```

## Tips and Best Practices

1. **Start with examples** - The `Examples/` directory contains working setups
2. **Check screen output** - HyPar prints useful diagnostic information
3. **Verify conservation** - Use `conservation_check yes` for conservation laws
4. **Test with exact solutions** - When available, use them to verify accuracy
5. **Parallel I/O** - For large simulations, use `parallel` or `mpi-io` modes
6. **GPU usage** - Set `use_gpu yes` if compiled with CUDA support

For detailed examples, see the [Examples](examples.md) section.
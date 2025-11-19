# Examples

HyPar includes over 100 comprehensive example test cases organized by dimensionality and physical model. Each example is a complete, ready-to-run simulation with:

- All required input files (`solver.inp`, `physics.inp`, `boundary.inp`)
- Initialization code (`aux/init.c`) to generate initial conditions
- MATLAB/Python visualization scripts
- Expected output and solution images
- Detailed documentation

This page provides an overview of available examples with representative results. For complete details on each example, see the [full Doxygen documentation](../doc/html/examples.html).

## Example Directory Structure

All examples are located in the `Examples/` directory:

```
Examples/
├── 1D/              # One-dimensional test cases
├── 2D/              # Two-dimensional test cases  
├── 3D/              # Three-dimensional test cases
├── 4D/              # Four-dimensional test cases
├── LaSDI/           # Latent space dynamics identification examples
├── Python/          # Python scripts for visualization
├── Matlab/          # MATLAB scripts for visualization
└── STLGeometries/   # STL geometry files for immersed boundary method
```

## Categories of Examples

### Basic Examples

Explicit time integration examples that do not require external libraries. These can be run with the base HyPar installation.

### PETSc Examples

Examples demonstrating implicit and IMEX (implicit-explicit) time integration using PETSc. Requires HyPar to be compiled with PETSc support (`--with-petsc`).

### GPU Examples

Examples that leverage NVIDIA GPUs for acceleration. Requires HyPar compiled with CUDA support (`--enable-cuda`).

### Immersed Boundary Examples

Simulations with complex geometries using the immersed boundary method. Requires STL geometry files.

### libROM Examples

Reduced-order modeling examples using libROM library for model reduction and acceleration. Requires HyPar compiled with libROM support.

### Multidomain Examples

Examples demonstrating multidomain functionality for coupled simulations.

### Sparse Grids Examples

High-dimensional problems solved using sparse grids spatial discretization.

---

## 1D Examples

### Linear Advection-Diffusion-Reaction

#### Sine Wave Advection

**Location:** `Examples/1D/LinearAdvection/SineWave`

Advects a sine wave with constant speed using 5th order CRWENO spatial discretization and SSPRK3 time integration.

**Initial condition:** $u(x,0) = \sin(2\pi x)$ on $0 \le x < 1$ (periodic)

**Features:**
- Exact solution available (periodic, exact match after one period)
- Tests spatial accuracy of high-order schemes
- Numerical errors computed and reported

![1D Linear Advection Sine Wave](../doc/html/Solution_1DLinearAdvSine.png)

---

#### Sine Wave with Spatially-Varying Advection Speed

**Location:** `Examples/1D/LinearAdvection/SineWave_NonConstantAdvection`

Tests spatially-varying advection field: $a(x) = 1 + \frac{1}{2}\sin^2(2\pi x)$

**Features:**
- Demonstrates handling of non-constant coefficients
- Writes out advection field at each output time
- Uses binary file input for advection field

![Varying Advection Animation](../doc/html/Solution_1DLinearAdvSine_VaryingAdv.gif)
![Advection Field](../doc/html/Solution_1DLinearAdvSine_VaryingAdv.png)

---

#### Discontinuous Waves

**Location:** `Examples/1D/LinearAdvection/DiscontinuousWaves`

Advects four different discontinuous features simultaneously:
- Smooth Gaussian bell
- Square wave
- Triangle wave  
- Semi-circle

**Features:**
- Tests shock-capturing and oscillation control
- CRWENO scheme maintains sharp discontinuities
- Multiple wave types in single simulation

![Discontinuous Waves](../doc/html/Solution_1DLinearAdvDisc.png)

**Reference:** Ghosh & Baeder (2012), "Compact Reconstruction Schemes with Weighted ENO Limiting"

---

#### Linear Diffusion - Sine Wave

**Location:** `Examples/1D/LinearDiffusion/SineWave`

Pure diffusion problem: $\frac{\partial u}{\partial t} = \nu \frac{\partial^2 u}{\partial x^2}$

**Initial condition:** $u(x,0) = \sin(2\pi x)$

**Exact solution:** $u(x,t) = e^{-4\pi^2 \nu t}\sin(2\pi x)$

**Features:**
- Tests parabolic term discretization
- 2nd order conservative scheme
- Exponential decay validation

![Linear Diffusion](../doc/html/Solution_1DLinearDiffSine.png)

---

### Burgers Equation

#### Sine Wave

**Location:** `Examples/1D/Burgers/SineWave`

Nonlinear Burgers equation: $\frac{\partial u}{\partial t} + \frac{\partial}{\partial x}\left(\frac{1}{2}u^2\right) = 0$

**Initial condition:** $u(x,0) = \frac{1}{2\pi t_s}\sin(2\pi x)$ where $t_s=2$ is time to shock formation

**Features:**
- Tests shock-capturing for nonlinear hyperbolic equation
- Solution evolution before shock formation  
- Demonstrates steepening of wave profile

![Burgers Sine Wave](../doc/html/Solution_1DBurgersSine.png)

---

### 1D Euler Equations

#### Sod Shock Tube

**Location:** `Examples/1D/Euler1D/SodShockTube`

Classic 1D Riemann problem with left and right states:
- Left ($x < 0.5$): $\rho=1, u=0, p=1$
- Right ($x \ge 0.5$): $\rho=0.125, u=0, p=0.1$

**Governing equations:** 1D compressible Euler equations with $\gamma=1.4$

**Features:**
- Contact discontinuity, shock wave, and expansion fan
- Tests WENO shock-capturing capability
- Benchmark problem for compressible flow solvers

![Sod Shock Tube](../doc/html/Solution_1DSodShockTube.png)

**Reference:** Sod, G.A. (1978), "A survey of several finite difference methods for systems of nonlinear hyperbolic conservation laws"

---

#### Lax Shock Tube

**Location:** `Examples/1D/Euler1D/LaxShockTube`

Strong shock test case with conserved variable initialization:
- Left: $\rho=0.445, \rho u=0.311, e=8.928$
- Right: $\rho=0.5, \rho u=0, e=1.4275$

**Features:**
- Tests robustness for strong shocks
- Characteristic-based WENO scheme
- Run with 2 MPI ranks (parallel example)

![Lax Shock Tube](../doc/html/Solution_1DLaxShockTube.png)

**Reference:** Lax, P.D. (1954), "Weak solutions of nonlinear hyperbolic equations"

---

#### Shu-Osher Problem

**Location:** `Examples/1D/Euler1D/ShuOsherProblem`

Shock-density wave interaction problem. Tests ability to capture both shocks and smooth waves.

**Initial conditions:**
- Left ($x < -4$): Mach 3 shock: $\rho=27/7, u=4\sqrt{35}/7, p=31/3$
- Right ($x \ge -4$): Density wave: $\rho=1+0.2\sin(5x), u=0, p=1$

**Features:**
- High-frequency waves behind shock
- Tests resolution of high-order schemes
- Challenging problem for numerical methods

![Shu-Osher Problem](../doc/html/Solution_1DShuOsherProblem.png)

**Reference:** Shu, C.-W. & Osher, S. (1989), "Efficient implementation of essentially non-oscillatory schemes"

---

#### Sod Shock Tube with Gravitational Force

**Location:** `Examples/1D/Euler1D/SodShockTubeWithGravity`

Same as Sod problem but with uniform gravitational force $g=1$ and slip-wall boundaries.

**Features:**
- Tests well-balanced schemes for gravitational source terms
- Source term splitting and upwinding
- Preserves hydrostatic equilibrium

![Sod with Gravity](../doc/html/Solution_1DSodShockTubeWithGravity.png)

**Reference:** Xing & Shu (2013), "High Order Well-Balanced WENO Scheme for Gas Dynamics Under Gravitational Fields"

---

### 1D Shallow Water Equations

#### Dam Break over Rectangular Bump

**Location:** `Examples/1D/ShallowWater1D/DamBreakingRectangularBump`

Dam break problem over a rectangular bottom topography bump.

**Initial conditions:**
- $h(x) = \begin{cases}20-b(x) & x \le 750 \\ 15-b(x) & x > 750\end{cases}$, $u(x)=0$
- Bottom topography: $b(x) = \begin{cases}8 & |x-750| \le 1500/8 \\ 0 & \text{otherwise}\end{cases}$

**Features:**
- Tests well-balanced source term treatment
- Preserves "lake at rest" steady state exactly
- Demonstrates wetting/drying over topography

![Dam Break over Bump](../doc/html/Solution_1DSWDamBreak.png)

**Reference:** Xing & Shu (2005), "High order finite difference WENO schemes with exact conservation property for shallow water equations"

---

## 2D Examples

### Linear Advection-Diffusion

#### Gaussian Pulse Advection

**Location:** `Examples/2D/LinearAdvection/GaussianPulse`

Advects a 2D Gaussian pulse: $u(x,y,0) = \exp\left[-(x^2+y^2)/2\right]$ with constant velocity $(a_x, a_y) = (1, 1)$.

**Domain:** $-6 \le x,y < 6$ with periodic boundaries

**Features:**
- Smooth solution for error analysis
- 5th order CRWENO spatial discretization
- Parallel example (8 MPI ranks: 4×2)

![2D Gaussian Pulse](../doc/html/Solution_2DLinearAdvGauss.gif)

---

#### Sine Wave with Rotating Velocity Field

**Location:** `Examples/2D/LinearAdvection/SineWave_NonConstantAdvection/case_01`

Tests spatially-varying advection field with rotation:
- $a_x(x,y) = \sin(4\pi y)$
- $a_y(x,y) = -\cos(4\pi x)$

**Initial condition:** $u(x,y,0) = \cos(4\pi y)$

![Rotating Velocity Field Animation](../doc/html/Solution_2DLinearAdvSine_VaryingAdv.gif)
![Velocity Field Vectors](../doc/html/Solution_2DLinearAdvSine_VaryingAdv.png)

---

#### 2D Linear Diffusion - Sine Wave

**Location:** `Examples/2D/LinearDiffusion/SineWave`

Pure 2D diffusion with exact solution:
- Initial: $u(x,y,0) = \sin(2\pi x)\sin(2\pi y)$
- Exact: $u(x,y,t) = \exp[-\pi^2(4\nu_x+4\nu_y)t]\sin(2\pi x)\sin(2\pi y)$

![2D Diffusion](../doc/html/Solution_2DLinearDiffSine.gif)

---

### 2D Burgers Equation

#### Sine Wave  

**Location:** `Examples/2D/Burgers/SineWave`

2D nonlinear Burgers: $\frac{\partial u}{\partial t} + \frac{\partial}{\partial x}\left(\frac{u^2}{2}\right) + \frac{\partial}{\partial y}\left(\frac{u^2}{2}\right) = 0$

**Solution evolution:**

| t=0 (initial) | t=0.4 | t=0.8 |
|---------------|-------|-------|
| ![](../doc/html/Solution_2DBurgersSineWave_0.png) | ![](../doc/html/Solution_2DBurgersSineWave_1.png) | ![](../doc/html/Solution_2DBurgersSineWave_2.png) |

| t=1.2 | t=1.6 (final) |
|-------|---------------|
| ![](../doc/html/Solution_2DBurgersSineWave_3.png) | ![](../doc/html/Solution_2DBurgersSineWave_4.png) |

---

### 2D Vlasov Equation (1D-1V)

#### Two-Stream Instability

**Location:** `Examples/2D/Vlasov1D1V/TwoStreamInstability`

**Requires:** FFTW library for self-consistent electric field computation

Kinetic plasma simulation with two counter-streaming beams. The distribution function evolves in 2D phase space (x, v).

**Initial condition:** Two Gaussian beams centered at $v=\pm 2$ with small spatial perturbation

**Features:**
- Self-consistent Poisson solve via FFT: $\nabla^2\phi = -\rho_e$
- Electric field: $E = -\nabla\phi$
- Tests phase-space advection and field coupling

![Two-Stream Instability](../doc/html/Solution_1D1VVlasov_TwoStreamInstability.gif)
![Electric Field Evolution](../doc/html/Solution_E_1D1VVlasov_TwoStreamInstability.png)

---

#### Landau Damping

**Location:** `Examples/2D/Vlasov1D1V/LandauDamping`

**Requires:** FFTW library

Classic collisionless damping of plasma waves. Initial perturbation decays due to phase mixing.

**Initial:** $f(x,v,0) = \frac{1}{\sqrt{2\pi v_{th}^2}}\exp\left(-\frac{v^2}{2v_{th}^2}\right)[1 + \alpha\cos(kx)]$

**Features:**
- Benchmark for kinetic solvers
- Electric field amplitude decays exponentially  
- Validates phase-space resolution

![Landau Damping Rate](../doc/html/Solution_E_1D1VVlasov_LandauDamping.png)

**Reference:** Finn et al. (2023), "A Numerical Study of Landau Damping with PETSc-PIC"

---

### 2D Euler Equations

#### Isentropic Vortex Convection

**Location:** `Examples/2D/NavierStokes2D/InviscidVortexConvection`

Smooth vortex convection in uniform flow. Tests accuracy without shocks.

**Setup:**
- Freestream: $\rho_\infty=1, u_\infty=0.1, v_\infty=0, p_\infty=1$
- Vortex perturbation with strength $b=0.5$ added at center

**Features:**
- Exact solution available (vortex returns after one period)
- Tests dispersion and dissipation errors
- 5th order CRWENO + SSPRK3

![Isentropic Vortex](../doc/html/Solution_2DNavStokVortex.gif)

**Reference:** Shu (1997), "Essentially Non-oscillatory and Weighted Essentially Non-oscillatory Schemes"

---

#### 2D Riemann Problem (Case 4)

**Location:** `Examples/2D/NavierStokes2D/Riemann2DCase4`

Four-quadrant Riemann problem with complex wave interactions.

**Features:**
- Multiple shocks, contact discontinuities, and rarefactions  
- Tests 2D shock-capturing
- Characteristic-based WENO scheme

![Riemann Case 4](../doc/html/Solution_2DNavStokRiemann4.png)

**Reference:** Lax & Liu (1998), "Solution of two-dimensional Riemann problems by positive schemes"

---

#### Rising Thermal Bubble (Euler)

**Location:** `Examples/2D/NavierStokes2D/RisingThermalBubble`

Warm bubble rises due to buoyancy in stratified atmosphere.

**Domain:** $1000\times 1000$ m with slip-wall boundaries

**Features:**
- Gravitational source term with well-balanced scheme
- Hydrostatic equilibrium preservation (HB type 2)
- Atmospheric flow simulation

![Rising Thermal Bubble](../doc/html/Solution_2DNavStokRTB.png)
![Animation](../doc/html/Solution_2DNavStokRTB.gif)

**Reference:** Giraldo & Restelli (2008), "Spectral element and DG methods for Navier-Stokes equations in atmospheric modeling"

---

### 2D Navier-Stokes Equations

#### Lid-Driven Square Cavity

**Location:** `Examples/2D/NavierStokes2D/LidDrivenCavity`

Classic benchmark for incompressible viscous flow. Square cavity with moving top lid.

**Reynolds numbers tested:** Re = 100, 1000, 3200

**Features:**
- Tests viscous term discretization and boundary conditions
- Steady-state solution at various Reynolds numbers  
- Comparison with benchmark data (Erturk et al., 2005)

| Re = 100 | Re = 1000 | Re = 3200 |
|----------|-----------|-----------|
| ![](../doc/html/Solution_2DNavStokLDSC_Re0100.png) | ![](../doc/html/Solution_2DNavStokLDSC_Re1000.png) | ![](../doc/html/Solution_2DNavStokLDSC_Re3200.png) |

**Reference:** Erturk, Corke & Gokcol (2005), "Numerical Solutions of 2-D Steady Incompressible Driven Cavity Flow at High Reynolds Numbers"

---

#### Laminar Flow Past Flat Plate

**Location:** `Examples/2D/NavierStokes2D/FlatPlateBoundaryLayer`

Compressible laminar boundary layer development over flat plate.

**Features:**
- Tests boundary layer resolution
- Skin friction coefficient computation  
- Comparison with Blasius solution

![Flat Plate Flow](../doc/html/Solution_2DNavStokFlatPlate.png)
![Magnified View](../doc/html/Solution_2DNavStokFlatPlateMagnified.png)
![Skin Friction](../doc/html/Solution_2DNavStokFlatPlateSkinFriction.png)

---

### 2D Shallow Water Equations

| Example | Location | Description |
|---------|----------|-------------|
| **Circular Dam Break** | `2D/ShallowWater2D/CircularDamBreak` | Radially symmetric dam break |
| **Flow Over Bump** | `2D/ShallowWater2D/LatitudinalBeltFlow` | Geophysical flow test |

---

## 3D Examples

### 3D Navier-Stokes Equations

#### Direct Numerical Simulation of Isotropic Turbulence

**Location:** `Examples/3D/NavierStokes3D/DNS_IsotropicTurbulence`

High-resolution DNS of decaying isotropic turbulence on $512^3$ grid.

**Features:**
- Tests 3D solver performance and scaling
- Turbulent kinetic energy spectrum computation
- Energy decay monitoring
- GPU acceleration available

![Energy Spectrum](../doc/html/Solution_3DNavStok_IsoTurb_Spectrum.png)
![Energy Decay](../doc/html/Solution_3DNavStok_IsoTurb_Energy.png)
![Density Isosurface](../doc/html/Solution_3DNavStok_IsoTurb.png)

---

#### Rising Thermal Bubble (3D)

**Location:** `Examples/3D/NavierStokes3D/RisingThermalBubble_terrain`

3D atmospheric convection with terrain-following coordinates.

**Domain:** $1000\times 1000\times 1000$ m

**Features:**
- 3D well-balanced atmospheric dynamics
- Terrain-following coordinate system
- Gravitational source terms

![3D Bubble](../doc/html/Solution_3DNavStok_Bubble3D.png)
![Animation](../doc/html/Solution_3DNavStok_Bubble.gif)

---

#### Flow Past Sphere with Immersed Boundary

**Location:** `Examples/3D/NavierStokes3D/SphereViscous_Adiabatic`

Viscous flow past sphere using immersed boundary method.

**Reynolds numbers:** Re_D = 100, 1000

**Features:**
- STL geometry input (sphere)
- No-slip boundary conditions enforced on immersed surface
- Surface forces and heat flux computation
- Benchmark validation data

![Sphere Domain](../doc/html/Surface3D_Sphere.png)
![ReD=100 Flow](../doc/html/Solution_3DNavStokSphere_ReD100.png)
![Surface Pressure](../doc/html/IBSurface_3DNavStokSphereAdiabatic_Pressure.png)

---

### 3D NUMA

| Example | Location | Description |
|---------|----------|-------------|
| **Rising Thermal Bubble** | `3D/Numa3D/RisingThermalBubble` | Non-hydrostatic atmospheric model |

---

## 4D Examples

High-dimensional test cases for sparse grids and uncertainty quantification:

- **4D Advection-Diffusion**: Multi-dimensional transport
- **Fokker-Planck Equations**: Probability distribution evolution

---

## Running an Example

### Basic Workflow

1. **Navigate to example directory:**
   ```bash
   cd Examples/1D/LinearAdvection/SineWave
   ```

2. **Generate initial condition:**
   ```bash
   gcc -o aux/init aux/init.c -lm
   ./aux/init
   ```
   This creates `initial.inp` (and possibly `exact.inp` for error analysis).

3. **Run HyPar:**
   ```bash
   ../../../../bin/HyPar
   ```
   For parallel runs:
   ```bash
   mpiexec -n 4 ../../../../bin/HyPar
   ```

4. **Visualize results:**
   - Use provided MATLAB script: `matlab < Run.m`
   - Use Python scripts in `Examples/Python/`
   - Use Tecplot/VisIt for Tecplot format outputs
   - Plot text files directly

### Common Input Files

Every example requires:

- **`solver.inp`**: Numerical method parameters (scheme, time step, output frequency)
- **`physics.inp`**: Physical model parameters (e.g., gamma, viscosity, gravity)
- **`boundary.inp`**: Boundary condition specifications
- **`initial.inp`**: Initial solution (generated by `aux/init.c`)

Optional files:
- **`weno.inp`**: WENO scheme parameters
- **`lusolver.inp`**: Compact scheme parameters  
- **`exact.inp`**: Exact solution for error calculation
- **`topography.inp`**: Bottom topography (shallow water)
- **`advection.inp`**: Spatially-varying advection field

### Example: 1D Linear Advection Sine Wave

**Location:** `Examples/1D/LinearAdvection/SineWave`

**Quick start:**
```bash
cd Examples/1D/LinearAdvection/SineWave
gcc -o aux/init aux/init.c -lm
./aux/init
../../../../bin/HyPar
```

**Output:** 11 files `op_00000.dat` through `op_00010.dat` containing the solution at different times.

**Visualization (if MATLAB is available):**
```bash
matlab < Run.m
```

---

## Visualiza tion Tools

### MATLAB Scripts

Many examples include `Run.m` which:
- Generates initial conditions
- Runs HyPar
- Plots results

Usage:
```bash
matlab < Run.m
```

### Python Scripts

Located in `Examples/Python/`:
- `plotSolution1D.py`: Plot 1D solution files
- `plotSolution2D.py`: Plot 2D solution files
- `animate2D.py`: Create animations from solution sequence

### Tecplot/VisIt

For `op_file_format tecplot2d` or `tecplot3d`, open `.dat` files directly in Tecplot or VisIt.

---

## Advanced Examples

### GPU Acceleration

**Requirements:**
- CUDA-capable GPU
- HyPar compiled with `--enable-cuda`

**Examples:**
- `2D/NavierStokes2D/*_GPU/`
- `3D/NavierStokes3D/*_GPU/`

**Usage:** Same as regular examples; GPU utilization is automatic.

### Immersed Boundary Method

**Requirements:**
- STL geometry file
- `immersed_body <geometry.stl>` in `solver.inp`

**Example:**
```bash
cd Examples/2D/NavierStokes2D/FlowPastCylinder_IB
# STL file: ../../STLGeometries/cylinder.stl
../../../../bin/HyPar
```

**Visualizing IB:**
- Surface data written to `surface*.dat`
- Forces written to `forces*.dat`

### Reduced-Order Modeling with libROM

**Requirements:**
- HyPar compiled with libROM (`--with-librom`)

**Workflow:**
1. **Train ROM:**
   ```bash
   mpiexec -n 4 ../../../../bin/HyPar librom.inp train
   ```
   Generates ROM database files.

2. **Predict with ROM:**
   ```bash
   ../../../../bin/HyPar librom.inp predict
   ```
   Fast prediction using ROM.

**Examples:**
- `2D/NavierStokes2D/*_libROM_DMD/`
- `3D/NavierStokes3D/*_libROM_DMD/`

### Sparse Grids

**Requirements:**
- Multidimensional problem (typically 3D+)

**Example:**
```bash
cd Examples/SparseGrids/
../../../../bin/HyPar
```

Uses combination technique to reduce computational cost in high dimensions.

---

## Example Output

### Standard Output Files

- **`op_*.dat`**: Solution at output times
  - Format controlled by `op_file_format` in `solver.inp`
  - Text, binary, or Tecplot formats

- **`errors.dat`**: Numerical errors (if exact solution provided)
  - Columns: grid size, processors, dt, L1 error, L2 error, L∞ error, solver time, total time

- **`conservation.dat`**: Conservation errors
  - Monitors conservation properties

- **`function_counts.dat`**: Function evaluation counts (profiling)

### Model-Specific Output

- **Shallow Water:** `topography_*.dat`
- **Vlasov:** `efield_*.dat`, `potential_*.dat`
- **LinearADR:** `advection_*.dat` (for spatially-varying advection)
- **Immersed Boundary:** `surface*.dat`, `forces*.dat`

---

## Testing and Verification

### Convergence Tests

Run examples with different grid resolutions to verify spatial accuracy order:

```bash
cd Examples/ConvergenceTests/
./run_convergence_test.sh
```

Generates convergence plots showing error vs. grid spacing.

### Benchmark Cases

Many examples have exact or reference solutions for verification:
- Sod shock tube
- Isentropic vortex
- Sine waves with diffusion
- Landau damping

Compare `op_*.dat` with `exact.inp` or reference data.

---

## Tips and Best Practices

### Running Examples Efficiently

1. **Start small:** Begin with coarse grids to verify setup
2. **Check conservation:** Enable `conservation_check yes` in `solver.inp`
3. **Use MPI:** Examples with `iproc` > 1 require MPI
4. **Monitor CFL:** Check screen output for CFL number
5. **Visualize early:** Plot initial condition to verify setup

### Modifying Examples

To adapt an example:

1. **Copy the directory:**
   ```bash
   cp -r Examples/1D/LinearAdvection/SineWave MyTest/
   ```

2. **Modify input files:** Edit `solver.inp`, `physics.inp`, `boundary.inp`

3. **Update initialization:** Modify `aux/init.c` and recompile

4. **Run and compare:** Test against original example

### Common Issues

| Issue | Solution |
|-------|----------|
| **Segmentation fault** | Check grid size matches in all input files |
| **CFL violation** | Reduce time step size in `solver.inp` |
| **Non-convergence** | Try more robust scheme (e.g., WENO instead of upwind) |
| **Parallel I/O errors** | Check `iproc` matches MPI ranks |
| **Missing initial.inp** | Compile and run `aux/init.c` |

---

## Additional Resources

### Example Documentation

Detailed documentation for each example is available:
- **Doxygen HTML docs:** `doc/html/examples.html` (after building docs)
- **Markdown file:** `doc/Examples.md`

### Helper Scripts

- **`Examples/BashTools/`**: Bash utilities for batch runs
- **`Examples/Python/`**: Python visualization scripts  
- **`Examples/Matlab/`**: MATLAB plotting tools
- **`Examples/SpectralAnalysis/`**: Fourier analysis tools

### Online Resources

- **GitHub Repository:** https://github.com/debog/hypar
- **Main Documentation:** http://hypar.github.io
- **Issue Tracker:** Report problems or ask questions

---

## Summary

HyPar's extensive example suite provides:
- **Verification cases** with exact solutions
- **Benchmark problems** from literature
- **Application examples** for various physics
- **Tutorial cases** for learning features
- **Scaling tests** for performance analysis

Start with simple 1D examples and progressively explore more complex cases as you become familiar with HyPar's capabilities.

**Next:** [Physical Models](physical_models.md) | [Numerical Methods](numerical_methods.md) | [Usage Guide](usage.md)

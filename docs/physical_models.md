# Physical Models

HyPar provides implementations of several physical models representing different types of hyperbolic-parabolic PDEs. This section provides comprehensive documentation for each implemented model, including governing equations, parameters, and usage guidelines.

## Overview

Physical models in HyPar define:
- **Governing equations** - Mathematical formulation of the PDE system
- **Flux functions** - Hyperbolic (advection) and parabolic (diffusion) terms
- **Source terms** - External forces and reaction terms
- **Upwinding schemes** - Numerical flux computation methods
- **Auxiliary functions** - Model-specific calculations (eigenvalues, transformations, etc.)

## Available Models

### 1D Models
- [Linear Advection-Diffusion-Reaction](#linear-advection-diffusion-reaction-linearadr)
- [Burgers Equation](#burgers-equation)
- [1D Euler Equations](#1d-euler-equations)
- [1D Shallow Water Equations](#1d-shallow-water-equations)

### 2D Models
- [2D Euler Equations](#2d-euler-equations)
- [2D Navier-Stokes Equations](#2d-navier-stokes-equations)
- [2D Shallow Water Equations](#2d-shallow-water-equations)
- [2D NUMA (Non-hydrostatic Unified Model)](#2d-numa)
- [Vlasov Equation (1D-1V)](#vlasov-equation)

### 3D Models
- [3D Navier-Stokes Equations](#3d-navier-stokes-equations)
- [3D NUMA](#3d-numa)

### Special Models
- [Fokker-Planck Models](#fokker-planck-models)

---

## Linear Advection-Diffusion-Reaction (LinearADR)

**Model identifier:** `linear-advection-diffusion-reaction`

### Governing Equations

$$
\frac{\partial u}{\partial t} + \sum_d \frac{\partial}{\partial x_d}(a_d u) = \sum_d \nu_d \frac{\partial^2 u}{\partial x_d^2} + k u
$$

where:
- $u$ is the scalar solution variable
- $a_d$ are the advection speeds in each direction
- $\nu_d$ are the diffusion coefficients in each direction
- $k$ is the reaction rate

### Features
- **Dimensions:** Arbitrary (1D, 2D, 3D, etc.)
- **Variables:** 1 (scalar equation)
- **Type:** Linear hyperbolic-parabolic with source term
- **GPU Support:** No

### Parameters

Configure in `physics.inp`:

```
begin
  advection     <a_x> <a_y> <a_z> ...  # Advection speeds (one per dimension)
  diffusion     <d_x> <d_y> <d_z> ...  # Diffusion coefficients
  reaction      <k>                     # Reaction rate (optional, default: 0)
  constant_advection  <0 or 1>         # 1=constant, 0=spatially varying
  adv_filename  <filename>             # File for spatially-varying advection field
  centered_flux <yes/no>               # Use centered flux (no upwinding)
end
```

**Example:**
```
begin
  advection     1.0              # 1D advection with speed 1.0
  diffusion     0.01             # Diffusion coefficient
  reaction      0.0              # No reaction
  constant_advection  1          # Constant advection
  centered_flux no               # Use upwinding
end
```

### Spatially-Varying Advection

For spatially-varying advection fields, set `constant_advection 0` and provide a binary file containing the advection field at each grid point.

**File format:** Binary file with values at all grid points (including ghosts):
```
[a_0_x, a_0_y, a_0_z, a_1_x, a_1_y, a_1_z, ...]
```

### Upwinding Schemes

Available upwinding schemes (set in `solver.inp`):
- First-order upwind
- WENO3, WENO5
- CRWENO5, HCWENO5

### Use Cases

- Transport problems
- Pollution dispersion
- Heat conduction with advection
- Chemical reactions with transport

### Example Test Cases

See `Examples/` directory:
- `1D/LinearAdvection/` - 1D advection test cases
- `2D/LinearAdvection/` - 2D advection test cases
- `1D/LinearDiffusion/` - 1D diffusion test cases

---

## Burgers Equation

**Model identifier:** `burgers`

### Governing Equations

$$
\frac{\partial u}{\partial t} + \sum_{i=1}^{d} \frac{\partial}{\partial x_i}\left(\frac{1}{2}u^2\right) = 0
$$

### Features
- **Dimensions:** Arbitrary (typically 1D or 2D)
- **Variables:** 1 (scalar)
- **Type:** Nonlinear hyperbolic
- **GPU Support:** No

### Parameters

No parameters required in `physics.inp`. Create an empty file:

```
begin
end
```

### Characteristics

- Simplest nonlinear hyperbolic equation
- Models shock formation and propagation
- Often used for testing numerical schemes
- Admits exact solutions for many initial conditions

### Use Cases

- Testing shock-capturing schemes
- Verification of high-order methods
- Educational purposes
- Traffic flow modeling

### Example Test Cases

- `Examples/1D/Burgers/` - Various initial conditions
  - Sine wave
  - Traveling wave
  - Shock formation

---

## 1D Euler Equations

**Model identifier:** `euler1d`

### Governing Equations

$$
\frac{\partial}{\partial t}\begin{bmatrix}\rho \\ \rho u \\ e\end{bmatrix}
+ \frac{\partial}{\partial x}\begin{bmatrix}\rho u \\ \rho u^2 + p \\ (e+p)u\end{bmatrix}
= \begin{bmatrix}0 \\ -\rho g \\ -\rho u g\end{bmatrix}
$$

**Equation of state:**
$$
e = \frac{p}{\gamma - 1} + \frac{1}{2}\rho u^2
$$

where:
- $\rho$ = density
- $u$ = velocity
- $e$ = total energy per unit volume
- $p$ = pressure
- $\gamma$ = ratio of specific heats
- $g$ = gravitational acceleration

### Features
- **Dimensions:** 1
- **Variables:** 3 (density, momentum, energy)
- **Type:** Hyperbolic with source term (gravity)
- **GPU Support:** Yes
- **Supports:** Well-balanced schemes for gravitational source terms

### Parameters

Configure in `physics.inp`:

```
begin
  gamma        <value>        # Ratio of specific heats (typically 1.4 for air)
  gravity      <value>        # Gravitational acceleration (default: 0.0)
  grav_type    <0, 1, or 2>   # Type of gravitational field
  upwinding    <scheme>       # Upwinding scheme choice
end
```

**Gravitational field types:**
- `0`: Constant gravitational field
- `1`: Isothermal atmosphere
- `2`: Sinusoidal potential

**Upwinding schemes:**
- `roe` - Roe's approximate Riemann solver
- `rf-char` - Roe-Fixed scheme with entropy fix
- `llf-char` - Local Lax-Friedrichs
- `steger-warming` - Steger-Warming flux splitting
- `rusanov` - Rusanov scheme

### Example Configuration

```
begin
  gamma        1.4            # Ideal gas
  gravity      0.0            # No gravity
  upwinding    roe            # Roe scheme
end
```

### Flux Partitioning (IMEX)

For implicit-explicit time integration, the flux can be split into:
- **Stiff flux** (acoustic modes) - integrated implicitly
- **Non-stiff flux** (entropy mode) - integrated explicitly

This enables efficient time-stepping for flows with disparate time scales.

### Well-Balanced Schemes

For flows with gravity, HyPar implements well-balanced schemes that preserve hydrostatic equilibrium to machine precision.

**Reference:** Xing & Shu (2013), "High Order Well-Balanced WENO Scheme for the Gas Dynamics Equations Under Gravitational Fields"

### Use Cases

- Compressible gas dynamics
- Shock tube problems
- Astrophysical flows with gravity
- Acoustic wave propagation
- Supersonic flows

### Example Test Cases

- `Examples/1D/Euler1D/` directory:
  - Sod shock tube
  - Lax problem
  - Shu-Osher problem
  - Blast wave interactions
  - Gravitational flows

---

## 2D Euler Equations

**Model identifier:** `euler2d`

### Governing Equations

$$
\frac{\partial}{\partial t}\begin{bmatrix}\rho \\ \rho u \\ \rho v \\ e\end{bmatrix}
+ \frac{\partial}{\partial x}\begin{bmatrix}\rho u \\ \rho u^2 + p \\ \rho uv \\ (e+p)u\end{bmatrix}
+ \frac{\partial}{\partial y}\begin{bmatrix}\rho v \\ \rho uv \\ \rho v^2 + p \\ (e+p)v\end{bmatrix} = \mathbf{0}
$$

**Equation of state:**
$$
e = \frac{p}{\gamma - 1} + \frac{1}{2}\rho(u^2 + v^2)
$$

### Features
- **Dimensions:** 2
- **Variables:** 4 (density, x-momentum, y-momentum, energy)
- **Type:** Hyperbolic
- **GPU Support:** No

### Parameters

```
begin
  gamma        <value>        # Ratio of specific heats
  upwinding    <scheme>       # Upwinding scheme
end
```

**Upwinding schemes:** Same as 1D Euler

### Use Cases

- 2D shock interactions
- Vortex flows
- Rayleigh-Taylor instability
- Richtmyer-Meshkov instability
- Supersonic flow over bodies

### Example Test Cases

- `Examples/2D/Euler2D/`:
  - Isentropic vortex
  - Riemann problems
  - Double Mach reflection
  - Forward-facing step
  - Shock-vortex interaction

---

## 2D Navier-Stokes Equations

**Model identifier:** `navierstokes2d`

### Governing Equations

$$
\frac{\partial}{\partial t}\begin{bmatrix}\rho \\ \rho u \\ \rho v \\ e\end{bmatrix}
+ \frac{\partial}{\partial x}\begin{bmatrix}\rho u \\ \rho u^2 + p \\ \rho uv \\ (e+p)u\end{bmatrix}
+ \frac{\partial}{\partial y}\begin{bmatrix}\rho v \\ \rho uv \\ \rho v^2 + p \\ (e+p)v\end{bmatrix}
$$
$$
= \frac{\partial}{\partial x}\begin{bmatrix}0 \\ \tau_{xx} \\ \tau_{yx} \\ u\tau_{xx} + v\tau_{yx} - q_x\end{bmatrix}
+ \frac{\partial}{\partial y}\begin{bmatrix}0 \\ \tau_{xy} \\ \tau_{yy} \\ u\tau_{xy} + v\tau_{yy} - q_y\end{bmatrix}
+ \begin{bmatrix}0 \\ -\rho g_x \\ -\rho g_y \\ -\rho ug_x - \rho vg_y\end{bmatrix}
$$

**Viscous stress tensor:**
$$
\tau_{ij} = \frac{\mu}{Re_\infty}\left[\left(\frac{\partial u_i}{\partial x_j} + \frac{\partial u_j}{\partial x_i}\right) - \frac{2}{3}\frac{\partial u_k}{\partial x_k}\delta_{ij}\right]
$$

**Heat flux:**
$$
q_i = -\frac{\mu}{(\gamma-1)Re_\infty Pr}\frac{\partial T}{\partial x_i}
$$

where $\mu$ is computed using Sutherland's law.

### Features
- **Dimensions:** 2
- **Variables:** 4
- **Type:** Hyperbolic-parabolic with source terms
- **GPU Support:** Yes
- **Special features:** 
  - Immersed boundary method support
  - Well-balanced gravitational source terms
  - Flux partitioning for IMEX

### Parameters

```
begin
  gamma        <value>        # Ratio of specific heats
  upwinding    <scheme>       # Upwinding scheme
  gravity      <gx> <gy>      # Gravitational acceleration components
  rho0         <value>        # Reference density at zero altitude
  p0           <value>        # Reference pressure at zero altitude
  Re           <value>        # Reynolds number
  Pr           <value>        # Prandtl number
  Minf         <value>        # Freestream Mach number
  HB           <0,1,2,3>      # Hydrostatic balance type
  N_bv         <value>        # Brunt-Vaisala frequency (for HB=3)
end
```

**Hydrostatic balance types:**
- `1`: Isothermal equilibrium
- `2`: Constant potential temperature
- `3`: Stratified atmosphere with Brunt-Vaisala frequency

### Non-dimensionalization

The equations are non-dimensionalized using:
- Velocity scale: speed of sound at reference conditions
- Reynolds number based on speed of sound
- Reference values: $\rho_0$, $p_0$, $T_0$

### Sutherland's Law

Dynamic viscosity:
$$
\mu(T) = C_1 \frac{T^{3/2}}{T + C_2}
$$

### Use Cases

- Viscous flows
- Boundary layer flows
- Atmospheric flows with gravity
- Thermal convection
- Rayleigh-Benard convection
- Flows with immersed boundaries

### Example Test Cases

- `Examples/2D/NavierStokes2D/`:
  - Viscous shock tube
  - Lid-driven cavity
  - Rising thermal bubble
  - Density current
  - Flow past cylinder (immersed boundary)

---

## 1D Shallow Water Equations

**Model identifier:** `shallow-water-1d`

### Governing Equations

$$
\frac{\partial}{\partial t}\begin{bmatrix}h \\ hu\end{bmatrix}
+ \frac{\partial}{\partial x}\begin{bmatrix}hu \\ hu^2 + \frac{1}{2}gh^2\end{bmatrix}
= \begin{bmatrix}0 \\ -ghb_x\end{bmatrix}
$$

where:
- $h$ = water height
- $u$ = velocity
- $g$ = gravitational constant
- $b(x)$ = bottom topography

### Features
- **Dimensions:** 1
- **Variables:** 2 (height, momentum)
- **Type:** Hyperbolic with source term
- **GPU Support:** No
- **Special features:** Well-balanced schemes for topography

### Parameters

```
begin
  gravity      <value>        # Gravitational constant (typically 9.8)
  upwinding    <scheme>       # roe or llf-char
  bottomtopo_from_file  <filename>  # Optional: topography file
end
```

### Well-Balanced Treatment

The source term is discretized to preserve the "lake at rest" steady state exactly:

$$
hu = 0, \quad h + b = \text{constant}
$$

**Reference:** Xing & Shu (2005), "High order finite difference WENO schemes with the exact conservation property for the shallow water equations"

### Topography Input

Provide topography in a binary file with values at all grid points (including ghosts).

### Use Cases

- Dam break problems
- Tsunami modeling
- River flows
- Hydraulic jumps
- Verification of well-balanced schemes

### Example Test Cases

- `Examples/1D/ShallowWater1D/`:
  - Dam break
  - Lake at rest
  - Flow over a bump
  - Hydraulic jump

---

## 2D Shallow Water Equations

**Model identifier:** `shallow-water-2d`

### Governing Equations

$$
\frac{\partial}{\partial t}\begin{bmatrix}h \\ hu \\ hv\end{bmatrix}
+ \frac{\partial}{\partial x}\begin{bmatrix}hu \\ hu^2 + \frac{1}{2}gh^2 \\ huv\end{bmatrix}
+ \frac{\partial}{\partial y}\begin{bmatrix}hv \\ huv \\ hv^2 + \frac{1}{2}gh^2\end{bmatrix}
= \begin{bmatrix}0 \\ -ghb_x \\ -ghb_y\end{bmatrix}
$$

### Features
- **Dimensions:** 2
- **Variables:** 3
- **Type:** Hyperbolic with source terms
- **GPU Support:** No
- **Special features:** Well-balanced schemes

### Use Cases

- 2D dam breaks
- Geophysical flows
- Tsunami propagation
- Tidal flows
- Flow over complex topography

---

## 3D Navier-Stokes Equations

**Model identifier:** `navierstokes3d`

### Governing Equations

3D extension of 2D Navier-Stokes with additional z-direction terms.

### Features
- **Dimensions:** 3
- **Variables:** 5 (density, 3 momenta, energy)
- **Type:** Hyperbolic-parabolic with source terms
- **GPU Support:** Yes
- **Special features:**
  - Immersed boundary method
  - Well-balanced schemes

### Use Cases

- 3D flows
- Turbulence simulations
- Atmospheric flows
- Flow past 3D bodies
- Complex geometry flows

---

## 2D NUMA

**Model identifier:** `numa2d`

### Governing Equations

2D Non-hydrostatic Unified Model of the Atmosphere:

$$
\frac{\partial}{\partial t}\begin{bmatrix}\delta\rho \\ (\rho_0+\delta\rho)u \\ (\rho_0+\delta\rho)v \\ \delta\theta\end{bmatrix}
+ \nabla\cdot\mathbf{F} = \mathbf{S}
$$

where $\delta$ denotes perturbations from a reference hydrostatic state.

### Features
- **Dimensions:** 2
- **Variables:** 4
- **Type:** Hyperbolic-parabolic
- **GPU Support:** No
- **Special features:** Designed for atmospheric flows

### Parameters

```
begin
  gamma        <value>        # Ratio of specific heats
  R            <value>        # Universal gas constant
  g            <value>        # Gravitational acceleration
  Pref         <value>        # Reference pressure
  Tref         <value>        # Reference temperature
  init_atmos   <type>         # Initial atmosphere type
  upwind       <scheme>       # Upwinding (rusanov)
end
```

### Use Cases

- Atmospheric dynamics
- Weather modeling
- Gravity waves
- Density currents
- Mountain waves

**Reference:** Giraldo, Restelli & Laeuter (2010), "Semi-Implicit Formulations of the Euler Equations"

---

## 3D NUMA

**Model identifier:** `numa3d`

3D extension of NUMA2D model.

### Features
- **Dimensions:** 3
- **Variables:** 5
- **Type:** Hyperbolic-parabolic
- **GPU Support:** No

---

## Vlasov Equation

**Model identifier:** `vlasov`

### Governing Equations

1D-1V Vlasov equation (kinetic plasma physics):

$$
\frac{\partial f}{\partial t} + v\frac{\partial f}{\partial x} + E\frac{\partial f}{\partial v} = 0
$$

where:
- $f(x,v,t)$ = distribution function
- $x$ = spatial coordinate
- $v$ = velocity coordinate
- $E(x,t)$ = electric field

### Features
- **Dimensions:** 2 (1 spatial + 1 velocity)
- **Variables:** 1 (distribution function)
- **Type:** Hyperbolic
- **GPU Support:** No
- **Special features:** Self-consistent electric field via Poisson equation

### Self-Consistent Electric Field

When enabled, solves Poisson's equation:

$$
\frac{\partial^2\phi}{\partial x^2} = -\rho_e(x), \quad E = -\frac{\partial\phi}{\partial x}
$$

where $\rho_e(x) = \int f(x,v)\,dv$ is the charge density.

**Requires:** FFTW library for FFT-based Poisson solver

### Parameters

```
begin
  self_consistent_electric_field  <0 or 1>  # Enable self-consistent E-field
  use_log_form                    <0 or 1>  # Solve in logarithmic form
end
```

### Use Cases

- Plasma physics
- Landau damping
- Two-stream instability
- Kinetic simulations
- Collisionless plasmas

**Reference:** Henon (1982), "Vlasov equation?", Astronomy and Astrophysics

---

## Fokker-Planck Models

HyPar includes several Fokker-Planck equation models for uncertainty quantification and stochastic dynamics.

### FP Double Well

**Model identifier:** `fp-double-well`

Fokker-Planck equation for a particle in a double-well potential.

### FP Power System (1-Bus)

**Model identifier:** `fp-power-system-1bus`

Fokker-Planck model for a single-bus power system with stochastic perturbations.

### FP Power System (3-Bus)

**Model identifier:** `fp-power-system-3bus`

Fokker-Planck model for a three-bus power system.

### Common Features

- Model probability distribution evolution
- Source of stochasticity
- Useful for uncertainty quantification

---

## Summary Table

| Model | Dimensions | Variables | Type | GPU | Immersed Boundary |
|-------|------------|-----------|------|-----|-------------------|
| LinearADR | Any | 1 | Hyp-Par | No | No |
| Burgers | Any | 1 | Hyp | No | No |
| Euler1D | 1 | 3 | Hyp | Yes | No |
| Euler2D | 2 | 4 | Hyp | No | No |
| NavierStokes2D | 2 | 4 | Hyp-Par | Yes | Yes |
| NavierStokes3D | 3 | 5 | Hyp-Par | Yes | Yes |
| ShallowWater1D | 1 | 2 | Hyp | No | No |
| ShallowWater2D | 2 | 3 | Hyp | No | No |
| NUMA2D | 2 | 4 | Hyp-Par | No | No |
| NUMA3D | 3 | 5 | Hyp-Par | No | No |
| Vlasov | 2 | 1 | Hyp | No | No |

---

## Selecting a Model

In `solver.inp`, specify:

```
model     <model_identifier>
```

For example:
```
model     navierstokes2d
```

Then create the corresponding `physics.inp` with model-specific parameters.

---

## Advanced Topics

### Implementing Custom Models

See the [Developer Guide](developer_guide.md) for detailed instructions on implementing new physical models.

### Flux Partitioning

Several models support splitting the flux into stiff and non-stiff components for IMEX time integration:
- Euler1D
- NavierStokes2D
- NavierStokes3D

Set in `solver.inp`:
```
use_fast_jac           yes
use_petsc_ts_type      arkimex
```

### Immersed Boundary Method

Available for:
- NavierStokes2D
- NavierStokes3D

Provide STL geometry file and enable in `solver.inp`:
```
immersed_body          <geometry.stl>
ib_write_surface_data  yes
```

### GPU Acceleration

Models with GPU support can be run on CUDA-enabled GPUs. Compile with `--enable-cuda` and use GPU-compatible time integrators.

---

## References

1. **Euler/Navier-Stokes:** Tannehill, Anderson & Pletcher, "Computational Fluid Mechanics and Heat Transfer"
2. **Well-balanced schemes:** Ghosh & Constantinescu (2018), AIAA Journal, 54(4), 1370-1385
3. **IMEX splitting:** Ghosh & Constantinescu (2016), SIAM J. Sci. Comput., 38(3), A1848-A1875
4. **Shallow water:** Xing & Shu (2005), J. Comput. Phys., 208, 206-227
5. **NUMA:** Giraldo, Restelli & Laeuter (2010), SIAM J. Sci. Comput., 32(6), 3394-3425
6. **Vlasov:** Henon (1982), Astronomy and Astrophysics, 114

---

**Next:** [Numerical Methods](numerical_methods.md) | [Usage Guide](usage.md) | [Developer Guide](developer_guide.md)

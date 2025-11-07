# Numerical Methods

HyPar solves hyperbolic-parabolic partial differential equations using conservative finite-difference methods on Cartesian grids.

## Governing Equation

HyPar solves the following PDE:

$$
\frac{\partial \mathbf{u}}{\partial t} = \mathbf{F}_{\rm hyp}(\mathbf{u}) + \mathbf{F}_{\rm par}(\mathbf{u}) + \mathbf{F}_{\rm sou}(\mathbf{u})
$$

where:
- $\mathbf{u}$ is the solution vector
- $\mathbf{F}_{\rm hyp}$ is the **hyperbolic term** (advection/convection)
- $\mathbf{F}_{\rm par}$ is the **parabolic term** (diffusion)
- $\mathbf{F}_{\rm sou}$ is the **source term**

Each term is discretized in space to obtain a semi-discrete ODE:

$$
\frac{d\mathbf{u}}{dt} = \hat{\mathbf{F}}_{\rm hyp}(\mathbf{u}) + \hat{\mathbf{F}}_{\rm par}(\mathbf{u}) + \hat{\mathbf{F}}_{\rm sou}(\mathbf{u})
$$

where $\hat{(\cdot)}$ represents spatially discretized terms.

## Spatial Discretization

### Hyperbolic Term

The hyperbolic term has the form:

$$
\mathbf{F}_{\rm hyp}(\mathbf{u}) = -\sum_{d=0}^{D-1} \frac{\partial \mathbf{f}_d(\mathbf{u})}{\partial x_d}
$$

where $d$ is the spatial dimension and $D$ is the total number of dimensions.

**Finite-difference discretization:**

$$
\mathbf{F}_{\rm hyp}(\mathbf{u}) \approx -\sum_{d=0}^{D-1} \frac{1}{\Delta x_d} \left[\hat{\mathbf{f}}_{d,j+1/2} - \hat{\mathbf{f}}_{d,j-1/2}\right]
$$

where $j$ denotes the grid index along dimension $d$.

**Interface flux reconstruction:**

The numerical flux at interfaces $\hat{\mathbf{f}}_{d,j+1/2}$ is computed using:

$$
\hat{\mathbf{f}}_{d,j+1/2} = \mathcal{U}\left(\hat{\mathbf{f}}^L_{d,j+1/2}, \hat{\mathbf{f}}^R_{d,j+1/2}, \hat{\mathbf{u}}^L_{d,j+1/2}, \hat{\mathbf{u}}^R_{d,j+1/2}\right)
$$

where:
- $\mathcal{U}$ is an upwinding function (e.g., Roe, HLLC, Rusanov)
- $\hat{\mathbf{f}}^{L,R}_{d,j+1/2}$ are left- and right-biased flux reconstructions
- $\hat{\mathbf{u}}^{L,R}_{d,j+1/2}$ are left- and right-biased solution reconstructions

**Available spatial schemes:**
- **1st order**: First-order upwind
- **2nd order**: MUSCL reconstruction
- **3rd order**: WENO3
- **5th order**: WENO5, CRWENO5, HCWENO5
- **Non-compact schemes**: Up to 5th order

### Parabolic Term

The parabolic term can take different forms depending on whether cross-derivatives are present.

#### No Cross-Derivatives

For parabolic terms without cross-derivatives:

$$
\mathbf{F}_{\rm par}(\mathbf{u}) = \sum_{d=0}^{D-1} \frac{\partial^2 \mathbf{g}_d(\mathbf{u})}{\partial x_d^2}
$$

**Conservative discretization (1-stage):**

$$
\mathbf{F}_{\rm par}(\mathbf{u}) \approx \sum_{d=0}^{D-1} \frac{1}{\Delta x_d^2} \left[\hat{\mathbf{g}}_{d,j+1/2} - \hat{\mathbf{g}}_{d,j-1/2}\right]
$$

**Non-conservative discretization (1-stage):**

$$
\mathbf{F}_{\rm par}(\mathbf{u}) \approx \sum_{d=0}^{D-1} \frac{1}{\Delta x_d^2} \left[\mathcal{L}_d(\mathbf{g}_d)\right]
$$

where $\mathcal{L}$ represents the finite-difference Laplacian operator.

#### With Cross-Derivatives

For parabolic terms with cross-derivatives:

$$
\mathbf{F}_{\rm par}(\mathbf{u}) = \sum_{d_1=0}^{D-1}\sum_{d_2=0}^{D-1} \frac{\partial^2 \mathbf{h}_{d_1,d_2}(\mathbf{u})}{\partial x_{d_1} \partial x_{d_2}}
$$

**Non-conservative 2-stage discretization:**

$$
\mathbf{F}_{\rm par}(\mathbf{u}) \approx \sum_{d_1=0}^{D-1}\sum_{d_2=0}^{D-1} \frac{1}{\Delta x_{d_1} \Delta x_{d_2}} \left[\mathcal{D}_{d_1}(\mathcal{D}_{d_2}(\mathbf{g}_d))\right]
$$

where $\mathcal{D}_d$ denotes the finite-difference first derivative operator along dimension $d$.

**Available schemes:**
- 2nd order central
- 4th order central
- Higher-order non-compact schemes

### Source Term

Source terms are discretized pointwise - no spatial derivatives are involved:

$$
\mathbf{F}_{\rm sou}(\mathbf{u}) = \mathbf{S}(\mathbf{u})
$$

The source function is provided by the specific physical model.

## Time Integration

### Native Time Integrators

HyPar implements several explicit time integration methods:

#### Forward Euler

First-order explicit method:

$$
\mathbf{u}^{n+1} = \mathbf{u}^n + \Delta t \mathbf{F}(\mathbf{u}^n)
$$

**Usage:** Set `time_scheme = euler` in `solver.inp`

#### Runge-Kutta Methods

Multi-stage explicit RK methods:

**RK2 (Midpoint):**
- `time_scheme = rk`
- `time_scheme_type = 22`

**RK4 (Classic 4th order):**
- `time_scheme = rk`
- `time_scheme_type = 44`

**SSPRK3 (Strong Stability Preserving RK3):**
- `time_scheme = rk`
- `time_scheme_type = ssprk3`

The general RK form:

$$
\begin{aligned}
\mathbf{u}^{(0)} &= \mathbf{u}^n \\
\mathbf{u}^{(i)} &= \sum_{j=0}^{i-1} \left(a_{ij} \mathbf{u}^{(j)} + \Delta t \, b_{ij} \mathbf{F}(\mathbf{u}^{(j)})\right) \\
\mathbf{u}^{n+1} &= \mathbf{u}^{(s)}
\end{aligned}
$$

where $s$ is the number of stages.

#### GLM-GEE (General Linear Methods with Global Error Estimation)

Advanced time integrators with adaptive error control:

- `time_scheme = glm_gee`
- `time_scheme_type = <method_name>`

**Available methods:**
- `exrk2a` - 2nd order, 2 stages
- `ark2a` - 2nd order IMEX
- `ars3` - 3rd order IMEX, stiffly accurate

### PETSc Time Integrators

When compiled with PETSc, HyPar can use PETSc's TS (Time Stepper) module, enabling:

#### Explicit Methods
- `TSRK` - Runge-Kutta schemes
- `TSSSP` - Strong Stability Preserving schemes

#### Implicit Methods
- `TSBEULER` - Backward Euler
- `TSCN` - Crank-Nicolson
- `TSTHETA` - Theta method

#### IMEX (Implicit-Explicit) Methods
- `TSARKIMEX` - Additive Runge-Kutta IMEX schemes

**Usage:**

Create a `.petscrc` file or use command-line flags:

```
-use-petscts
-ts_type arkimex
-ts_arkimex_type 3
-ts_adapt_type basic
-ts_max_time 1.0
-ts_dt 0.001
```

**IMEX partitioning flags:**

Control which terms are treated explicitly vs. implicitly:

| Term | Explicit Flag | Implicit Flag |
|------|---------------|---------------|
| Hyperbolic | `-hyperbolic_explicit` | `-hyperbolic_implicit` |
| Parabolic | `-parabolic_explicit` | `-parabolic_implicit` |
| Source | `-source_explicit` | `-source_implicit` |

**Example IMEX setup (treat diffusion implicitly, advection explicitly):**

```
-use-petscts
-ts_type arkimex
-ts_arkimex_type 3
-hyperbolic_explicit
-parabolic_implicit
```

#### Jacobian-Free Newton-Krylov (JFNK)

For implicit methods, HyPar uses JFNK by default:
- Jacobian action approximated via directional derivatives
- No explicit Jacobian matrix required

Control the directional derivative parameter:
```
-jfnk_epsilon 1e-6
```

#### Preconditioning

If the physical model provides a Jacobian function, preconditioning can be enabled:

```
-with_pc
-pc_type lu
```

or

```
-with_pc
-pc_type bjacobi
-sub_pc_type ilu
```

## Upwinding Methods

For hyperbolic problems, various upwinding schemes are available:

### Characteristic-Based Upwinding

Uses local characteristic decomposition:
- Roe's approximate Riemann solver
- Local Lax-Friedrichs (Rusanov)
- HLLC (Harten-Lax-van Leer Contact)

### Component-Wise Upwinding

Simpler upwinding without characteristic decomposition (for scalar equations or when characteristics are not available).

## Spatial Schemes Summary

### Hyperbolic Schemes

| Scheme | Order | Type | Keyword |
|--------|-------|------|---------|
| First Order | 1 | Upwind | `1` |
| MUSCL | 2 | TVD | `muscl` or `2` |
| WENO3 | 3 | WENO | `weno3` or `3` |
| WENO5 | 5 | WENO | `weno5` or `5` |
| CRWENO5 | 5 | Compact WENO | `crweno5` |
| HCWENO5 | 5 | Hybrid Compact WENO | `hcweno5` |

### Parabolic Schemes

| Scheme | Order | Type | Keyword |
|--------|-------|------|---------|
| 2nd Order | 2 | Central | `2` |
| 4th Order | 4 | Central | `4` |

**Configuration in solver.inp:**

```
begin
    ...
    hyp_space_scheme   weno5
    par_space_scheme   2
    ...
end
```

## Stability and CFL Condition

For explicit time integration, stability requires:

$$
\Delta t \leq \text{CFL} \times \min_d \left(\frac{\Delta x_d}{\max(|\lambda|)}\right)
$$

where:
- $\text{CFL}$ is the CFL number (typically ≤ 1 for explicit schemes)
- $\lambda$ are the eigenvalues of the Jacobian
- $\Delta x_d$ is the grid spacing in dimension $d$

HyPar computes and reports the maximum CFL during runtime.

## References

1. **WENO Schemes:**
   - Jiang, G.-S., & Shu, C.-W. (1996). Efficient implementation of weighted ENO schemes. *Journal of Computational Physics*, 126(1), 202-228.

2. **CRWENO:**
   - Ghosh, D., & Baeder, J. D. (2012). Compact reconstruction schemes with weighted ENO limiting for hyperbolic conservation laws. *SIAM Journal on Scientific Computing*, 34(3), A1678-A1706.

3. **Time Integration:**
   - Shu, C.-W., & Osher, S. (1988). Efficient implementation of essentially non-oscillatory shock-capturing schemes. *Journal of Computational Physics*, 77(2), 439-471.

4. **PETSc:**
   - Balay, S., et al. PETSc Users Manual. Argonne National Laboratory. https://petsc.org/release/docs/manual.pdf

For implementation details, see the source code in `src/HyParFunctions/` and the physical model directories.

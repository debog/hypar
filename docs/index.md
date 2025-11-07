# HyPar Documentation

Welcome to the official documentation for **HyPar** - a finite-difference solver for hyperbolic-parabolic partial differential equations on Cartesian grids.

## What is HyPar?

HyPar is a unified computational framework for solving systems of hyperbolic-parabolic PDEs with source terms. It is designed to be:

- **Flexible**: Supports arbitrary spatial dimensions and solution components
- **Modular**: Easy to add new physical models and numerical methods
- **Parallel**: MPI-based parallelization for large-scale simulations
- **Feature-rich**: Includes multiple spatial schemes (WENO, MUSCL), time integrators, and advanced features
- **Well-tested**: Extensive examples covering 1D, 2D, 3D, and 4D problems

### Key Features

- **Spatial Discretization**: 1st to 5th order schemes including WENO, CRWENO, MUSCL
- **Time Integration**: Explicit RK methods, implicit/IMEX schemes (via PETSc)
- **Physical Models**: Euler equations, Navier-Stokes, shallow water, Vlasov, and more
- **Advanced Capabilities**:
  - GPU acceleration (CUDA)
  - Immersed boundary method
  - Reduced-order modeling (libROM)
  - Sparse grids
  - Ensemble simulations
- **Parallel I/O**: Efficient parallel file operations for large datasets

### Mathematical Formulation

HyPar solves PDEs of the form:

$$
\frac{\partial \mathbf{u}}{\partial t} = \mathbf{F}_{\text{hyp}}(\mathbf{u}) + \mathbf{F}_{\text{par}}(\mathbf{u}) + \mathbf{F}_{\text{sou}}(\mathbf{u})
$$

where:
- $\mathbf{u}$ is the solution vector
- $\mathbf{F}_{\text{hyp}}$ is the hyperbolic (advection) term
- $\mathbf{F}_{\text{par}}$ is the parabolic (diffusion) term  
- $\mathbf{F}_{\text{sou}}$ is the source term

## Quick Links

```{toctree}
:maxdepth: 2

installation
usage
numerical_methods
developer_guide
api
```

## Getting Started

### For Users

1. **[Installation](installation.md)** - Build and install HyPar
2. **[Usage Guide](usage.md)** - Learn how to set up and run simulations
3. **[Numerical Methods](numerical_methods.md)** - Understand the algorithms

### For Developers

1. **[Developer Guide](developer_guide.md)** - Add new physical models
2. **[API Reference](api.md)** - Code documentation

### Example Workflow

Here's a minimal example to get started with a 1D linear advection problem:

**1. Install HyPar:**
```bash
git clone https://github.com/debog/hypar.git
cd hypar
./configure
make
```

**2. Navigate to an example:**
```bash
cd Examples/1D/LinearAdvection/SineWave
```

**3. Run the simulation:**
```bash
../../../../bin/HyPar
```

**4. Visualize results:**
Use MATLAB, Python, or any plotting tool to visualize the `op_*.dat` output files.

## Supported Physical Models

HyPar includes implementations for:

### 1D Models
- Linear advection-diffusion-reaction
- Burgers equation
- Euler equations (gas dynamics)
- Shallow water equations
- Vlasov equation (kinetic plasma physics)

### 2D Models
- Linear advection and diffusion
- 2D Burgers equation
- 2D Euler equations (compressible flow)
- 2D Navier-Stokes equations (viscous flow)
- 2D Shallow water equations

### 3D Models
- 3D Euler equations
- 3D Navier-Stokes equations
- 3D Non-hydrostatic atmospheric flow (NUMA)

## Available Spatial Schemes

| Scheme | Order | Type | Best For |
|--------|-------|------|----------|
| First-order upwind | 1 | FD | Debugging |
| MUSCL | 2 | TVD | General purpose |
| WENO3 | 3 | WENO | Smooth + shocks |
| WENO5 | 5 | WENO | High accuracy |
| CRWENO5 | 5 | Compact WENO | High accuracy (compact) |
| HCWENO5 | 5 | Hybrid Compact WENO | Very high accuracy |

## Time Integration Methods

### Native Schemes
- Forward Euler
- RK2, RK4
- Strong Stability Preserving RK3 (SSPRK3)
- General Linear Methods (GLM-GEE)

### PETSc Integration
When compiled with PETSc, access to:
- Implicit methods (backward Euler, Crank-Nicolson)
- IMEX (Implicit-Explicit) schemes
- Adaptive time-stepping
- Nonlinear solvers (Newton-Krylov)

## Examples Directory

The `Examples/` directory contains over 100 test cases organized by dimension:

- **1D/**: Linear advection, Burgers, Euler, shallow water, Vlasov
- **2D/**: 2D versions + vortex flows, shock interactions, boundary layers
- **3D/**: 3D flows including turbulence, shock-body interactions
- **4D/**: High-dimensional test cases

Each example includes:
- Complete input files (`solver.inp`, `physics.inp`, etc.)
- Initialization code (`aux/init.c`)
- Expected output or reference solutions
- MATLAB/Python visualization scripts

## Performance Features

### Parallelization
- **MPI**: Domain decomposition for distributed memory
- **OpenMP**: Shared-memory parallelization (optional)
- **GPU**: CUDA support for accelerated computations

### I/O Modes
- **Serial I/O**: Simple, single-process file operations
- **Parallel I/O**: Multiple MPI ranks write simultaneously
- **MPI-IO**: Collective I/O for best performance on large systems

### Output Formats
- **Text**: Human-readable ASCII
- **Binary**: Compact binary format
- **Tecplot**: Direct Tecplot format support

## Advanced Features

### Immersed Boundary Method
Simulate complex geometries without body-fitted grids:
- Reads STL geometry files
- Enforces boundary conditions on immersed surfaces
- Examples: flow over cylinders, spheres, airfoils

### Reduced-Order Modeling (libROM)
- Dynamic Mode Decomposition (DMD)
- Train ROMs from full-order simulations
- Accelerated predictions for parameter studies

### Sparse Grids
- Combination technique for high-dimensional problems
- Reduced computational cost
- Maintains accuracy in lower dimensions

## Building Documentation Locally

To build this documentation on your machine:

```bash
cd docs
pip install -r requirements-docs.txt
make html
```

Then open `_build/html/index.html` in your browser.

## Citation

If you use HyPar in your research, please cite:

```bibtex
@article{ghosh2018hypar,
  title={Well-balanced, conservative finite-difference algorithm for atmospheric flows},
  author={Ghosh, Debojyoti and Constantinescu, Emil M},
  journal={AIAA Journal},
  volume={56},
  number={4},
  pages={1370--1385},
  year={2018},
  publisher={American Institute of Aeronautics and Astronautics}
}
```

## License

HyPar is released under the BSD 3-Clause License. See `License.md` in the repository for details.

## Contact & Support

- **GitHub**: https://github.com/debog/hypar
- **Issues**: https://github.com/debog/hypar/issues
- **Email**: debojyoti.ghosh@gmail.com
- **Website**: http://hypar.github.io/

## Contributing

Contributions are welcome! Whether it's:
- Bug reports and fixes
- New physical models
- Documentation improvements
- Example test cases
- Performance optimizations

Please see the [Developer Guide](developer_guide.md) for information on contributing code.

---

**Version**: 4.1  
**Last Updated**: 2025
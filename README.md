HyPar - Hyperbolic-Parabolic Partial Differential Equations Solver
------------------------------------------------------------------

`HyPar` is a finite-difference algorithm to solve hyperbolic-parabolic partial differential
equations (with source terms) on Cartesian grids. It is a unified framework that can handle 
systems of PDEs with arbitrary number of spatial dimensions and solution components. It 
provides the spatial discretization and time integration functions, functions to read and 
write solutions from/to files, as well as functions required to solve the system on parallel 
(MPI) platforms. The physical models define the physics-specific functions such as the exact 
forms of the hyperbolic flux, parabolic flux, source terms, upwinding functions, etc.

Documentation
-------------

See documentation on how to compile and run along with examples:
+ https://hypar.readthedocs.io/en/latest/
+ http://hypar.github.io/

Building HyPar
--------------

HyPar supports two build systems: **Autotools** (traditional) and **CMake** (modern alternative).

### Building with Autotools

The traditional build system uses GNU Autotools. Quick start:

```bash
autoreconf -i
./configure
make -j$(nproc)
```

The executable `HyPar` will be created in the `bin/` directory.

**Common configure options:**

```bash
# Serial build (without MPI)
./configure --enable-serial

# Enable CUDA support
./configure --enable-cuda --with-cuda-dir=/path/to/cuda

# Enable OpenMP
./configure --enable-omp

# Enable ScaLAPACK
./configure --enable-scalapack --with-blas-dir=/path/to/blas --with-lapack-dir=/path/to/lapack

# Enable FFTW (requires MPI)
./configure --enable-fftw --with-fftw-dir=/path/to/fftw

# Specify MPI installation
./configure --with-mpi-dir=/path/to/mpi

# Enable Python features
./configure --enable-python
```

**Installation:**
```bash
make install
```

### Building with CMake

HyPar now also supports building with CMake. For detailed instructions, see [CMAKE.md](CMAKE.md).

Quick start:
```bash
mkdir build && cd build
cmake ..
cmake --build . -j$(nproc)
```

The executable `HyPar` will be created in the `build/src/` directory.

**Common CMake options:**

```bash
# Serial build
cmake -DENABLE_SERIAL=ON ..

# Enable CUDA support
cmake -DENABLE_CUDA=ON ..

# Enable OpenMP
cmake -DENABLE_OMP=ON ..

# Enable ScaLAPACK
cmake -DENABLE_SCALAPACK=ON ..

# Enable FFTW
cmake -DENABLE_FFTW=ON ..

# Enable Python features
cmake -DENABLE_PYTHON=ON ..
```



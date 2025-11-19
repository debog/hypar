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

## Testing

HyPar includes comprehensive unit tests for core numerical functions.

### Running Tests

**With Autotools:**
```bash
make check
```

**With CMake:**
```bash
cd build
ctest
# Or for verbose output:
ctest --output-on-failure
```

### Test Coverage

The test suite includes **58 tests** across **9 modules**:
- InterpolationFunctions: 6 tests
- ArrayFunctions: 9 tests
- FirstDerivative: 4 tests
- SecondDerivative: 4 tests
- TridiagLU: 3 tests
- CommonFunctions: 4 tests
- LimiterFunctions: 5 tests
- MathFunctions: 8 tests
- MPIFunctions: 15 tests (12 parallel MPI + 3 unit tests)

**Note**: MPIFunctions tests require MPI and run with 2 processes. They are automatically disabled in serial builds.

For detailed test documentation, see [tests/README.md](tests/README.md).



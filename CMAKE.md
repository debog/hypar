# Building HyPar with CMake

This document describes how to build HyPar using CMake as an alternative to the traditional autotools-based build system.

## Quick Start

```bash
mkdir build
cd build
cmake ..
cmake --build . -j$(nproc)
```

The executable `HyPar` will be built in the `build/src/` directory.

## Running Tests

HyPar includes comprehensive unit tests that can be run with CMake/CTest:

```bash
cd build
ctest
```

Or for more verbose output:
```bash
ctest --output-on-failure
```

To disable building tests:
```bash
cmake -DBUILD_TESTING=OFF ..
```

**Note**: In serial mode (`-DENABLE_SERIAL=ON`), the MPIFunctions tests are automatically disabled. All other tests (8 modules, 46 tests) will run.

## Build Options

CMake provides several options to customize the build:

### Serial Mode
Build without MPI support:
```bash
cmake -DENABLE_SERIAL=ON ..
```

### OpenMP Support
Enable OpenMP threads:
```bash
cmake -DENABLE_OMP=ON ..
```

### CUDA Support
Enable CUDA for GPU acceleration:
```bash
cmake -DENABLE_CUDA=ON ..
```

### ScaLAPACK Support
Enable ScaLAPACK for parallel linear algebra:
```bash
cmake -DENABLE_SCALAPACK=ON ..
```

### FFTW Support
Enable FFTW for FFT operations (requires MPI, not available in serial mode):
```bash
cmake -DENABLE_FFTW=ON ..
```

### Python Support
Enable Python features (requires environment variables to be set):
```bash
cmake -DENABLE_PYTHON=ON ..
```

Required environment variables for Python support:
- `PYTHON_INCLUDE_PATH`: Path to Python include directory
- `PYTHON_LIB_PATH`: Path to Python library directory
- `PYTHON_LIB_NAME`: Python library name (e.g., `python3.11`)
- `PYTHON_BIN_PATH`: Path to Python binary directory
- `NUMPY_INCLUDE_PATH` (optional): Path to NumPy include directory

## Optional Dependencies

### PETSc Integration
If the `PETSC_DIR` and `PETSC_ARCH` environment variables are set and point to a valid PETSc installation, CMake will automatically enable PETSc support.

### libROM Integration
If the `LIBROM_DIR` environment variable is set and points to a valid libROM installation, CMake will automatically enable libROM support.

## Installation

To install HyPar to a custom location:

```bash
cmake -DCMAKE_INSTALL_PREFIX=/path/to/install ..
cmake --build .
cmake --install .
```

This will install:
- The `HyPar` executable to `<prefix>/bin/`
- Header files to `<prefix>/include/`
- Documentation to `<prefix>/share/doc/hypar/`

## Example Builds

### Basic parallel build with MPI:
```bash
mkdir build && cd build
cmake ..
cmake --build . -j8
```

### Serial build with OpenMP:
```bash
mkdir build && cd build
cmake -DENABLE_SERIAL=ON -DENABLE_OMP=ON ..
cmake --build . -j8
```

### Parallel build with CUDA and PETSc:
```bash
export PETSC_DIR=/path/to/petsc
export PETSC_ARCH=arch-linux-c-opt
mkdir build && cd build
cmake -DENABLE_CUDA=ON ..
cmake --build . -j8
```

### Build with all features:
```bash
export PETSC_DIR=/path/to/petsc
export PETSC_ARCH=arch-linux-c-opt
export LIBROM_DIR=/path/to/librom
mkdir build && cd build
cmake -DENABLE_OMP=ON -DENABLE_CUDA=ON -DENABLE_SCALAPACK=ON -DENABLE_FFTW=ON ..
cmake --build . -j8
```

## Comparison with Autotools

The CMake build system provides the same functionality as the autotools-based system:

| Autotools | CMake |
|-----------|-------|
| `./configure --enable-serial` | `cmake -DENABLE_SERIAL=ON ..` |
| `./configure --enable-cuda` | `cmake -DENABLE_CUDA=ON ..` |
| `./configure --enable-omp` | `cmake -DENABLE_OMP=ON ..` |
| `./configure --enable-scalapack` | `cmake -DENABLE_SCALAPACK=ON ..` |
| `./configure --enable-fftw` | `cmake -DENABLE_FFTW=ON ..` |
| `./configure --enable-python` | `cmake -DENABLE_PYTHON=ON ..` |
| `./configure --with-mpi-dir=<path>` | Use `CMAKE_PREFIX_PATH` or let CMake find MPI automatically |
| `./configure --with-cuda-dir=<path>` | Use `CMAKE_CUDA_COMPILER` or let CMake find CUDA automatically |
| `make` | `cmake --build .` |
| `make install` | `cmake --install .` |

## Troubleshooting

### MPI not found
If CMake cannot find MPI, you can specify the MPI compiler wrappers:
```bash
cmake -DMPI_C_COMPILER=/path/to/mpicc -DMPI_CXX_COMPILER=/path/to/mpicxx ..
```

### CUDA not found
Specify the CUDA toolkit path:
```bash
cmake -DCUDAToolkit_ROOT=/path/to/cuda ..
```

### PETSc not detected
Ensure `PETSC_DIR` and `PETSC_ARCH` environment variables are correctly set:
```bash
export PETSC_DIR=/path/to/petsc
export PETSC_ARCH=arch-linux-c-opt
cmake ..
```

### Clean rebuild
To perform a clean rebuild:
```bash
rm -rf build
mkdir build
cd build
cmake ..
cmake --build .
```

## Advanced Options

### Build type
Set the build type (Debug, Release, RelWithDebInfo, MinSizeRel):
```bash
cmake -DCMAKE_BUILD_TYPE=Release ..
```

### Compiler flags
Customize compiler flags:
```bash
cmake -DCMAKE_C_FLAGS="-O3 -march=native" -DCMAKE_CXX_FLAGS="-O3 -march=native" ..
```

### Verbose build
See the actual compiler commands:
```bash
cmake --build . --verbose
```

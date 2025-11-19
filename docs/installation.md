# Installation

HyPar is a C-based finite-difference solver that supports two build systems: **GNU Autotools** (traditional) and **CMake** (modern). It supports MPI for parallel computing and can optionally be compiled with PETSc, CUDA, and other libraries.

## Prerequisites

### Required
- C compiler (GCC recommended)
- C++ compiler (G++ recommended, for C++11 support)
- Standard math library
- **For Autotools**: GNU Autotools (autoconf, automake, libtool)
- **For CMake**: CMake 3.12 or higher

### Optional
- **MPI**: For parallel computing (OpenMPI or MPICH)
- **PETSc**: For implicit/IMEX time integration methods
- **CUDA**: For GPU acceleration
- **OpenMP**: For shared-memory parallelization
- **BLAS/LAPACK**: For linear algebra operations
- **ScaLAPACK**: For parallel linear algebra
- **FFTW**: For Fourier transforms
- **libROM**: For reduced-order modeling capabilities
- **Python**: For Python interface features

## Basic Installation

HyPar can be built using either Autotools or CMake. Choose the method that best suits your workflow.

### Option 1: Building with CMake (Recommended)

#### 1. Clone the Repository

```bash
git clone https://github.com/debog/hypar.git
cd hypar
```

#### 2. Create Build Directory

```bash
mkdir build
cd build
```

#### 3. Configure

For a basic serial build:

```bash
cmake -DENABLE_SERIAL=ON ..
```

For a basic parallel (MPI) build (CMake will auto-detect MPI):

```bash
cmake ..
```

#### 4. Compile

```bash
cmake --build . -j$(nproc)
```

The compiled executable will be created at `build/src/HyPar`.

### Option 2: Building with Autotools

#### 1. Clone the Repository

```bash
git clone https://github.com/debog/hypar.git
cd hypar
```

#### 2. Generate Configure Script

If the `configure` script is not present, generate it using:

```bash
autoreconf -i
```

#### 3. Configure

For a basic serial build:

```bash
./configure --enable-serial
```

For a basic parallel (MPI) build:

```bash
./configure
```

If MPI is installed in a non-standard location:

```bash
./configure --with-mpi-dir=/path/to/mpi
```

#### 4. Compile

```bash
make
```

The compiled executable `bin/HyPar` will be created in the source directory.

## Advanced Configuration Options

### CMake Options

#### With PETSc Support

PETSc provides implicit and IMEX time integration methods. Set environment variables before running CMake:

```bash
export PETSC_DIR=/path/to/petsc
export PETSC_ARCH=arch-linux-c-opt
cmake ..
```

CMake will automatically detect and configure PETSc if these variables are set.

#### With CUDA Support

For GPU acceleration:

```bash
cmake -DENABLE_CUDA=ON ..
```

To specify CUDA location:

```bash
cmake -DENABLE_CUDA=ON -DCUDAToolkit_ROOT=/path/to/cuda ..
```

#### With OpenMP Support

For shared-memory parallelization:

```bash
cmake -DENABLE_OMP=ON ..
```

#### With FFTW

For Fourier transform capabilities (requires MPI):

```bash
cmake -DENABLE_FFTW=ON ..
```

#### With ScaLAPACK

For parallel linear algebra:

```bash
cmake -DENABLE_SCALAPACK=ON ..
```

#### With libROM Support

Set the environment variable before running CMake:

```bash
export LIBROM_DIR=/path/to/librom
cmake ..
```

#### With Python Interface

```bash
cmake -DENABLE_PYTHON=ON ..
```

Ensure the following environment variables are set:
- `PYTHON_LIB_PATH`
- `PYTHON_BIN_PATH`
- `PYTHON_INCLUDE_PATH`
- `NUMPY_INCLUDE_PATH` (optional)
- `PYTHON_LIB_NAME`

#### CMake Complete Example

A full-featured CMake build with multiple options:

```bash
export PETSC_DIR=/opt/petsc
export PETSC_ARCH=arch-linux-c-opt
export LIBROM_DIR=/opt/librom

mkdir build && cd build
cmake \
  -DENABLE_CUDA=ON \
  -DCUDAToolkit_ROOT=/usr/local/cuda \
  -DENABLE_OMP=ON \
  -DENABLE_FFTW=ON \
  -DENABLE_SCALAPACK=ON \
  ..
cmake --build . -j$(nproc)
```

For more detailed CMake options and troubleshooting, see [CMAKE.md](https://github.com/debog/hypar/blob/master/CMAKE.md) in the repository.

### Autotools Options

#### Parallel with MPI

```bash
./configure --with-mpi-dir=/path/to/mpi
```

#### With PETSc Support

PETSc provides implicit and IMEX time integration methods:

```bash
./configure --with-petsc-dir=/path/to/petsc --with-petsc-arch=arch-name
```

Make sure the `PETSC_DIR` and `PETSC_ARCH` environment variables are set.

#### With CUDA Support

For GPU acceleration:

```bash
./configure --enable-cuda --with-cuda-dir=/path/to/cuda
```

#### With OpenMP Support

For shared-memory parallelization:

```bash
./configure --enable-omp
```

#### With FFTW

For Fourier transform capabilities:

```bash
./configure --enable-fftw --with-fftw-dir=/path/to/fftw
```

#### With ScaLAPACK

For parallel linear algebra:

```bash
./configure --enable-scalapack \
            --with-blas-dir=/path/to/blas \
            --with-lapack-dir=/path/to/lapack \
            --with-scalapack-dir=/path/to/scalapack
```

#### With Python Interface

```bash
./configure --enable-python
```

Ensure the following environment variables are set:
- `PYTHON_LIB_PATH`
- `PYTHON_BIN_PATH`
- `PYTHON_INCLUDE_PATH`
- `NUMPY_INCLUDE_PATH`
- `PYTHON_LIB_NAME`

#### Autotools Complete Example

A full-featured Autotools build with multiple options:

```bash
./configure \
  --with-mpi-dir=/opt/openmpi \
  --with-petsc-dir=/opt/petsc \
  --with-petsc-arch=arch-linux-c-opt \
  --enable-cuda \
  --with-cuda-dir=/usr/local/cuda \
  --enable-omp \
  --enable-fftw \
  --with-fftw-dir=/usr/local
make
```

## Installation

After compilation, you can optionally install the executable to a system location:

```bash
make install
```

By default, files are installed in the source directory. To change the installation prefix:

```bash
./configure --prefix=/usr/local
make
make install
```

## Verification

After compilation, test the installation by running one of the provided examples:

```bash
cd Examples/1D/LinearAdvection/SineWave
../../../../bin/HyPar
```

You should see output indicating the simulation is running successfully.

## Troubleshooting

### Configure Issues

If configure fails to find libraries, ensure:
- Library paths are correctly specified with `--with-*-dir` flags
- Environment variables like `MPI_DIR`, `PETSC_DIR`, `PETSC_ARCH` are set
- Required development packages are installed (e.g., `libopenmpi-dev`)

### Compilation Issues

- Ensure C/C++ compilers support C99 and C++11 standards
- Check that all dependent libraries are compiled with compatible compilers
- For CUDA, ensure the CUDA toolkit version is compatible with your GPU

### Runtime Issues

- For MPI builds, use the appropriate MPI launcher: `mpiexec -n 4 bin/HyPar`
- Ensure all required input files are present in the run directory
- Check that library paths are in `LD_LIBRARY_PATH` (Linux) or `DYLD_LIBRARY_PATH` (macOS)

## Building Documentation

To build this documentation locally:

```bash
cd docs
pip install -r requirements-docs.txt
make html
```

Open `_build/html/index.html` in your browser to view the documentation.
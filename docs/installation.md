# Installation

HyPar is a C-based finite-difference solver that uses GNU Autotools for configuration and compilation. It supports MPI for parallel computing and can optionally be compiled with PETSc, CUDA, and other libraries.

## Prerequisites

### Required
- C compiler (GCC recommended)
- C++ compiler (G++ recommended, for C++11 support)
- Standard math library

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

### 1. Clone the Repository

```bash
git clone https://github.com/debog/hypar.git
cd hypar
```

### 2. Generate Configure Script

If the `configure` script is not present, generate it using:

```bash
autoreconf -i
```

### 3. Configure

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

### 4. Compile

```bash
make
```

The compiled executable `bin/HyPar` will be created in the source directory.

## Advanced Configuration Options

### Parallel with MPI

```bash
./configure --with-mpi-dir=/path/to/mpi
```

### With PETSc Support

PETSc provides implicit and IMEX time integration methods:

```bash
./configure --with-petsc-dir=/path/to/petsc --with-petsc-arch=arch-name
```

Make sure the `PETSC_DIR` and `PETSC_ARCH` environment variables are set.

### With CUDA Support

For GPU acceleration:

```bash
./configure --enable-cuda --with-cuda-dir=/path/to/cuda
```

### With OpenMP Support

For shared-memory parallelization:

```bash
./configure --enable-omp
```

### With FFTW

For Fourier transform capabilities:

```bash
./configure --enable-fftw --with-fftw-dir=/path/to/fftw
```

### With ScaLAPACK

For parallel linear algebra:

```bash
./configure --enable-scalapack \
            --with-blas-dir=/path/to/blas \
            --with-lapack-dir=/path/to/lapack \
            --with-scalapack-dir=/path/to/scalapack
```

### With Python Interface

```bash
./configure --enable-python
```

Ensure the following environment variables are set:
- `PYTHON_LIB_PATH`
- `PYTHON_BIN_PATH`
- `PYTHON_INCLUDE_PATH`
- `NUMPY_INCLUDE_PATH`
- `PYTHON_LIB_NAME`

### Complete Example

A full-featured build with multiple options:

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
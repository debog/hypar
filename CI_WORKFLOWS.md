# GitHub Actions CI Workflows for CMake

This document describes the GitHub Actions CI workflows that test HyPar CMake builds.

## Workflow Files

### 1. cmake-serial.yml
**Purpose:** Tests serial (non-MPI) build
**Configuration:**
- Serial mode enabled (`-DENABLE_SERIAL=ON`)
- No MPI dependencies
- Minimal configuration

### 2. cmake-mpich.yml
**Purpose:** Tests parallel build with MPICH
**Configuration:**
- Default parallel build (MPI auto-detected)
- Uses MPICH implementation
- Tests basic MPI functionality

### 3. cmake-openmpi.yml
**Purpose:** Tests parallel build with OpenMPI
**Configuration:**
- Default parallel build (MPI auto-detected)
- Uses OpenMPI implementation
- Alternative to MPICH for MPI testing

### 4. cmake-openmpi-with-petsc.yml
**Purpose:** Tests build with PETSc integration
**Configuration:**
- Builds PETSc from source
- Sets PETSC_DIR and PETSC_ARCH environment variables
- Tests CMake's automatic PETSc detection

### 5. cmake-openmpi-with-fftw.yml
**Purpose:** Tests build with FFTW support
**Configuration:**
- Builds FFTW from source with MPI support
- Enables FFTW option (`-DENABLE_FFTW=ON`)
- Tests FFTW integration

### 6. cmake-openmpi-with-openmp.yml
**Purpose:** Tests build with OpenMP threading
**Configuration:**
- Enables OpenMP (`-DENABLE_OMP=ON`)
- Tests hybrid MPI+OpenMP build

### 7. cmake-matrix.yml
**Purpose:** Comprehensive matrix build testing
**Configuration:**
Tests multiple configurations in parallel:
- Serial build
- MPI with MPICH
- MPI with OpenMPI
- MPI + OpenMP
- MPI + ScaLAPACK

This workflow provides broad coverage of different build options.

### 8. cmake-build-types.yml
**Purpose:** Tests different CMake build types
**Configuration:**
Tests three build types:
- Debug
- Release
- RelWithDebInfo

Ensures the project builds correctly with different optimization levels.

### 9. cmake-install.yml
**Purpose:** Tests full build and installation process
**Configuration:**
- Builds with custom install prefix
- Tests `cmake --install`
- Verifies installed files (binary, headers, documentation)

### 10. cmake-mpich-with-cuda.yml
**Purpose:** Tests CUDA build with MPICH
**Configuration:**
- Installs CUDA Toolkit 11.2
- Enables CUDA support (`-DENABLE_CUDA=ON`)
- Uses GCC 10 for compatibility
- Tests CUDA compilation (runtime requires GPU hardware)

### 11. cmake-openmpi-with-cuda.yml
**Purpose:** Tests CUDA build with OpenMPI
**Configuration:**
- Installs CUDA Toolkit 11.2
- Enables CUDA support with OpenMPI
- Uses GCC 10 for compatibility
- Verifies CUDA libraries are linked

### 12. cmake-cuda-compile-test.yml
**Purpose:** Matrix test of multiple CUDA versions
**Configuration:**
- Tests CUDA 11.2 and 11.8
- Matrix build for version compatibility
- Verbose build output for debugging
- Verifies CUDA object files are created

## Workflow Triggers

All workflows are triggered on:
- `push` events to any branch
- `pull_request` events

## Dependencies

Each workflow installs necessary dependencies:
- **Base:** cmake, build-essential, g++, gfortran
- **MPI:** mpich/libmpich-dev OR openmpi-bin/libopenmpi-dev
- **PETSc:** Built from GitLab source (release branch)
- **FFTW:** Downloaded and built from source (version 3.3.10)
- **ScaLAPACK:** libscalapack-mpi-dev, libblas-dev, liblapack-dev
- **CUDA:** CUDA Toolkit 11.2 or 11.8 from NVIDIA repositories (with GCC 10)

## Build Process

Standard workflow steps:
1. Checkout code
2. Install dependencies
3. Configure CMake (`mkdir build && cmake [options] ..`)
4. Build (`cmake --build build -j$(nproc)`)
5. Verify executable exists

## Comparison with Autotools Workflows

Existing Autotools workflows:
- `c-cpp-serial.yml` → `cmake-serial.yml`
- `c-cpp-mpich.yml` → `cmake-mpich.yml`
- `c-cpp-openmpi.yml` → `cmake-openmpi.yml`
- `c-cpp-mpich-with-petsc.yml` → `cmake-openmpi-with-petsc.yml`
- `c-cpp-mpich-with-fftw.yml` → `cmake-openmpi-with-fftw.yml`
- `c-cpp-mpich-nvcc11.yml` → `cmake-mpich-with-cuda.yml` / `cmake-openmpi-with-cuda.yml`

The CMake workflows provide equivalent coverage plus additional testing (matrix builds, build types, installation, multiple CUDA versions).

## Testing Strategy

The workflows use a layered approach:
1. **Basic builds** (serial, mpich, openmpi) - Test core functionality
2. **Feature builds** (petsc, fftw, openmp) - Test optional features
3. **CUDA builds** (mpich+cuda, openmpi+cuda, cuda matrix) - Test GPU compilation
4. **Matrix build** - Comprehensive configuration coverage
5. **Build types** - Test optimization levels
6. **Installation** - Verify install process

### CUDA Testing Notes

CUDA workflows test compilation only, as GitHub Actions runners don't have GPU hardware:
- Installs CUDA Toolkit (nvcc compiler and libraries)
- Compiles all CUDA source files (`.cu` files)
- Verifies CUDA libraries are linked
- Cannot run the executable (requires actual GPU)
- Tests multiple CUDA versions for compatibility

## Adding New Workflows

To add a new CMake workflow:

1. Create `.github/workflows/cmake-<feature>.yml`
2. Use existing workflows as templates
3. Install required dependencies
4. Set appropriate CMake options
5. Test locally before committing

Example:
```yaml
name: CMake Build - <Description>

on: [push, pull_request]

jobs:
  build:
    runs-on: ubuntu-latest
    steps:
    - uses: actions/checkout@v3
    - name: Install dependencies
      run: |
        sudo apt-get update
        sudo apt-get install -y cmake <other-deps>
    - name: Configure CMake
      run: |
        mkdir build
        cd build
        cmake <options> ..
    - name: Build
      run: cmake --build build -j$(nproc)
    - name: Verify
      run: test -f build/src/HyPar
```

## Maintenance

These workflows should be updated when:
- New CMake options are added
- Dependency versions need updating
- Build process changes
- New features require testing

## Status Badges

To add status badges to README.md:

```markdown
![CMake Build Matrix](https://github.com/<owner>/hypar/workflows/CMake%20Build%20Matrix/badge.svg)
![CMake Serial](https://github.com/<owner>/hypar/workflows/CMake%20Build%20-%20Linux%20GCC%20Serial/badge.svg)
```

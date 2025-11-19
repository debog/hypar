# GitHub Actions CI Workflows for CMake

This document describes the GitHub Actions CI workflows that test HyPar CMake builds.

## Workflow Files

### 1. cmake-serial.yml
**Purpose:** Tests serial (non-MPI) build
**Configuration:**
- Serial mode enabled (`-DENABLE_SERIAL=ON`)
- No MPI dependencies
- Minimal configuration

### 2. cmake-matrix.yml
**Purpose:** Comprehensive matrix build testing
**Configuration:**
Tests multiple MPI configurations in parallel:
- MPI with MPICH
- MPI with OpenMPI
- MPI + OpenMP
- MPI + ScaLAPACK

This workflow provides broad coverage of different MPI-based build options.

### 3. cmake-openmpi-with-petsc.yml
**Purpose:** Tests build with PETSc integration
**Configuration:**
- Builds PETSc from source
- Sets PETSC_DIR and PETSC_ARCH environment variables
- Tests CMake's automatic PETSc detection

### 4. cmake-openmpi-with-fftw.yml
**Purpose:** Tests build with FFTW support
**Configuration:**
- Builds FFTW from source with MPI support
- Enables FFTW option (`-DENABLE_FFTW=ON`)
- Tests FFTW integration

### 5. cmake-mpich-with-cuda.yml
**Purpose:** Tests CUDA build with MPICH
**Configuration:**
- Installs CUDA Toolkit 11.2
- Enables CUDA support (`-DENABLE_CUDA=ON`)
- Uses GCC 10 for compatibility
- Tests CUDA compilation (runtime requires GPU hardware)
- Limited parallelism (-j2) to prevent memory issues

### 6. cmake-build-types.yml
**Purpose:** Tests different CMake build types
**Configuration:**
Tests three build types:
- Debug
- Release
- RelWithDebInfo

Ensures the project builds correctly with different optimization levels.

### 7. cmake-install.yml
**Purpose:** Tests full build and installation process
**Configuration:**
- Builds with custom install prefix
- Tests `cmake --install`
- Verifies installed files (binary, headers, documentation)

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

The workflows use an efficient layered approach:
1. **Serial build** (`cmake-serial.yml`) - Tests non-MPI compilation
2. **Matrix build** (`cmake-matrix.yml`) - Tests 4 MPI configurations in parallel
3. **Feature builds** - Test complex integrations:
   - PETSc (`cmake-openmpi-with-petsc.yml`)
   - FFTW (`cmake-openmpi-with-fftw.yml`)
   - CUDA (`cmake-mpich-with-cuda.yml`)
4. **Build types** (`cmake-build-types.yml`) - Test Debug/Release optimization
5. **Installation** (`cmake-install.yml`) - Verify install process

### CUDA Testing Notes

CUDA workflow tests compilation only, as GitHub Actions runners don't have GPU hardware:
- Installs CUDA Toolkit 11.2 (nvcc compiler and libraries)
- Compiles all CUDA source files (`.cu` files)
- Verifies CUDA libraries are linked
- Uses limited parallelism (-j2) to prevent memory issues
- Cannot run the executable (requires actual GPU)

### Workflow Optimization

To minimize CI time and resource usage:
- **Matrix builds** cover multiple configurations in one workflow
- **Dedicated workflows** only for complex features (PETSc, FFTW, CUDA)
- **Removed redundancies**: Single serial test, single CUDA test
- **Total: 7 workflows** (reduced from 12) covering all features

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

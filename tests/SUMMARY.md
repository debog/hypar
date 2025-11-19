# GitHub CI Tests for Interpolation Functions - Summary

## Overview

This document summarizes the GitHub CI test infrastructure created for the InterpolationFunctions in HyPar.

## Files Created

### 1. Test Source Code
- **`tests/InterpolationFunctions/test_interpolation.c`** (11KB)
  - Comprehensive unit test suite for interpolation functions
  - Tests 6 different interpolation schemes
  - Includes 6 test cases covering various function types
  - Simple custom test harness (no external dependencies)

### 2. Build Configuration
- **`tests/InterpolationFunctions/Makefile.am`**
  - Autotools build configuration for tests
  - Links against necessary HyPar libraries
  - Defines `test_interpolation` as a check program

- **`tests/Makefile.am`**
  - Parent Makefile for tests directory
  - Includes InterpolationFunctions subdirectory

- **Modified: `configure.ac`**
  - Added tests/Makefile and tests/InterpolationFunctions/Makefile to AC_CONFIG_FILES
  - Enables autotools to generate Makefiles for tests

- **Modified: `Makefile.am`** (root)
  - Added `tests` to SUBDIRS
  - Integrates tests into main build system

### 3. GitHub Actions Workflow
- **`.github/workflows/unit-tests.yml`** (1.9KB)
  - Automated CI workflow for unit tests
  - Runs on every push and pull request
  - Two test configurations:
    - Serial mode (--enable-serial)
    - MPI mode (default)
  - Displays test logs on failure

### 4. Documentation
- **`tests/README.md`** (3.5KB)
  - Comprehensive documentation for the test suite
  - Explains test structure and coverage
  - Instructions for running and adding tests
  - Future test additions roadmap

- **`tests/SUMMARY.md`** (this file)
  - High-level overview of created infrastructure

### 5. Git Configuration
- **`tests/.gitignore`**
  - Excludes test executables, logs, and build artifacts
  - Keeps repository clean

## Test Coverage

### Current Tests

| Interpolation Scheme | Test Cases | Tested Functions |
|---------------------|------------|------------------|
| First Order Upwind | 2 | Constant, Linear |
| Second Order Central | 2 | Linear, Quadratic |
| Fourth Order Central | 1 | Cubic |
| Fifth Order WENO | 1 | Smooth Sine Wave |

### Tested Interpolation Functions

1. `Interp1PrimFirstOrderUpwind()` - Component-wise 1st order upwind
2. `Interp1PrimSecondOrderCentral()` - Component-wise 2nd order central
3. `Interp1PrimFourthOrderCentral()` - Component-wise 4th order central
4. `Interp1PrimFifthOrderWENO()` - Component-wise 5th order WENO

### Test Properties Verified

- **Polynomial Reproduction**: Schemes correctly reproduce polynomials up to their design order
- **Constant Preservation**: All schemes preserve constant functions exactly
- **Stability**: No spurious oscillations for smooth functions
- **Accuracy**: Interface values match analytical expectations within tolerance (1e-10)

## How to Use

### Running Tests Locally

```bash
# From HyPar root directory
autoreconf -i
./configure --enable-serial  # or without for MPI
make
make check
```

### Viewing Test Results

Tests output results in a clear format:
```
[PASS] First Order Upwind - Constant Function
[PASS] Second Order Central - Linear Function
...
Test Results:
  Passed: 6
  Failed: 0
```

### CI Integration

Tests automatically run on GitHub Actions:
1. Push code or create pull request
2. GitHub Actions triggers unit-tests.yml workflow
3. Tests run in both serial and MPI configurations
4. Results visible in Actions tab
5. Failures block merging (if configured)

## Architecture

### Test Design Principles

1. **Self-contained**: No external test frameworks required
2. **Fast**: All tests complete in < 1 second
3. **Deterministic**: Same input always produces same result
4. **Isolated**: Tests don't depend on each other
5. **Clear Output**: Pass/fail clearly indicated

### Test Structure

Each test function follows this pattern:
```c
int test_scheme_property() {
  // 1. Setup test data (arrays, solver structure)
  // 2. Call interpolation function
  // 3. Verify results against expected values
  // 4. Clean up allocated memory
  // 5. Return TEST_PASS or TEST_FAIL
}
```

## Future Enhancements

### Additional Tests to Add

1. **MUSCL Schemes**
   - Second Order MUSCL (with limiters)
   - Third Order MUSCL (with Koren limiter)

2. **High-Order WENO Variants**
   - Fifth Order CRWENO (Compact Reconstruction)
   - Fifth Order HCWENO (Hybrid Compact)
   - Fifth Order Upwind

3. **Compact Schemes**
   - Fifth Order Compact Upwind

4. **Characteristic-Based**
   - All schemes with characteristic decomposition
   - Requires defining physics (e.g., Euler equations)

5. **Edge Cases**
   - Discontinuities (shock capturing)
   - Multiple variables (nvars > 1)
   - Multi-dimensional tests
   - Ghost point handling

6. **Performance Tests**
   - Timing benchmarks
   - Memory usage verification

### Integration with Existing Tests

- Current tests are **unit tests** (test individual functions)
- Complement existing **regression tests** (test complete simulations)
- Both types of tests are valuable:
  - Unit tests: Fast, pinpoint failures
  - Regression tests: Verify overall correctness

## Dependencies

### Build Dependencies
- autoconf
- automake
- libtool
- C compiler (gcc)
- make

### Runtime Dependencies
- None (tests are self-contained)

### Optional Dependencies
- MPI (for MPI mode testing)
- MPICH or OpenMPI (GitHub Actions uses MPICH)

## Maintenance

### When to Update Tests

1. **Adding new interpolation schemes**: Add corresponding test cases
2. **Modifying existing schemes**: Update expected values if behavior changes
3. **Bug fixes**: Add regression test for the bug
4. **Performance improvements**: Verify correctness not compromised

### Test Failure Debugging

If tests fail:
1. Check test logs: `tests/InterpolationFunctions/test_interpolation.log`
2. Run test manually: `./tests/InterpolationFunctions/test_interpolation`
3. Use debugger: `gdb ./tests/InterpolationFunctions/test_interpolation`
4. Check recent code changes
5. Verify build configuration

## Integration with HyPar Build System

The tests are fully integrated into the autotools build system:

```
HyPar/
├── configure.ac          ← Modified (added tests)
├── Makefile.am           ← Modified (added tests)
├── src/                  ← Existing source code
└── tests/                ← New test directory
    ├── Makefile.am       ← New
    └── InterpolationFunctions/
        ├── Makefile.am   ← New
        └── test_interpolation.c ← New
```

Running `make check` will:
1. Build all HyPar libraries
2. Build test executables
3. Run all tests
4. Report pass/fail status
5. Generate test logs

## Contact

For questions or issues with the test infrastructure:
- Open an issue on GitHub
- Check existing documentation in `tests/README.md`
- Review test source code for examples

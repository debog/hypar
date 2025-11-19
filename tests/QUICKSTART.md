# Quick Start Guide - Running Interpolation Function Tests

## One-Command Test Run

```bash
# From HyPar root directory
autoreconf -i && ./configure --enable-serial && make && make check
```

## Step-by-Step Instructions

### 1. Configure the Build System

```bash
cd /path/to/hypar
autoreconf -i
```

This generates the configure script and Makefiles.

### 2. Configure HyPar

**Option A: Serial Mode (Faster, No MPI required)**
```bash
./configure --enable-serial
```

**Option B: MPI Mode (Requires MPI installation)**
```bash
./configure
```

### 3. Build HyPar

```bash
make -j$(nproc)
```

This compiles:
- All HyPar libraries in `src/`
- Test executables in `tests/`

### 4. Run Tests

```bash
make check
```

This runs all unit tests and displays results.

## Understanding Test Output

### Successful Test Run

```
========================================
Interpolation Functions Unit Tests
========================================

Testing First Order Upwind Scheme:
[PASS] First Order Upwind - Constant Function
[PASS] First Order Upwind - Linear Function

Testing Second Order Central Scheme:
[PASS] Second Order Central - Linear Function
[PASS] Second Order Central - Quadratic Function

Testing Fourth Order Central Scheme:
[PASS] Fourth Order Central - Cubic Function

Testing WENO Scheme:
[PASS] Fifth Order WENO - Smooth Sine Wave

========================================
Test Results:
  Passed: 6
  Failed: 0
========================================
```

Exit code: 0 (success)

### Failed Test Run

```
Testing Second Order Central Scheme:
[PASS] Second Order Central - Linear Function
[FAIL] Second Order Central - Quadratic Function
    At interface 3: got 9.500000, expected 9.250000
...
Test Results:
  Passed: 5
  Failed: 1
========================================

Last failed test: Second Order Central - Quadratic Function
```

Exit code: 1 (failure)

## Running Individual Tests

### Run Specific Test Executable

```bash
cd tests/InterpolationFunctions
./test_interpolation
```

### Debug with GDB

```bash
cd tests/InterpolationFunctions
gdb ./test_interpolation
(gdb) run
(gdb) bt  # if it crashes
```

### Check Test Logs

```bash
cat tests/InterpolationFunctions/test-suite.log
cat tests/InterpolationFunctions/test_interpolation.log
```

## Common Issues and Solutions

### Issue 1: `autoreconf: command not found`

**Solution:** Install autotools
```bash
# Ubuntu/Debian
sudo apt-get install autoconf automake libtool

# macOS
brew install autoconf automake libtool
```

### Issue 2: `configure: error: cannot find MPI`

**Solution:** Either install MPI or use serial mode
```bash
# Option A: Install MPI
sudo apt-get install mpich libmpich-dev

# Option B: Use serial mode
./configure --enable-serial
```

### Issue 3: `make: *** No rule to make target 'check'`

**Cause:** Makefiles not generated properly

**Solution:** Reconfigure
```bash
make distclean
autoreconf -i
./configure --enable-serial
make
make check
```

### Issue 4: Tests fail with "Segmentation fault"

**Cause:** Array bounds issue or null pointer

**Solution:** Run with debugger
```bash
cd tests/InterpolationFunctions
gdb ./test_interpolation
(gdb) run
(gdb) backtrace
```

### Issue 5: `undefined reference` linker errors

**Cause:** Missing library dependencies

**Solution:** Check Makefile.am has all required libraries
```makefile
test_interpolation_LDADD = \
  ../../src/InterpolationFunctions/libInterpolationFunctions.a \
  ../../src/ArrayFunctions/libArrayFunctions.a \
  ../../src/MathFunctions/libMathFunctions.a \
  ../../src/CommonFunctions/libCommonFunctions.a
```

## GitHub Actions CI

Tests automatically run on GitHub for:
- Every push to any branch
- Every pull request

### Viewing CI Results

1. Go to repository on GitHub
2. Click "Actions" tab
3. Select latest workflow run
4. View "Unit Tests (Serial)" and "Unit Tests (MPI)" jobs
5. Check test output in job logs

### CI Configuration

File: `.github/workflows/unit-tests.yml`

The workflow:
1. Checks out code
2. Installs dependencies
3. Configures HyPar (both serial and MPI)
4. Builds HyPar
5. Runs `make check`
6. Displays test logs on failure

## Clean Build

To start fresh:

```bash
# Clean build artifacts
make clean

# Clean everything including configure
make distclean

# Or manually remove generated files
rm -rf autom4te.cache config.h config.log config.status Makefile
find . -name Makefile.in -delete
find . -name Makefile -delete

# Then reconfigure
autoreconf -i
./configure --enable-serial
make
make check
```

## Test Development Workflow

### Adding a New Test

1. Edit `tests/InterpolationFunctions/test_interpolation.c`
2. Add test function:
```c
int test_my_new_feature() {
  // Setup
  // Execute
  // Verify
  return TEST_PASS or TEST_FAIL;
}
```
3. Register in `main()`:
```c
print_test_result(&stats, "My New Feature Test", test_my_new_feature());
```
4. Rebuild and test:
```bash
cd tests/InterpolationFunctions
make clean
make
./test_interpolation
```

### Testing Your Changes

Before committing:
```bash
# Run tests locally
make check

# Verify no memory leaks (if valgrind available)
cd tests/InterpolationFunctions
valgrind --leak-check=full ./test_interpolation

# Push and check CI results on GitHub
git push
```

## Performance Tips

### Faster Builds

```bash
# Use parallel compilation
make -j$(nproc)

# Only rebuild tests (after source change)
cd tests
make clean
make
```

### Faster Test Runs

Tests are already very fast (< 1 second). To run only interpolation tests:

```bash
cd tests/InterpolationFunctions
make check
```

## Getting Help

1. **Test Documentation**: `tests/README.md`
2. **Test Summary**: `tests/SUMMARY.md`
3. **GitHub Issues**: Open an issue for bugs or questions
4. **Code Comments**: Read `test_interpolation.c` for examples

## Next Steps

- Read `tests/README.md` for detailed documentation
- Read `tests/SUMMARY.md` for architecture overview
- Add tests for additional interpolation schemes
- Integrate with your development workflow

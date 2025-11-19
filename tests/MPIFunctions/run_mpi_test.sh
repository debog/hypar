#!/bin/bash
# Wrapper script to run MPI test with 2 ranks

set -e  # Exit on error

# Check if test executable exists
if [ ! -f ./test_mpi_parallel ]; then
    echo "Error: test_mpi_parallel executable not found in $(pwd)"
    exit 1
fi

# Try to find MPI executable in common locations
MPI_CMD=""
MPI_FLAGS=""

# First try standard commands
if command -v mpirun &> /dev/null; then
    MPI_CMD="mpirun"
    # Check if this is OpenMPI (supports --oversubscribe)
    if mpirun --version 2>&1 | grep -iq "open mpi"; then
        MPI_FLAGS="--oversubscribe"
    fi
elif command -v mpiexec &> /dev/null; then
    MPI_CMD="mpiexec"
# Try common installation paths
elif [ -x /usr/bin/mpirun ]; then
    MPI_CMD="/usr/bin/mpirun"
elif [ -x /usr/bin/mpiexec ]; then
    MPI_CMD="/usr/bin/mpiexec"
elif [ -x /usr/local/bin/mpirun ]; then
    MPI_CMD="/usr/local/bin/mpirun"
elif [ -x /usr/local/bin/mpiexec ]; then
    MPI_CMD="/usr/local/bin/mpiexec"
else
    echo "Error: MPI launcher not found"
    echo "Searched for: mpirun, mpiexec in PATH and common locations"
    echo "PATH: $PATH"
    echo "Available commands:"
    which -a mpirun mpiexec 2>/dev/null || echo "  None found"
    exit 1
fi

echo "Using MPI command: $MPI_CMD"
if [ -n "$MPI_FLAGS" ]; then
    echo "MPI flags: $MPI_FLAGS"
fi

# Run the MPI tests with 2 ranks
$MPI_CMD $MPI_FLAGS -n 2 ./test_mpi_parallel
exit $?

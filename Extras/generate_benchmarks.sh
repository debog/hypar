#!/bin/bash

#####################################################################
#
# README:
#
# This script generates benchmark solutions for HyPar. This is available
# at: https://github.com/debog/hypar_benchmarks_mpi
# These are generated using the master branch of the main HyPar repo.
# In general, these benchmarks should be used. If you know what you
# are doing, modify this script to generate a new set of benchmarks,
# but **do not** try to commit them to the repo.
#
# ** Please copy this script to a new location before using it! **
# It will create new copies of HyPar and benchmarks solutions at
# the location it is run from.
#
#####################################################################

clear
root_dir=$PWD

# Details about HyPar (repo and branch that we want to
# generate benchmarks from)
hypar_repo="https://github.com/debog/hypar.git"
hypar_branch="master"
# other HyPar-related stuff
hypar_dir="hypar"
hypar_exec="HyPar"
hypar_compile_log_file="compile_hypar.log"

#export env vars for other scripts to run HyPar
export HYPAR_EXEC_W_PATH="${root_dir}/${hypar_dir}/bin/${hypar_exec}"
export MPI_EXEC="mpiexec"
export HYPAR_EXEC_OTHER_ARGS=""

# HyPar benchmarks repo to update
hypar_benchmarks_repo="git@github.com:debog/hypar_benchmarks_mpi.git"
hypar_benchmarks_branch="master"
hypar_benchmarks_dir="benchmarks_mpi"

# other stuff
RUN_SCRIPT="run.sh"
NEEDS_PETSC="dep.PETSc"
NEEDS_FFTW="dep.fftw"
DISABLED=".disabled"

# Clone HyPar
echo "-------------------------"
echo "removing any previous HyPar copy..."
rm -rf $hypar_dir
echo "-------------------------"
echo "cloning HyPar from $hypar_repo ($hypar_branch branch)..."
git clone $hypar_repo $hypar_dir
cd $root_dir/$hypar_dir
git checkout $hypar_branch

# compile hypar
echo "-------------------------"
echo "compiling hypar..."
autoreconf -i 2>&1 > $root_dir/$hypar_compile_log_file
if [ -z "$PETSC_DIR" ]; then
  echo "Environment variable PETSC_DIR not found."
  echo "Compiling without PETSc; will not be able to run simulations that need PETSc."
  opt_with_petsc=false
else
  echo "PETSc found at $PETSC_DIR. Compiling with PETSc."
  opt_with_petsc=true
fi
if [ -z "$LIBROM_DIR" ]; then
  echo "Environment variable LIBROM_DIR not found."
  echo "Compiling without libROM; will not be able to run simulations that need libROM."
  opt_with_librom=false
else
  echo "libROM found at $LIBROM_DIR. Compiling with libROM."
  opt_with_librom=true
fi
if [ -z "$FFTW_DIR" ]; then
  echo "Environment variable FFTW_DIR not found."
  echo "Compiling without FFTW; will not be able to run simulations that need FFTW."
  opt_with_fftw=false
  ./configure 2>&1 >> $root_dir/$hypar_compile_log_file
else
  echo "FFTW found at $FFTW_DIR."
  echo "Compiling with FFTW."
  opt_with_fftw=true
  ./configure --enable-fftw --with-fftw-dir=${FFTW_DIR} 2>&1 >> $root_dir/$hypar_compile_log_file
fi
make -j all 2>&1 >> $root_dir/$hypar_compile_log_file
make install 2>&1 >> $root_dir/$hypar_compile_log_file

# done with compiling hypar
cd $root_dir
echo "finished compilation. see $hypar_compile_log_file for log."

if [ -f "$HYPAR_EXEC_W_PATH" ]; then

  echo "compilation successful."
  echo "-------------------------"

else

  echo "compilation failed!"
  echo "-------------------------"

fi

# clone benchmarks
if [ -d "$hypar_benchmarks_dir" ]; then
  cd $hypar_benchmarks_dir
  if [ -d ".git" ]; then
    echo "benchmarks directory already exists. refreshing it."
    git reset HEAD --hard
    git checkout $hypar_benchmarks_branch
    git pull --force
    cd ../
  else
    echo "benchmarks directory exists, but not a git repo. erasing..."
    cd $root_dir
    rm -rf $hypar_benchmarks_dir
    echo "getting HyPar benchmarks from $hypar_benchmarks_repo ($hypar_benchmarks_branch branch)"
    git clone $hypar_benchmarks_repo $hypar_benchmarks_dir
    cd $hypar_benchmarks_dir
    git checkout $hypar_benchmarks_branch
  fi
else
  echo "getting HyPar benchmarks from $hypar_benchmarks_repo ($hypar_benchmarks_branch branch)"
  git clone $hypar_benchmarks_repo $hypar_benchmarks_dir
  cd $hypar_benchmarks_dir
  git checkout $hypar_benchmarks_branch
fi
cd $root_dir
echo "-------------------------"

# run the baseline cases
cd $root_dir/$hypar_benchmarks_dir

echo "generating benchmarks ..."
echo " "
for f in *; do
  if [ -d "$f" ]; then
    echo "entering $f..."
    cd $f
    if [ -f "$NEEDS_FFTW" ] && [ "$opt_with_fftw" == "false" ]; then
      echo "Skipping; $f has unmet dependencies (FFTW)."
    elif [ -f "$NEEDS_PETSC" ] && [ "$opt_with_petsc" == "false" ]; then
      echo "Skipping; $f has unmet dependencies (PETSc)."
    elif [ -f "$DISABLED" ]; then
      echo "Skipping; $f test is disabled."
    else
      if [ -f "$RUN_SCRIPT" ]; then
        chmod +x $RUN_SCRIPT && ./$RUN_SCRIPT
      else
        echo "Error: $RUN_SCRIPT not found."
      fi
    fi
    echo " "
    cd ../
  fi
done
echo "done."
echo "-------------------------"
cd $root_dir

# check git status for benchmarks repo
cd $root_dir/$hypar_benchmarks_dir
echo "checking git status..."
git status
echo "if necessary, manually commit changes."
echo "-------------------------"


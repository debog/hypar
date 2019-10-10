#!/bin/bash

#####################################################################
#
# README:
#
# This script generates baseline solutions for HyPar. A set of baseline
# solutions is available at: 
# https://deboghosh@bitbucket.org/deboghosh/hypar_baselines.git
# These are generated using the master branch of the main HyPar repo.
# In general, these baselines should be used. If you know what you
# are doing, modify this script to generate a new set of baselines,
# but **do not** try to commit them to the main baselines repo.
#
# ** Please copy this script to a new location before using it! **
# It will create new copies of HyPar and baseline solutions at
# the location it is run from.
#
#####################################################################

clear
root_dir=$PWD

# Details about HyPar (repo and branch that we want to
# generate baselines from)
hypar_repo="https://deboghosh@bitbucket.org/deboghosh/hypar.git"
hypar_branch="master"
# other HyPar-related stuff
hypar_dir="hypar"
hypar_exec="HyPar"
hypar_compile_log_file="compile_hypar.log"

#export env vars for other scripts to run HyPar
export HYPAR_EXEC_W_PATH="${root_dir}/${hypar_dir}/bin/${hypar_exec}"
export MPI_EXEC="mpiexec"
export HYPAR_EXEC_OTHER_ARGS=""

# HyPar baselines repo to update
hypar_baselines_repo="git@bitbucket.org:deboghosh/hypar_baselines.git"
hypar_baselines_branch="master"
hypar_baselines_dir="baselines"

# other stuff
RUN_SCRIPT="run.sh"

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
./configure 2>&1 >> $root_dir/$hypar_compile_log_file
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

# clone baselines
if [ -d "$hypar_baselines_dir" ]; then
  cd $hypar_baselines_dir
  if [ -d ".git" ]; then
    echo "baselines directory already exists. refreshing it."
    git reset HEAD --hard
    git checkout $hypar_baselines_branch
    git pull --force
    cd ../
  else
    echo "baselines directory exists, but not a git repo. erasing..."
    cd $root_dir
    rm -rf $hypar_baselines_dir
    echo "getting HyPar baselines from $hypar_baselines_repo ($hypar_baselines_branch branch)"
    git clone $hypar_baselines_repo $hypar_baselines_dir
    cd $hypar_baselines_dir
    git checkout $hypar_baselines_branch
  fi
else
  echo "getting HyPar baselines from $hypar_baselines_repo ($hypar_baselines_branch branch)"
  git clone $hypar_baselines_repo $hypar_baselines_dir
  cd $hypar_baselines_dir
  git checkout $hypar_baselines_branch
fi
cd $root_dir
echo "-------------------------"

# run the baseline cases
cd $root_dir/$hypar_baselines_dir

echo "generating baselines ..."
echo " "
for f in *; do
  if [ -d "$f" ]; then
    echo "entering $f..."
    cd $f
    if [ -f "$RUN_SCRIPT" ]; then
      chmod +x $RUN_SCRIPT && ./$RUN_SCRIPT
    else
      echo "Error: $RUN_SCRIPT not found."
    fi
    echo " "
    cd ../
  fi
done
echo "done."
echo "-------------------------"
cd $root_dir

# check git status for baselines repo
cd $root_dir/$hypar_baselines_dir
echo "checking git status..."
git status
echo "if necessary, manually commit changes."
echo "-------------------------"


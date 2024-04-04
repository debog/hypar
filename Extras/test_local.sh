#!/bin/bash

#################################################################
#
# README:
#
# This script tests a local copy of HyPar (edit the variable
# hypar_dir to the location of the code you want to test) against
# a set of baselines available as:
# https://gitlab.com/debojyoti.ghosh/hypar_baselines
# The variables hypar_baselines_repo and hypar_baselines_branch can
# be changed below to test against a different set of baselines.
#
# ** Please copy this script to a new location before using it! **
# It will create new copies of HyPar and baseline solutions and
# create test directories at the location it is run from.
#
#################################################################

clear
root_dir=$PWD

# HyPar location
hypar_dir="/path/to/hypar"
# other HyPar-related stuff
hypar_exec="HyPar"

#export env vars for other scripts to run HyPar
export HYPAR_EXEC_W_PATH="${hypar_dir}/bin/${hypar_exec}"
export MPI_EXEC="mpiexec"
export HYPAR_EXEC_OTHER_ARGS=""

# some details about HyPar baselines
# (baseline solutions maintained on the public repository)
# do not change these, unless you know what you are doing
hypar_baselines_repo="https://gitlab.com/debojyoti.ghosh/hypar_baselines.git"
hypar_baselines_branch="master"
hypar_baselines_dir="baselines"

# stuff about test directory
hypar_test_dir="_test"
exclude_flag="--exclude={'op*','surface*','ibank*','initial*','out.log','README.md','.git*'}"
diff_filelistname="diff_file_list"
diff_file="diff.log"
report_filename="test_report.txt"

# other stuff
RUN_SCRIPT="run.sh"
NEEDS_PETSC="dep.PETSc"
NEEDS_FFTW="dep.fftw"
NEEDS_LIBROM="dep.libROM"
DISABLED=".disabled"

if [ -f "$HYPAR_EXEC_W_PATH" ]; then

  echo "HyPar binary found."
  echo "-------------------------"

else

  echo "---------------------------------"
  echo "ERROR !!!"
  echo " "
  echo "HyPar binary NOT FOUND !!!"
  echo " "
  echo "$HYPAR_EXEC_W_PATH does not exist"
  echo " "
  echo "---------------------------------"

fi

# Since this is using a pre-compiled binary, we can only
# guess whether it was compiled using various dependencies
# by checking for the environment variables
if [ -z "$PETSC_DIR" ]; then
  echo "Environment variable PETSC_DIR not found."
  echo "Will skip tests that need PETSc"
  opt_with_petsc=false
else
  echo "PETSc found at $PETSC_DIR."
  echo "HyPar was probably compiled with PETSc. Will attempt tests that need PETSc"
  opt_with_petsc=true
fi
if [ -z "$LIBROM_DIR" ]; then
  echo "Environment variable LIBROM_DIR not found."
  echo "Will skip tests that need libROM"
  opt_with_librom=false
else
  echo "libROM found at $LIBROM_DIR."
  echo "HyPar was probably compiled with libROM. Will attempt tests that need libROM"
  opt_with_librom=true
fi
if [ -z "$FFTW_DIR" ]; then
  echo "Environment variable FFTW_DIR not found."
  echo "Will skip tests that need FFTW"
  opt_with_fftw=false
else
  echo "FFTW found at $FFTW_DIR."
  echo "HyPar was probably compiled with FFTW. Will attempt tests that need FFTW"
  opt_with_fftw=true
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

# create test dir and copy input files
timestamp=`date | sed -e 's/ /_/g' -e 's/:/./g'`
test_dirname=${hypar_test_dir}_${timestamp}
rm -rf $test_dirname && mkdir $test_dirname
echo "copying test cases to $test_dirname ..."
rsync_cmd="rsync -a $exclude_flag $root_dir/$hypar_baselines_dir/ $root_dir/$test_dirname/"
eval $rsync_cmd

# run the cases
cd $root_dir/$test_dirname
report_file="${root_dir}/${test_dirname}/${report_filename}"

rm -rf $report_file && touch $report_file
echo "HyPar Tests" >> $report_file
echo "Date/Time       : $(date '+%d/%m/%Y %H:%M:%S')" >> $report_file
echo "HyPar repo      : $hypar_repo" >> $report_file
echo "HyPar branch    : $hypar_branch" >> $report_file
echo "Baselines repo  : $hypar_baselines_repo" >> $report_file
echo "Baselines branch: $hypar_baselines_branch" >> $report_file
echo " " >> $report_file
echo " " >> $report_file

n_pass=0
n_fail=0
n_skip=0
echo "running tests..."
echo " "
for f in *; do
  if [ -d "$f" ]; then
    echo "entering $f..."
    echo "$f" >> $report_file
    cd $f
    if [ -f "$NEEDS_LIBROM" ] && [ "$opt_with_librom" == "false" ]; then
      echo "Skipping; $f has unmet dependencies (liROM)."
      echo "Skipping; $f has unmet dependencies (liROM)." >> $report_file
      ((n_skip+=1))
    elif [ -f "$NEEDS_PETSC" ] && [ "$opt_with_petsc" == "false" ]; then
      echo "Skipping; $f has unmet dependencies (PETSc)."
      echo "Skipping; $f has unmet dependencies (PETSc)." >> $report_file
      ((n_skip+=1))
    elif [ -f "$NEEDS_FFTW" ] && [ "$opt_with_fftw" == "false" ]; then
      echo "Skipping; $f has unmet dependencies (FFTW)."
      echo "Skipping; $f has unmet dependencies (FFTW)." >> $report_file
      ((n_skip+=1))
    elif [ -f "$DISABLED" ]; then
      echo "Skipping; $f test is disabled."
      echo "Skipping; $f test is disabled." >> $report_file
      ((n_skip+=1))
    else
      if [ -f "$RUN_SCRIPT" ]; then
        chmod +x $RUN_SCRIPT && ./$RUN_SCRIPT
        while read F  ; do
          echo "    comparing $F ..."
          echo "  $F" >> $report_file
          result=$(diff $F $root_dir/$hypar_baselines_dir/$f/$F 2>&1 >> $diff_file)
          if [ -z "$result" ]; then
            if [ -s "$diff_file" ]; then
              ((n_fail+=1))
              echo "                        **DIFFERENCES FOUND**"
              echo "                        **DIFFERENCES FOUND**" >> $report_file
              echo "        check"
              echo "          $root_dir/$test_dirname/$f/$diff_file"
              echo "        check" >> $report_file
              echo "          $root_dir/$test_dirname/$f/$diff_file" >> $report_file
              echo " "
            else
              ((n_pass+=1))
              echo "                        **passed**"
              echo "                        **passed**" >> $report_file
            fi
          else
            ((n_fail+=1))
            echo $result >> $report_file
            echo "                        **FILE COMPARISON FAILED**"
            echo "                        **FILE COMPARISON FAILED**" >> $report_file
            echo "        check"
            echo "          $report_file"
            echo " "
          fi
        done <./$diff_filelistname
      else
        echo "Error: $RUN_SCRIPT not found."
        echo "Error: $RUN_SCRIPT not found." >> $report_file
      fi
    fi
    echo " "
    echo " " >> $report_file
    cd ../
  fi
done
echo "done."
echo "-------------------------"
cd $root_dir

echo "all done. Bye!"
echo "report in ${report_file}"
echo "$n_pass file comparisons **passed**"
echo "$n_fail file comparisons **failed**"
echo "$n_skip tests **skipped**"
echo "-------------------------"

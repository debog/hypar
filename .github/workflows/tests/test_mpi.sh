#!/bin/bash

root_dir=$PWD

# HyPar location (HYPAR_DIR should exist in the environment)
hypar_dir=$HYPAR_DIR
# other HyPar-related stuff
hypar_exec="HyPar"

#export env vars for other scripts to run HyPar
export HYPAR_EXEC_W_PATH="${hypar_dir}/bin/${hypar_exec}"
export MPI_EXEC="mpiexec"
export HYPAR_EXEC_OTHER_ARGS=""

# some details about HyPar benchmarks
# (benchmark solutions maintained on the public repository)
# do not change these, unless you know what you are doing
hypar_benchmarks_repo="https://github.com/debog/hypar_benchmarks_mpi.git"
hypar_benchmarks_branch="master"
hypar_benchmarks_dir="benchmarks"

# stuff about test directory
hypar_test_dir="_test"
exclude_flag="--exclude={'op*','surface*','ibank*','initial*','out.log','README.md','.git*'}"
diff_filelistname="diff_file_list"
diff_file="diff.log"

# other stuff
RUN_SCRIPT="run.sh"
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
if [ -z "$FFTW_DIR" ]; then
  echo "Environment variable FFTW_DIR not found."
  echo "Will skip tests that need FFTW"
  opt_with_fftw=false
else
  echo "FFTW found at $FFTW_DIR."
  echo "HyPar was probably compiled with FFTW. Will attempt tests that need FFTW"
  opt_with_fftw=true
fi

# compile the HyPar diff code
hypar_diff_srcname="Extras/hyparDiff_RegTests.c"
HYPAR_DIFF="HYPAR_DIFF"
if [ -f "$hypar_dir/$hypar_diff_srcname" ]; then
  echo "Compiling HyPar-diff."
  gcc $hypar_dir/$hypar_diff_srcname -lm -o $root_dir/$HYPAR_DIFF
else
  echo "---------------------------------"
  echo "ERROR !!!"
  echo " "
  echo "HyPar-Diff source NOT FOUND !!!"
  echo " "
  echo "$hypar_dir/$hypar_diff_srcname does not exist"
  echo " "
  echo "---------------------------------"
fi
HYPAR_DIFF_CMD="$root_dir/$HYPAR_DIFF -r 1.0e-14 "

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

# create test dir and copy input files
timestamp=`date | sed -e 's/ /_/g' -e 's/:/./g'`
test_dirname=${hypar_test_dir}_${timestamp}
rm -rf $test_dirname && mkdir $test_dirname
echo "copying test cases to $test_dirname ..."
rsync_cmd="rsync -a $exclude_flag $root_dir/$hypar_benchmarks_dir/ $root_dir/$test_dirname/"
eval $rsync_cmd

# run the cases
cd $root_dir/$test_dirname

echo "HyPar Tests"
echo "Date/Time       : $(date '+%d/%m/%Y %H:%M:%S')"
echo "HyPar repo      : $hypar_repo"
echo "HyPar branch    : $hypar_branch"
echo "benchmarks repo  : $hypar_benchmarks_repo"
echo "benchmarks branch: $hypar_benchmarks_branch"
echo " "
echo " "

n_pass=0
n_fail=0
n_skip=0
echo "running tests..."
echo " "
for f in *; do
  if [ -d "$f" ]; then
    echo "entering $f..."
    cd $f
    if [ -f "$DISABLED" ]; then
      echo "Skipping; $f test is disabled."
      ((n_skip+=1))
    else
      if [ -f "$RUN_SCRIPT" ]; then
        chmod +x $RUN_SCRIPT && ./$RUN_SCRIPT
        while read F  ; do
          echo "    comparing $F ..."
          result=$($HYPAR_DIFF_CMD $F $root_dir/$hypar_benchmarks_dir/$f/$F 2>&1 >> $diff_file)
          if [ -z "$result" ]; then
            if [ -s "$diff_file" ]; then
              ((n_fail+=1))
              echo "                        **DIFFERENCES FOUND**"
              echo "### Dump of HyPar-diff output ####"
              cat $diff_file
              echo "### End ####"
              echo " "
            else
              ((n_pass+=1))
              echo "                        **passed**"
            fi
          else
            ((n_fail+=1))
            echo "                        **FILE COMPARISON FAILED**"
            echo "Command output: $result "
            echo " "
            echo "### Dump of screen output   ####"
            cat out.log
            echo "### End ####"
            echo "### Directory contents   ####"
            ls -lh ./
            echo "### End ####"
            echo "### Benchmark directory contents   ####"
            ls -lh $root_dir/$hypar_benchmarks_dir/$f/
            echo "### End ####"
          fi
        done <./$diff_filelistname
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

echo "all done. Bye!"
echo "$n_pass file comparisons **passed**"
if [[ $n_fail -gt 0 ]]; then
  echo "$n_fail file comparisons **failed**"
  exit 1
fi
echo "-------------------------"

rm -rf $root_dir/$HYPAR_DIFF

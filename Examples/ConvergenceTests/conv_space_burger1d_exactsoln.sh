#!/bin/bash

#################################################################
#
# Bash script to run a spatial convergence test on a smooth 1D Burgers
# equation problem - evolution of a sine wave.
#
# ** Please copy this script to a new location before using it! **
# It will create test directories at the location it is run from.
#
# The exact solution is available in analytical form for the
# example being simulated.
#################################################################

clear

# Please adjust the variables below as needed.

# HyPar location
hypar_dir="/g/g92/ghosh5/Codes/hypar"
# HyPar binary name
hypar_exec="HyPar"
# therefore, HyPar binary with path
HYPAR_EXEC_W_PATH="${hypar_dir}/bin/${hypar_exec}"
# write HyPar screen output to this file
hypar_stdout="out.log"
# add path to HyPar's bash scripts to path
export PATH=$PATH:$hypar_dir/Examples/BashTools

if [ -f "$HYPAR_EXEC_W_PATH" ]; then

  echo "-------------------------"
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
  exit

fi
echo ""

# command to launch MPI jobs (mpiexec, srun, etc.)
mpicmd="mpiexec"
# other necessary args for launching MPI jobs,
# eg., queue/allocation specifications
mpiargs=""

# Convergence test parameters
num_pts_coarsest=16
num_refinements=5
nproc_coarsest=1

# set default parameters as env vars
source set_defaults.sh
# overwrite env vars as needed for this test case
export p_hyp_scheme="weno5"
export p_ts="rk"
export p_tstype="44"
export p_dt="0.0001"
export p_niter=10000
export p_model="burgers"

# path to initial solution code for this case
physics_name="Burgers"
case_name="SineWave"
init_c_file="${hypar_dir}/Examples/${p_ndims}D/${physics_name}/${case_name}/aux/init.c"
init_compile_cmd="gcc -lm"
init_exec="INIT"

# create run directories with names starting with this string
RUNDIRNAME_ROOT=run_grid_

# file to write errors in
ERRFILE="convergence.dat"

##################################################################
# DO NOT EDIT BELOW THIS LINE UNLESS YOU KNOW WHAT YOU ARE DOING #
##################################################################

############################################
# Function for floating point calculations #
calc(){ awk "BEGIN { print "$*" }"; }      #
#                                          #
############################################

root_dir=$PWD

# Clean up
rm -rf $ERRFILE ${RUNDIRNAME_ROOT}*

HERE=$PWD

REF_LEVEL=0
MAX_NUM_PTS=$num_pts_coarsest
while [ $REF_LEVEL -lt $num_refinements ]; do
  let REF_LEVEL=REF_LEVEL+1
  MAX_NUM_PTS=$(calc $MAX_NUM_PTS*2)
done
MAX_NUM_PTS=$(calc $MAX_NUM_PTS/2)
DIRINDEXLEN=${#MAX_NUM_PTS}

echo "--"
echo "Spatial convergence test on a smooth solution to the Burgers "
echo "equation (1D):"
echo "Grid sizes varying from $num_pts_coarsest to $MAX_NUM_PTS:"
echo "Spatial discretization (hyperbolic): ${p_hyp_scheme}"
echo "Spatial discretization (parabolic): ${p_par_scheme}"
echo "Time integration: ${p_ts} (${p_tstype})"
echo "--"
echo ""

REF_LEVEL=0
NUM_PTS=$num_pts_coarsest
NPROC=$nproc_coarsest
while [ $REF_LEVEL -lt $num_refinements ]; do

  echo "Running: ref_lvl=$REF_LEVEL, n=$NUM_PTS, dt=$p_dt, maxiter=$p_niter, $NPROC cores ..."

  DIRINDEX=$NUM_PTS
  len=${#DIRINDEX}
  while [ "$len" -lt $DIRINDEXLEN ]; do
    DIRINDEX=0$DIRINDEX
    len=${#DIRINDEX}
  done
  DIRNAME=${RUNDIRNAME_ROOT}${DIRINDEX}

  echo "  creating run directory ${DIRNAME}"
  mkdir $DIRNAME
  echo "  entering run directory ${DIRNAME}"
  cd $HERE/$DIRNAME

  # set grid size and number of MPI ranks
  export p_size=$NUM_PTS
  export p_iproc=$NPROC

  # write input files
  echo "    writing input files"
  write_solver_inp.sh
  write_boundary_inp.sh
  write_physics_inp_burgers.sh

  # compile the code to generate initial solution and run it
  echo "    compiling and running code to generate initial solution"
  eval '${init_compile_cmd} $init_c_file -o $init_exec'
  eval './${init_exec} 2>&1 > init.log'
  rm -rf $init_exec

  # run HyPar
  RUNCOMMAND="${mpicmd} -n ${NPROC} ${mpiargs} ${HYPAR_EXEC_W_PATH} 2>&1 > $hypar_stdout"
  START=$(date +%s.%N)
  echo "    running HyPar ($RUNCOMMAND)"
  eval $RUNCOMMAND
  END=$(date +%s.%N)
  RUNTIME=$(echo "$END - $START" | bc)
  echo "    wallclock time = $RUNTIME seconds"

  echo "  exiting run directory ${DIRNAME}"
  cd $HERE
  echo ""

  let REF_LEVEL=REF_LEVEL+1
  NUM_PTS=$(calc $NUM_PTS*2)
  NPROC=$(calc $NPROC*2)

done

# gather errors and wall times into one text file
cat ${RUNDIRNAME_ROOT}*/errors.dat > $ERRFILE
echo ""

#!/bin/bash

#################################################################
#
# Bash script to run a spatial convergence test on a smooth 2D linear
# advection problem - periodic advection of a sine wave.
#
# ** Please copy this script to a new location before using it! **
# It will create test directories at the location it is run from.
#
# The exact solution is not available in analytical form, so the
# solution of a finer grid is used as the "exact" solution to
# compute error norms.
#
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
num_refinements=4
nproc_coarsest=1

# set default parameters as env vars
source set_defaults.sh
# overwrite env vars as needed for this test case
export p_ndims=2
export p_hyp_scheme="upw5"
export p_ts="rk"
export p_tstype="44"
export p_dt="0.0025"
export p_niter=200000
export p_op_format="binary"
export p_screen_op_iter=$p_niter
export p_file_op_iter=$p_niter
export p_op_overwrite="yes"
export p_adv_file="advection"

# path to initial solution code for this case
physics_name="LinearAdvection"
case_name="SineWave_NonConstantAdvection/case_00"
init_c_file="${hypar_dir}/Examples/${p_ndims}D/${physics_name}/${case_name}/aux/init.c"
init_compile_cmd="gcc -lm"
init_exec="INIT"

# set the boundaries
export p_nb=4
BC_COUNTER=0
while [ $BC_COUNTER -lt $p_nb ]; do
  varname="p_bctype_${BC_COUNTER}"
  eval export "$varname"="periodic"
  let BC_COUNTER=BC_COUNTER+1
done
BC_COUNTER=0
while [ $BC_COUNTER -lt 2 ]; do
  varname="p_dim_${BC_COUNTER}"
  eval export "$varname"="0"
  varname="p_limits_${BC_COUNTER}"
  eval 'export "$varname"="0 0 0 1"'
  let BC_COUNTER=BC_COUNTER+1
done
BC_COUNTER=2
while [ $BC_COUNTER -lt $p_nb ]; do
  varname="p_dim_${BC_COUNTER}"
  eval export "$varname"="1"
  varname="p_limits_${BC_COUNTER}"
  eval 'export "$varname"="0 1 0 0"'
  let BC_COUNTER=BC_COUNTER+1
done
BC_COUNTER=0
while [ $BC_COUNTER -lt $p_nb ]; do
  varname="p_face_${BC_COUNTER}"
  eval export "$varname"="1"
  let BC_COUNTER=BC_COUNTER+2
done
BC_COUNTER=1
while [ $BC_COUNTER -lt $p_nb ]; do
  varname="p_face_${BC_COUNTER}"
  eval export "$varname"="-1"
  let BC_COUNTER=BC_COUNTER+2
done


# create run directories with names starting with this string
RUNDIRNAME_ROOT=run_grid_
REFERENCEDIR_NAME=reference

# file to write errors in
ERRFILE="convergence.dat"

# other stuff
output_filename="op.bin"
refsoln_filename="solution.inp"
binop2inp_compile_cmd="gcc"
binop2inp_c_file="${hypar_dir}/Extras/BinaryOPToInitialSolution.c"
binop2inp_exec="BINOP2INP"
binop2inp_inp="binop2inp.inp"

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
rm -rf $ERRFILE ${RUNDIRNAME_ROOT}* $REFERENCEDIR_NAME

HERE=$PWD

REF_LEVEL=0
MAX_NUM_PTS=$num_pts_coarsest
MAX_NPROC=$nproc_coarsest
while [ $REF_LEVEL -lt $num_refinements ]; do
  let REF_LEVEL=REF_LEVEL+1
  MAX_NUM_PTS=$(calc $MAX_NUM_PTS*2)
  MAX_NPROC=$(calc $MAX_NPROC*2)
done
MAX_NUM_PTS=$(calc $MAX_NUM_PTS/2)
MAX_NPROC=$(calc $MAX_NPROC/2)
DIRINDEXLEN=${#MAX_NUM_PTS}

NUM_PTS_REFERENCE=$(calc $MAX_NUM_PTS*2)
NPROC_REFERENCE=$MAX_NPROC

echo "--"
echo "Spatial convergence test on a smooth solution to the linear "
echo "advection equation (2D):"
echo "Spatial discretization (hyperbolic): ${p_hyp_scheme}"
echo "Spatial discretization (parabolic): ${p_par_scheme}"
echo "Time integration: ${p_ts} (${p_tstype})"
echo "--"
echo ""

### Generate reference solution ####

NUM_PTS=$NUM_PTS_REFERENCE
NPROC=$NPROC_REFERENCE
nproc_total=$(calc $NPROC*$NPROC)

echo ""
echo "Generating reference solution."
echo "Running: n=$NUM_PTS^${p_ndims}, dt=$p_dt, maxiter=$p_niter, $nproc_total cores ..."

DIRNAME=$REFERENCEDIR_NAME

echo "  creating run directory ${DIRNAME}"
mkdir $DIRNAME
echo "  entering run directory ${DIRNAME}"
cd $HERE/$DIRNAME
REFERENCE_SOLUTION_DIR=$PWD

# set grid size and number of MPI ranks
export p_size="$NUM_PTS $NUM_PTS"
export p_iproc="$NPROC $NPROC"

# write input files
echo "    writing input files"
write_solver_inp.sh
write_boundary_inp.sh
write_physics_inp_linadv.sh

# compile the code to generate initial solution and run it
echo "    compiling and running code to generate initial solution"
eval '${init_compile_cmd} $init_c_file -o $init_exec'
eval './${init_exec} 2>&1 > init.log'
rm -rf $init_exec

# run HyPar
RUNCOMMAND="${mpicmd} -n ${nproc_total} ${mpiargs} ${HYPAR_EXEC_W_PATH} 2>&1 > $hypar_stdout"
START=$(date +%s.%N)
echo "    running HyPar ($RUNCOMMAND)"
eval $RUNCOMMAND
END=$(date +%s.%N)
RUNTIME=$(echo "$END - $START" | bc)
echo "    wallclock time = $RUNTIME seconds"

# Convert output file to input format
echo "    compiling and running code to convert binary output to input format"
eval '${binop2inp_compile_cmd} $binop2inp_c_file -o $binop2inp_exec'
/bin/cat <<EOM > $binop2inp_inp
${output_filename}
EOM
eval './${binop2inp_exec} < ${binop2inp_inp} 2>&1 > binop2inp.log'
rm -rf $binop2inp_exec $binop2inp_inp *.dat

echo "  exiting run directory ${DIRNAME}"
cd $HERE
echo ""

REFERENCE_SOLUTION=${REFERENCE_SOLUTION_DIR}/${refsoln_filename}
if [[ -f "$REFERENCE_SOLUTION" ]]; then
  echo "Reference solution generated"
  echo ""
else
  echo "---------------------------------"
  echo "ERROR !!!"
  echo "Unable to generate reference solution."
  echo "---------------------------------"
  exit
fi

##### Now run the convergence test ####

echo "--"
echo "Running convergence test."
echo "Grid sizes varying from ${num_pts_coarsest}^${p_ndims} to ${MAX_NUM_PTS}^${p_ndims}:"
echo "--"
echo ""

REF_LEVEL=0
NUM_PTS=$num_pts_coarsest
NPROC=$nproc_coarsest
while [ $REF_LEVEL -lt $num_refinements ]; do

  nproc_total=$(calc $NPROC*$NPROC)
  echo "Running: ref_lvl=$REF_LEVEL, n=${NUM_PTS}^${p_ndims}, dt=$p_dt, maxiter=$p_niter, $nproc_total cores ..."

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
  export p_size="$NUM_PTS $NUM_PTS"
  export p_size_exact="$NUM_PTS_REFERENCE $NUM_PTS_REFERENCE"
  export p_iproc="$NPROC $NPROC"

  # write input files
  echo "    writing input files"
  write_solver_inp.sh
  write_boundary_inp.sh
  write_physics_inp_linadv.sh

  # compile the code to generate initial solution and run it
  echo "    compiling and running code to generate initial solution"
  eval '${init_compile_cmd} $init_c_file -o $init_exec'
  eval './${init_exec} 2>&1 > init.log'
  rm -rf $init_exec

  # set exact solution to point to reference solution
  ln -sf $REFERENCE_SOLUTION exact.inp

  # run HyPar
  RUNCOMMAND="${mpicmd} -n ${nproc_total} ${mpiargs} ${HYPAR_EXEC_W_PATH} 2>&1 > $hypar_stdout"
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

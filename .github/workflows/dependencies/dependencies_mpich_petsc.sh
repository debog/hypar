#!/usr/bin/env bash

set -eu -o pipefail

sudo apt-get update

sudo apt-get install -y --no-install-recommends\
    build-essential \
    g++ gfortran mpich libmpich-dev

git clone -b release https://gitlab.com/petsc/petsc.git ../petsc
pushd ../petsc
export PETSC_DIR=$PWD
export PETSC_ARCH=arch-linux-gcc
./configure --with-cc=mpicc --with-fc=mpif90 --with-cxx=mpicxx --with-shared-libraries --with-debugging=0
make -j4 all
popd

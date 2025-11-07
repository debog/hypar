#!/usr/bin/env bash

set -eu -o pipefail

sudo apt-get update

sudo apt-get install -y --no-install-recommends\
    build-essential \
    g++ gfortran openmpi-bin libopenmpi-dev \
    libblas-dev liblapack-dev

for i in {1..5}; do
  git clone --depth=1 -b release https://gitlab.com/petsc/petsc.git ../petsc && break
  echo "Clone failed, retrying ($i/5)..."
  sleep 10
done
pushd ../petsc
export PETSC_DIR=$PWD
export PETSC_ARCH=arch-linux-gcc
./configure --with-cc=mpicc --with-fc=mpif90 --with-cxx=mpicxx --with-shared-libraries --with-debugging=0
make -j4 all
popd

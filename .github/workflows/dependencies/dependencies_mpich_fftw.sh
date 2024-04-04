#!/usr/bin/env bash

set -eu -o pipefail

sudo apt-get update

sudo apt-get install -y --no-install-recommends\
    build-essential \
    g++ gfortran mpich libmpich-dev

pushd ../
wget https://www.fftw.org/fftw-3.3.10.tar.gz
tar xvof fftw-3.3.10.tar.gz
rm -rf fftw-3.3.10.tar.gz
cd fftw-3.3.10
CC=mpicc MPICC=mpicc ./configure --enable-mpi --disable-fortran --prefix=$PWD
make -j4
make install
popd

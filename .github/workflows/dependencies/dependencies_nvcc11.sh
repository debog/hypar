#!/usr/bin/env bash

set -eu -o pipefail

sudo apt-get -qqq update

sudo apt-get install -y --no-install-recommends\
    build-essential \
    ca-certificates \
    cmake

sudo apt-get install -y --no-install-recommends gcc-10 g++10
sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-10 1
sudo update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-10 1

sudo apt-get install -y --no-install-recommends\
    mpich libmpich-dev \
    gnupg \
    pkg-config \
    wget

sudo apt-key adv --fetch-keys https://developer.download.nvidia.com/compute/cuda/repos/ubuntu1804/x86_64/3bf863cc.pub
echo "deb https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2004/x86_64 /" \
    | sudo tee /etc/apt/sources.list.d/cuda.list
sudo apt-get update
sudo apt-get install -y \
    cuda-command-line-tools-11-2 \
    cuda-compiler-11-2           \
    cuda-cupti-dev-11-2          \
    cuda-minimal-build-11-2      \
    cuda-nvml-dev-11-2           \
    cuda-nvtx-11-2               \
    libcurand-dev-11-2
sudo ln -s cuda-11.2 /usr/local/cuda

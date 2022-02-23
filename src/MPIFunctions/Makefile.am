noinst_LIBRARIES = libMPIFunctions.a
libMPIFunctions_a_SOURCES = \
  MPIBroadcast.c \
  MPICommunicators.c \
  MPIExchangeBoundaries1D.c \
  MPIExchangeBoundariesnD.c \
  MPIGatherArray1D.c \
  MPIGatherArraynD.c \
  MPIGetArrayDatanD.c \
  MPIGetFilename.c \
	MPIIOGroups.c \
  MPILocalDomainLimits.c \
  MPIMax.c \
  MPIMin.c \
  MPIPartition1D.c \
  MPIPartitionArray1D.c \
  MPIPartitionArraynD.c \
  MPIRank1D.c \
  MPIRanknD.c \
  MPISum.c

if ENABLE_CUDA
noinst_LIBRARIES += libMPIFunctions_GPU.a
libMPIFunctions_GPU_a_SOURCES = \
  MPIExchangeBoundariesnD_GPU.cu
libMPIFunctions_GPU_a_LIBADD = \
  MPIExchangeBoundariesnD_GPU.dlink.o
libMPIFunctions_GPU_a_CFLAGS = $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) $(CUDA_CFLAGS)
libMPIFunctions_GPU_a_AR = $(NVCC) -arch=sm_70 -lib -o

.cu.o:
	$(NVCC) $(libMPIFunctions_GPU_a_CFLAGS) -arch=sm_70 -dc -o $@ $<
	$(NVCC) $(libMPIFunctions_GPU_a_CFLAGS) -arch=sm_70 -dlink -o $*.dlink.o $@
endif

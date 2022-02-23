noinst_LIBRARIES = libArrayFunctions.a
libArrayFunctions_a_SOURCES = \
  ArrayImplementations.c

if ENABLE_CUDA
noinst_LIBRARIES += libArrayFunctions_GPU.a
libArrayFunctions_GPU_a_SOURCES = \
  ArrayImplementations_GPU.cu
libArrayFunctions_GPU_a_LIBADD = \
  ArrayImplementations_GPU.dlink.o
libArrayFunctions_GPU_a_CFLAGS = $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) $(CUDA_CFLAGS)
libArrayFunctions_GPU_a_AR = $(NVCC) -arch=sm_70 -lib -o

.cu.o:
	$(NVCC) $(libArrayFunctions_GPU_a_CFLAGS) -arch=sm_70 -dc -o $@ $<
	$(NVCC) $(libArrayFunctions_GPU_a_CFLAGS) -arch=sm_70 -dlink -o $*.dlink.o $@
endif

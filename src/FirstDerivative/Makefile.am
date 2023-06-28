noinst_LIBRARIES = libFirstDerivative.a
libFirstDerivative_a_SOURCES = \
  FirstDerivativeFirstOrder.c \
  FirstDerivativeSecondOrder.c \
  FirstDerivativeSecondOrderNoGhosts.c \
  FirstDerivativeFourthOrder.c

if ENABLE_CUDA
noinst_LIBRARIES += libFirstDerivative_GPU.a
libFirstDerivative_GPU_a_SOURCES = \
  FirstDerivativeFourthOrder_GPU.cu
libFirstDerivative_GPU_a_LIBADD = \
  FirstDerivativeFourthOrder_GPU.dlink.o
libFirstDerivative_GPU_a_CFLAGS = $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) $(CUDA_CFLAGS)
libFirstDerivative_GPU_a_AR = $(NVCC) -arch=sm_70 -lib -o

.cu.o:
	$(NVCC) $(libFirstDerivative_GPU_a_CFLAGS) -arch=sm_70 -dc -o $@ $<
	$(NVCC) $(libFirstDerivative_GPU_a_CFLAGS) -arch=sm_70 -dlink -o $*.dlink.o $@
endif

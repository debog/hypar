noinst_LIBRARIES = libHyParFunctions.a
libHyParFunctions_a_SOURCES = \
  ApplyBoundaryConditions.c \
  ApplyIBConditions.c \
  BoundaryIntegral.c \
  CalculateError.c \
  CalculateConservationError.c \
  CalculateROMDiff.c \
	ComputeRHSOperators.c \
  ExactSolution.c \
  HyperbolicFunction.c \
  IncrementFilename.c \
  NonLinearInterpolation.c \
  ParabolicFunctionCons1Stage.c \
  ParabolicFunctionNC1Stage.c \
  ParabolicFunctionNC2Stage.c \
  ParabolicFunctionNC1.5Stage.c \
  SourceFunction.c \
  VolumeIntegral.c

if ENABLE_CUDA
noinst_LIBRARIES += libHyParFunctions_GPU.a
libHyParFunctions_GPU_a_SOURCES = \
  HyperbolicFunction_GPU.cu \
  ParabolicFunctionNC2Stage_GPU.cu
libHyParFunctions_GPU_a_LIBADD = \
  HyperbolicFunction_GPU.dlink.o \
  ParabolicFunctionNC2Stage_GPU.dlink.o
libHyParFunctions_GPU_a_CFLAGS = $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) $(CUDA_CFLAGS)
libHyParFunctions_GPU_a_AR = $(NVCC) -arch=sm_70 -lib -o

.cu.o:
	$(NVCC) $(libHyParFunctions_GPU_a_CFLAGS) -arch=sm_70 -dc -o $@ $<
	$(NVCC) $(libHyParFunctions_GPU_a_CFLAGS) -arch=sm_70 -dlink -o $*.dlink.o $@
endif

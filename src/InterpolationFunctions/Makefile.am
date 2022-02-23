noinst_LIBRARIES = libInterpolationFunctions.a
libInterpolationFunctions_a_SOURCES = \
  Interp1PrimFifthOrderCRWENO.c \
  Interp1PrimFifthOrderCRWENOChar.c \
  Interp1PrimFifthOrderCompactUpwind.c \
  Interp1PrimFifthOrderCompactUpwindChar.c \
  Interp1PrimFifthOrderHCWENO.c \
  Interp1PrimFifthOrderHCWENOChar.c \
  Interp1PrimFifthOrderUpwind.c \
  Interp1PrimFifthOrderUpwindChar.c \
  Interp1PrimFifthOrderWENO.c \
  Interp1PrimFifthOrderWENOChar.c \
  Interp1PrimFirstOrderUpwind.c \
  Interp1PrimFirstOrderUpwindChar.c \
  Interp1PrimFourthOrderCentral.c \
  Interp1PrimFourthOrderCentralChar.c \
  Interp1PrimSecondOrderMUSCL.c \
  Interp1PrimSecondOrderMUSCLChar.c \
  Interp1PrimThirdOrderMUSCL.c \
  Interp1PrimThirdOrderMUSCLChar.c \
  Interp1PrimSecondOrderCentral.c \
  Interp1PrimSecondOrderCentralChar.c \
  Interp2PrimSecondOrder.c \
  CompactSchemeCleanup.c \
  CompactSchemeInitialize.c \
  MUSCLInitialize.c \
  WENOCleanup.c \
  WENOFifthOrderCalculateWeights.c \
  WENOFifthOrderInitializeWeights.c \
  WENOInitialize.c

if ENABLE_CUDA
noinst_LIBRARIES += libInterpolationFunctions_GPU.a
libInterpolationFunctions_GPU_a_SOURCES = \
  WENOFifthOrderCalculateWeights_GPU.cu \
  Interp1PrimFifthOrderWENO_GPU.cu
libInterpolationFunctions_GPU_a_LIBADD = \
  WENOFifthOrderCalculateWeights_GPU.dlink.o \
  Interp1PrimFifthOrderWENO_GPU.dlink.o
libInterpolationFunctions_GPU_a_CFLAGS = $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) $(CUDA_CFLAGS)
libInterpolationFunctions_GPU_a_AR = $(NVCC) -arch=sm_70 -lib -o

.cu.o:
	$(NVCC) $(libInterpolationFunctions_GPU_a_CFLAGS) -arch=sm_70 -dc -o $@ $<
	$(NVCC) $(libInterpolationFunctions_GPU_a_CFLAGS) -arch=sm_70 -dlink -o $*.dlink.o $@
endif

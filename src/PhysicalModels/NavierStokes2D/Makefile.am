noinst_LIBRARIES = libNavierStokes2D.a
libNavierStokes2D_a_SOURCES = \
  NavierStokes2DCleanup.c \
  NavierStokes2DComputeCFL.c \
  NavierStokes2DEigen.c \
  NavierStokes2DFlux.c \
  NavierStokes2DFunctions.c \
	NavierStokes2DGravityField.c \
  NavierStokes2DInitialize.c \
  NavierStokes2DJacobian.c \
	NavierStokes2DModifiedSolution.c \
	NavierStokes2DParabolicFunction.c \
	NavierStokes2DPreStep.c \
	NavierStokes2DSource.c \
  NavierStokes2DUpwind.c

if ENABLE_CUDA
libNavierStokes2D_a_SOURCES += \
  NavierStokes2DInitialize_GPU.c
noinst_LIBRARIES += libNavierStokes2D_GPU.a
libNavierStokes2D_GPU_a_SOURCES = \
  NavierStokes2DFlux_GPU.cu \
  NavierStokes2DModifiedSolution_GPU.cu \
  NavierStokes2DParabolicFunction_GPU.cu \
  NavierStokes2DPreStep_GPU.cu \
  NavierStokes2DSource_GPU.cu \
  NavierStokes2DUpwind_GPU.cu
libNavierStokes2D_GPU_a_LIBADD = \
  NavierStokes2DFlux_GPU.dlink.o \
  NavierStokes2DModifiedSolution_GPU.dlink.o \
  NavierStokes2DParabolicFunction_GPU.dlink.o \
  NavierStokes2DPreStep_GPU.dlink.o \
  NavierStokes2DSource_GPU.dlink.o \
  NavierStokes2DUpwind_GPU.dlink.o
libNavierStokes2D_GPU_a_CFLAGS = $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) $(CUDA_CFLAGS)
libNavierStokes2D_GPU_a_AR = $(NVCC) -arch=sm_70 -lib -o

.cu.o:
	$(NVCC) $(libNavierStokes2D_GPU_a_CFLAGS) -arch=sm_70 -dc -o $@ $<
	$(NVCC) $(libNavierStokes2D_GPU_a_CFLAGS) -arch=sm_70 -dlink -o $*.dlink.o $@
endif

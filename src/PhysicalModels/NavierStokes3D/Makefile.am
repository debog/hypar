noinst_LIBRARIES = libNavierStokes3D.a
libNavierStokes3D_a_SOURCES = \
  NavierStokes3DCleanup.c \
  NavierStokes3DComputeCFL.c \
  NavierStokes3DEigen.c \
  NavierStokes3DFlux.c \
  NavierStokes3DFunctions.c \
  NavierStokes3DGravityField.c \
  NavierStokes3DIBForces.c \
  NavierStokes3DImmersedBoundary.c \
  NavierStokes3DInitialize.c \
	NavierStokes3DJacobian.c \
  NavierStokes3DModifiedSolution.c \
  NavierStokes3DParabolicFunction.c \
	NavierStokes3DPreStep.c \
  NavierStokes3DSource.c \
  NavierStokes3DUpwind.c

if ENABLE_CUDA
libNavierStokes3D_a_SOURCES += \
  NavierStokes3DInitialize_GPU.c \
  NavierStokes3DCleanup_GPU.c
noinst_LIBRARIES += libNavierStokes3D_GPU.a
libNavierStokes3D_GPU_a_SOURCES = \
  NavierStokes3DFlux_GPU.cu \
  NavierStokes3DModifiedSolution_GPU.cu \
  NavierStokes3DParabolicFunction_GPU.cu \
  NavierStokes3DPreStep_GPU.cu \
  NavierStokes3DSource_GPU.cu \
  NavierStokes3DUpwind_GPU.cu
libNavierStokes3D_GPU_a_LIBADD = \
  NavierStokes3DFlux_GPU.dlink.o \
  NavierStokes3DModifiedSolution_GPU.dlink.o \
  NavierStokes3DParabolicFunction_GPU.dlink.o \
  NavierStokes3DPreStep_GPU.dlink.o \
  NavierStokes3DSource_GPU.dlink.o \
  NavierStokes3DUpwind_GPU.dlink.o
libNavierStokes3D_GPU_a_CFLAGS = $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) $(CUDA_CFLAGS)
libNavierStokes3D_GPU_a_AR = $(NVCC) -arch=sm_70 -lib -o

.cu.o:
	$(NVCC) $(libNavierStokes3D_GPU_a_CFLAGS) -arch=sm_70 -dc -o $@ $<
	$(NVCC) $(libNavierStokes3D_GPU_a_CFLAGS) -arch=sm_70 -dlink -o $*.dlink.o $@
endif

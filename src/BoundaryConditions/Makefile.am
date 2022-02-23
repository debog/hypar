noinst_LIBRARIES = libBoundaryConditions.a
libBoundaryConditions_a_SOURCES = \
  BCCleanup.c \
  BCDirichlet.c \
  BCExtrapolate.c \
  BCInitialize.c \
  BCIO.c \
  BCNoFlux.c \
  BCNoslipWall.c \
  BCPeriodic.c \
  BCReflect.c \
  BCSlipWall.c \
  BCSponge.c \
	BCSubsonicAmbivalent.c \
	BCSubsonicInflow.c \
	BCSubsonicOutflow.c \
	BCSupersonicInflow.c \
	BCSupersonicOutflow.c \
  BCSWSlipWall.c \
  BCThermalNoslipWall.c \
  BCThermalSlipWall.c \
	BCTurbulentSupersonicInflow.c

if ENABLE_CUDA
noinst_LIBRARIES += libBoundaryConditions_GPU.a
libBoundaryConditions_GPU_a_SOURCES = \
  BCSlipWall_GPU.cu \
  BCPeriodic_GPU.cu
libBoundaryConditions_GPU_a_LIBADD = \
  BCSlipWall_GPU.dlink.o \
  BCPeriodic_GPU.dlink.o
libBoundaryConditions_GPU_a_CFLAGS = $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) $(CUDA_CFLAGS)
libBoundaryConditions_GPU_a_AR = $(NVCC) -arch=sm_70 -lib -o

.cu.o:
	$(NVCC) $(libBoundaryConditions_GPU_a_CFLAGS) -arch=sm_70 -dc -o $@ $<
	$(NVCC) $(libBoundaryConditions_GPU_a_CFLAGS) -arch=sm_70 -dlink -o $*.dlink.o $@
endif

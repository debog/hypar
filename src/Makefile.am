SUBDIRS = \
  ArrayFunctions \
  BandedMatrix \
  BoundaryConditions \
  CommonFunctions \
  FirstDerivative \
  HyParFunctions \
  ImmersedBoundaries \
  InterpolationFunctions \
  IOFunctions \
  LimiterFunctions \
  PETScFunctions \
  PhysicalModels \
  PlottingFunctions \
  MathFunctions \
  MPIFunctions \
	ROM \
  SecondDerivative \
  Simulation \
  SparseGrids \
  TimeIntegration \
  TridiagLU

bin_PROGRAMS = HyPar
HyPar_SOURCES = main.cpp
HyPar_LDADD = \
  SparseGrids/libSparseGrids.a \
  Simulation/libSimulation.a \
  HyParFunctions/libHyParFunctions.a \
  TimeIntegration/libTimeIntegration.a \
  PETScFunctions/libPETScFunctions.a \
  ROM/libROM.a \
  BoundaryConditions/libBoundaryConditions.a \
  PhysicalModels/LinearADR/libLinearADR.a \
  PhysicalModels/Burgers/libBurgers.a \
  PhysicalModels/FPDoubleWell/libFPDoubleWell.a \
  PhysicalModels/FPPowerSystem/libFPPowerSystem.a \
  PhysicalModels/FPPowerSystem1Bus/libFPPowerSystem1Bus.a \
  PhysicalModels/FPPowerSystem3Bus/libFPPowerSystem3Bus.a \
  PhysicalModels/Euler1D/libEuler1D.a \
  PhysicalModels/Euler2D/libEuler2D.a \
  PhysicalModels/NavierStokes2D/libNavierStokes2D.a \
  PhysicalModels/NavierStokes3D/libNavierStokes3D.a \
  PhysicalModels/Numa2D/libNuma2D.a \
  PhysicalModels/Numa3D/libNuma3D.a \
  PhysicalModels/ShallowWater1D/libShallowWater1D.a \
  PhysicalModels/ShallowWater2D/libShallowWater2D.a \
  PhysicalModels/Vlasov/libVlasov.a \
  ImmersedBoundaries/libImmersedBoundaries.a \
  InterpolationFunctions/libInterpolationFunctions.a \
  LimiterFunctions/libLimiterFunctions.a \
  FirstDerivative/libFirstDerivative.a \
  SecondDerivative/libSecondDerivative.a \
  TridiagLU/libTridiagLU.a \
  PlottingFunctions/libPlottingFunctions.a \
  IOFunctions/libIOFunctions.a \
  MPIFunctions/libMPIFunctions.a \
  BandedMatrix/libBandedMatrix.a \
  ArrayFunctions/libArrayFunctions.a \
  MathFunctions/libMathFunctions.a \
  CommonFunctions/libCommonFunctions.a

if ENABLE_CUDA
HyPar_LDADD += HyParFunctions/libHyParFunctions_GPU.a
HyPar_LDADD += BoundaryConditions/libBoundaryConditions_GPU.a
HyPar_LDADD += PhysicalModels/NavierStokes2D/libNavierStokes2D_GPU.a
HyPar_LDADD += PhysicalModels/NavierStokes3D/libNavierStokes3D_GPU.a
HyPar_LDADD += InterpolationFunctions/libInterpolationFunctions_GPU.a
HyPar_LDADD += FirstDerivative/libFirstDerivative_GPU.a
HyPar_LDADD += MPIFunctions/libMPIFunctions_GPU.a
HyPar_LDADD += ArrayFunctions/libArrayFunctions_GPU.a
endif

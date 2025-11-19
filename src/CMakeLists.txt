function(add_hypar_library lib_name source_dir)
    file(GLOB_RECURSE c_sources "${source_dir}/*.c")
    file(GLOB_RECURSE cpp_sources "${source_dir}/*.cpp")
    list(FILTER c_sources EXCLUDE REGEX ".*_GPU\\.c$")
    set(sources ${c_sources} ${cpp_sources})
    if(sources)
        add_library(${lib_name} STATIC ${sources})
        target_include_directories(${lib_name} PRIVATE ${CMAKE_SOURCE_DIR}/include)
    endif()
endfunction()

function(add_hypar_library_with_cuda lib_name source_dir)
    file(GLOB_RECURSE c_sources "${source_dir}/*.c")
    file(GLOB_RECURSE cpp_sources "${source_dir}/*.cpp")
    list(FILTER c_sources EXCLUDE REGEX ".*_GPU\\.c$")
    set(sources ${c_sources} ${cpp_sources})
    if(sources)
        add_library(${lib_name} STATIC ${sources})
        target_include_directories(${lib_name} PRIVATE ${CMAKE_SOURCE_DIR}/include)
    endif()
    
    if(ENABLE_CUDA)
        file(GLOB_RECURSE cuda_c_sources "${source_dir}/*_GPU.c")
        if(cuda_c_sources)
            target_sources(${lib_name} PRIVATE ${cuda_c_sources})
        endif()
    endif()
endfunction()

function(add_hypar_cuda_library lib_name source_dir)
    if(ENABLE_CUDA)
        file(GLOB_RECURSE cuda_sources "${source_dir}/*.cu")
        if(cuda_sources)
            add_library(${lib_name} STATIC ${cuda_sources})
            target_include_directories(${lib_name} PRIVATE ${CMAKE_SOURCE_DIR}/include)
            set_target_properties(${lib_name} PROPERTIES
                CUDA_SEPARABLE_COMPILATION ON
                CUDA_ARCHITECTURES "70"
            )
        endif()
    endif()
endfunction()

add_hypar_library(CommonFunctions CommonFunctions)
add_hypar_library(MathFunctions MathFunctions)
add_hypar_library(ArrayFunctions ArrayFunctions)
add_hypar_library(BandedMatrix BandedMatrix)
add_hypar_library(MPIFunctions MPIFunctions)
add_hypar_library(IOFunctions IOFunctions)
add_hypar_library(PlottingFunctions PlottingFunctions)
add_hypar_library(TridiagLU TridiagLU)
add_hypar_library(SecondDerivative SecondDerivative)
add_hypar_library(FirstDerivative FirstDerivative)
add_hypar_library(LimiterFunctions LimiterFunctions)
add_hypar_library(InterpolationFunctions InterpolationFunctions)
add_hypar_library(ImmersedBoundaries ImmersedBoundaries)
add_hypar_library(BoundaryConditions BoundaryConditions)
add_hypar_library(PETScFunctions PETScFunctions)
add_hypar_library(ROM ROM)
add_hypar_library(TimeIntegration TimeIntegration)
add_hypar_library(HyParFunctions HyParFunctions)
add_hypar_library(Simulation Simulation)
add_hypar_library(SparseGrids SparseGrids)

add_hypar_library(LinearADR PhysicalModels/LinearADR)
add_hypar_library(Burgers PhysicalModels/Burgers)
add_hypar_library(FPDoubleWell PhysicalModels/FPDoubleWell)
add_hypar_library(FPPowerSystem PhysicalModels/FPPowerSystem)
add_hypar_library(FPPowerSystem1Bus PhysicalModels/FPPowerSystem1Bus)
add_hypar_library(FPPowerSystem3Bus PhysicalModels/FPPowerSystem3Bus)
add_hypar_library(Euler1D PhysicalModels/Euler1D)
add_hypar_library(Euler2D PhysicalModels/Euler2D)
add_hypar_library_with_cuda(NavierStokes2D PhysicalModels/NavierStokes2D)
add_hypar_library_with_cuda(NavierStokes3D PhysicalModels/NavierStokes3D)
add_hypar_library(Numa2D PhysicalModels/Numa2D)
add_hypar_library(Numa3D PhysicalModels/Numa3D)
add_hypar_library(ShallowWater1D PhysicalModels/ShallowWater1D)
add_hypar_library(ShallowWater2D PhysicalModels/ShallowWater2D)
add_hypar_library(Vlasov PhysicalModels/Vlasov)

if(ENABLE_CUDA)
    add_hypar_cuda_library(HyParFunctions_GPU HyParFunctions)
    add_hypar_cuda_library(BoundaryConditions_GPU BoundaryConditions)
    add_hypar_cuda_library(NavierStokes2D_GPU PhysicalModels/NavierStokes2D)
    add_hypar_cuda_library(NavierStokes3D_GPU PhysicalModels/NavierStokes3D)
    add_hypar_cuda_library(InterpolationFunctions_GPU InterpolationFunctions)
    add_hypar_cuda_library(FirstDerivative_GPU FirstDerivative)
    add_hypar_cuda_library(MPIFunctions_GPU MPIFunctions)
    add_hypar_cuda_library(ArrayFunctions_GPU ArrayFunctions)
endif()

add_executable(HyPar main.cpp)

target_link_libraries(HyPar
    SparseGrids
    Simulation
    HyParFunctions
    TimeIntegration
    PETScFunctions
    ROM
    BoundaryConditions
    LinearADR
    Burgers
    FPDoubleWell
    FPPowerSystem
    FPPowerSystem1Bus
    FPPowerSystem3Bus
    Euler1D
    Euler2D
    NavierStokes2D
    NavierStokes3D
    Numa2D
    Numa3D
    ShallowWater1D
    ShallowWater2D
    Vlasov
    ImmersedBoundaries
    InterpolationFunctions
    LimiterFunctions
    FirstDerivative
    SecondDerivative
    TridiagLU
    PlottingFunctions
    IOFunctions
    MPIFunctions
    BandedMatrix
    ArrayFunctions
    MathFunctions
    CommonFunctions
)

if(ENABLE_CUDA)
    target_link_libraries(HyPar
        HyParFunctions_GPU
        BoundaryConditions_GPU
        NavierStokes2D_GPU
        NavierStokes3D_GPU
        InterpolationFunctions_GPU
        FirstDerivative_GPU
        MPIFunctions_GPU
        ArrayFunctions_GPU
    )
endif()

install(TARGETS HyPar DESTINATION bin)

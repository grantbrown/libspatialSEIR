cmake_minimum_required(VERSION 2.8.7)


project(libSpatialSEIR)
include(ExternalProject)
ExternalProject_Add(clBLAS
                    SOURCE_DIR "${CMAKE_SOURCE_DIR}/clBLAS/src"
                    BUILD_COMMAND make
                    UPDATE_COMMAND ""
                    INSTALL_DIR lib
                    )

install (FILES "${CMAKE_BINARY_DIR}/clBLAS-prefix/src/clBLAS-build/package/lib64/libclBLAS.so"
               "${CMAKE_BINARY_DIR}/clBLAS-prefix/src/clBLAS-build/package/lib64/libclBLAS.so.2"
               "${CMAKE_BINARY_DIR}/clBLAS-prefix/src/clBLAS-build/package/lib64/libclBLAS.so.2.1.0"
               DESTINATION lib)


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic -Wextra")
set(Boost_USE_STATIC_LIBS OFF) 
set(Boost_USE_MULTITHREADED ON)  
set(Boost_USE_STATIC_RUNTIME OFF) 
find_package(Boost 1.49.0 COMPONENTS random) 

set(SpatialSEIR_Boost_LIBRARIES   )
set(SpatialSEIR_Boost_INCLUDE_DIRS   )
if(Boost_FOUND)
	set(SpatialSEIR_Boost_LIBRARIES ${Boost_LIBRARIES})
	set(SpatialSEIR_Boost_INCLUDE_DIRS ${Boost_INCLUDE_DIRS})
	message(${SpatialSEIR_Boost_LIBRARIES})
	message(${SpatialSEIR_Boost_INCLUDE_DIRS})
else()
    message("Error: Boost was not found on the system. This dependency 
             will be relaxed in the release version of libSpatialSEIR,
             which will allow the required random number generation
             to be supplied by R, or potentially an embedded version
             of the boost random library.")
endif()





include_directories("${CMAKE_SOURCE_DIR}/libSpatialSEIR/include"
                    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/FullConditionals"
                    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/DataStructures"
                    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/clBLAS"
                    "."
					"${SpatialSEIR_Boost_INCLUDE_DIRS}")



add_library(SpatialSEIR SHARED 
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/src/RandomNumberProvider.cpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/src/CompartmentalModelMatrix.cpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/src/CovariateMatrix.cpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/src/OCLProvider.cpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/src/OCLFunctions.cpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/src/ModelContext.cpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/src/ModelContextSummaryFunctions.cpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/src/FullConditionals/FullConditional.cpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/src/FullConditionals/FC_S0.cpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/src/FullConditionals/FC_E0.cpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/src/FullConditionals/FC_I0.cpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/src/FullConditionals/FC_R0.cpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/src/FullConditionals/FC_S_star.cpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/src/FullConditionals/FC_E_star.cpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/src/FullConditionals/FC_R_star.cpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/src/FullConditionals/FC_Beta.cpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/src/FullConditionals/FC_Beta_P_RS.cpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/src/FullConditionals/FC_Rho.cpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/src/FullConditionals/FC_Gamma_EI.cpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/src/FullConditionals/FC_Gamma_IR.cpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/src/FullConditionals/FC_Gamma.cpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/src/Samplers/CompartmentSamplers/CompartmentMetropolisSampler.cpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/src/Samplers/CompartmentSamplers/CompartmentMetropolisSampler.cpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/src/Samplers/CompartmentSamplers/IndexedCompartmentMetropolisSampler.cpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/src/Samplers/CompartmentSamplers/IndexedCompartmentSliceSampler.cpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/src/Samplers/InitCompartmentSamplers/InitCompartmentMetropolisSampler.cpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/src/Samplers/InitCompartmentSamplers/IndexedInitCompartmentMetropolisSampler.cpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/src/Samplers/InitCompartmentSamplers/IndexedInitCompartmentSliceSampler.cpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/src/Samplers/ParameterSamplers/ParameterSingleMetropolisSampler.cpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/src/Samplers/ParameterSamplers/ParameterJointMetropolisSampler.cpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/src/Samplers/ParameterSamplers/ParameterJointSliceSampler.cpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/src/Samplers/ParameterSamplers/ParameterDecorrelationSampler.cpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/src/Samplers.cpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/src/DistanceMatrix.cpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/clBLAS/clBLAS.h"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/clBLAS/clBLAS-complex.h"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/clBLAS/clBLAS.version.h"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/clBLAS/clAmdBlas.h"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/clBLAS/clAmdBlas.version.h" 
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/src/IOProvider.cpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/RandomNumberProvider.hpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/DataStructures/CompartmentalModelMatrix.hpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/DataStructures/CovariateMatrix.hpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/DataStructures/DistanceMatrix.hpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/OCLProvider.hpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/ModelContext.hpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/LSS_Samplers.hpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/FullConditionals/LSS_FullConditional.hpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/FullConditionals/LSS_FullConditionalList.hpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/FullConditionals/LSS_FC_S0.hpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/FullConditionals/LSS_FC_E0.hpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/FullConditionals/LSS_FC_I0.hpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/FullConditionals/LSS_FC_R0.hpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/FullConditionals/LSS_FC_S_star.hpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/FullConditionals/LSS_FC_E_star.hpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/FullConditionals/LSS_FC_R_star.hpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/FullConditionals/LSS_FC_Beta.hpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/FullConditionals/LSS_FC_Beta_P_RS.hpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/FullConditionals/LSS_FC_Rho.hpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/FullConditionals/LSS_FC_Gamma_EI.hpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/FullConditionals/LSS_FC_Gamma_IR.hpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/FullConditionals/LSS_FC_Gamma.hpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/IOProvider.hpp"
	"${SpatialSEIR_Boost_LIBRARIES}"
	)

add_library(CLBLAS1 SHARED IMPORTED)
add_library(CLBLAS2 SHARED IMPORTED)
add_library(CLBLAS3 SHARED IMPORTED)
set_property(TARGET CLBLAS1 PROPERTY IMPORTED_LOCATION ${CMAKE_BINARY_DIR}/clBLAS-prefix/src/clBLAS-build/package/lib64/libclBLAS.so)
set_property(TARGET CLBLAS2 PROPERTY IMPORTED_LOCATION ${CMAKE_BINARY_DIR}/clBLAS-prefix/src/clBLAS-build/package/lib64/libclBLAS.so.2)
set_property(TARGET CLBLAS3 PROPERTY IMPORTED_LOCATION ${CMAKE_BINARY_DIR}/clBLAS-prefix/src/clBLAS-build/package/lib64/libclBLAS.so.2.1.0)

target_link_libraries(SpatialSEIR CLBLAS1)
target_link_libraries(SpatialSEIR CLBLAS2)
target_link_libraries(SpatialSEIR CLBLAS3)

add_dependencies(SpatialSEIR clBLAS1)
add_dependencies(SpatialSEIR clBLAS2)
add_dependencies(SpatialSEIR clBLAS3)
add_dependencies(CLBLAS1 clBLAS)
add_dependencies(CLBLAS2 clBLAS)
add_dependencies(CLBLAS3 clBLAS)





add_definitions( -DLSS_KERNEL_DIRECTORY="${CMAKE_BINARY_DIR}/libSpatialSEIR_Kernels/" )

install (TARGETS SpatialSEIR DESTINATION lib)

install (FILES 
            "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/RandomNumberProvider.hpp"
            "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/DataStructures/CompartmentalModelMatrix.hpp"
            "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/DataStructures/CovariateMatrix.hpp"
            "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/DataStructures/DistanceMatrix.hpp"
            "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/clBLAS/clBLAS.h"
            "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/clBLAS/clBLAS-complex.h"
            "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/clBLAS/clBLAS.version.h"
            "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/clBLAS/clAmdBlas.h"
            "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/clBLAS/clAmdBlas.version.h"
            "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/OCLProvider.hpp"
            "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/ModelContext.hpp" 
            "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/LSS_Samplers.hpp"
            "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/IOProvider.hpp"
            "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/FullConditionals/LSS_FullConditional.hpp"
            "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/FullConditionals/LSS_FullConditionalList.hpp"
            "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/FullConditionals/LSS_FC_S0.hpp"
            "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/FullConditionals/LSS_FC_E0.hpp"
            "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/FullConditionals/LSS_FC_I0.hpp"
            "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/FullConditionals/LSS_FC_R0.hpp"
            "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/FullConditionals/LSS_FC_S_star.hpp"
            "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/FullConditionals/LSS_FC_E_star.hpp"
            "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/FullConditionals/LSS_FC_R_star.hpp"
            "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/FullConditionals/LSS_FC_Beta.hpp"
            "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/FullConditionals/LSS_FC_Beta_P_RS.hpp"
            "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/FullConditionals/LSS_FC_Rho.hpp"
            "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/FullConditionals/LSS_FC_Gamma_EI.hpp"
            "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/FullConditionals/LSS_FC_Gamma_IR.hpp"
            "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/FullConditionals/LSS_FC_Gamma.hpp"

            DESTINATION include)

add_custom_command(TARGET SpatialSEIR PRE_BUILD
                    COMMAND ${CMAKE_COMMAND} -E copy_directory
                    ${CMAKE_CURRENT_SOURCE_DIR}/libSpatialSEIR/src/kernels/
                    ${CMAKE_BINARY_DIR}/libSpatialSEIR_Kernels/
                 )




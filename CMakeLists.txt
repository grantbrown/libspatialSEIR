cmake_minimum_required(VERSION 2.8.7)

####################################
## General Configuration Options  ##
####################################
if (WIN32)
	SET(CMAKE_C_COMPILER "gcc")
	SET(CMAKE_CXX_COMPILER "g++")
	SET(BOOST_LOCATION "./LSS_Embedded_Boost")
	set(USE_BLAS "false")
	SET(R_INSTALL_DIRECTORY "C:/Program Files/R/R-3.1.1")
	SET(R_INCLUDE_DIR "${R_INSTALL_DIRECTORY}/include" )
	SET(RCPP_INCLUDE_DIR "~/Documents/R/win-library/3.1/Rcpp/include")
	#add_definitions(-D__cplusplus="true")
else()
	SET(BOOST_LOCATION "NOT_SET")
	SET(R_INSTALL_DIRECTORY "NOT_SET")
	set(USE_BLAS "true")
	add_definitions(-DLSS_USE_BLAS=1)
	SET(R_INCLUDE_DIR "NOT_SET")
	# Only windows users and those not using boost need to worry about the following two lines:
	SET(RCPP_INCLUDE_DIR "~/R/x86_64-pc-linux-gnu-library/3.0/Rcpp/include/")
	SET(R_INCLUDE_DIR "/usr/share/R/include/")
endif()


# Compile libspatialSEIR with OpenCL support?
SET(OPEN_CL_SUPPORT false)

# Use Boost? Alternative is to fall back and try to link to R. This can be difficult, 
# and is considerably slower. You probably want to use boost. 
SET(USE_BOOST true)

# libspatialSEIR requires the eigen3 library, so here we must 
# specify the location. By default, Eigen is supplied as a submodule.
SET(EIGEN3_INCLUDE_DIR "./eigen/")

project(libSpatialSEIR)

# Use stdio instead of Rcpp::Rcout
add_definitions( -DLSS_USE_STDIO="true" )

if (OPEN_CL_SUPPORT)
	message("Building libSpatialSEIR with OpenCL Support")
    add_definitions( -DENABLE_OPEN_CL="true" )
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
endif()

if (USE_BOOST)
    add_definitions( -DLSS_USE_BOOST="true" )
    ## Configure Boost
    SET(Boost_USE_STATIC_LIBS OFF) 
    SET(Boost_USE_MULTITHREADED ON)  
    SET(Boost_USE_STATIC_RUNTIME OFF) 
	
	if ("${BOOST_LOCATION}" STREQUAL "NOT_SET")
		SET(RCPP_INCLUDE_DIR ".")
		SET(R_INCLUDE_DIR ".")
		find_package(Boost 1.49.0 COMPONENTS random REQUIRED) 
		set(SpatialSEIR_Boost_LIBRARIES)
		set(SpatialSEIR_Boost_INCLUDE_DIRS)
		if(Boost_FOUND)
			set(SpatialSEIR_Boost_LIBRARIES ${Boost_LIBRARIES})
			set(SpatialSEIR_Boost_INCLUDE_DIRS ${Boost_INCLUDE_DIRS})
		else()
			message("Error: Boost was not found on the system. If you think this 
			is a mistake, you might try setting BOOST_ROOT etc (there are probably other 
			messages on your screen right now to the same effect). If you think 
			this is a mistake and are using windows, please either move the boost_1_56_0 
			folder to the libSpatialSEIR root, or update CMakeLists.txt to reflect its location.")
		endif()
	else()
		message("Setting boost location manually: ${BOOST_LOCATION}")
		SET(BOOST_ROOT "${BOOST_LOCATION}")
		SET(SpatialSEIR_Boost_LIBRARIES "${BOOST_LOCATION}/libs")
		SET(SpatialSEIR_Boost_INCLUDE_DIRS "${BOOST_LOCATION}")
	endif()
   
else()
    if ("${RCPP_INCLUDE_DIR}" STREQUAL "")
        message("ERROR: if not using boost, you must set the location of Rcpp.h via RCPP_INCLUDE_DIR")
    endif()
    if ("${R_INCLUDE_DIR}" STREQUAL "")
        message("ERROR: if not using boost, you must set the location of R.h via R_INCLUDE_DIR")
    endif()
    set(SpatialSEIR_Boost_LIBRARIES ".")
    set(SpatialSEIR_Boost_INCLUDE_DIRS ".")
endif()


if (WIN32)
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic -Wextra -O2")
endif()
## Configure BLAS
set(BLA_STATIC OFF)
if ("${USE_BLAS}" STREQUAL "true")

	find_package(BLAS)
	add_definitions(-DLSS_CBLAS_DGEMV_NAME=cblas_dgemv)
	add_definitions(-DLSS_CBLAS_DGEMM_NAME=cblas_dgemm)
	if (BLAS_FOUND)
		message("BLAS found: ${BLAS_LIBRARIES}")  
		message("            ${BLAS95_LIBRARIES}")
		message("BLAS linker flags: ${BLAS_LINKER_FLAGS}")
	else()
		message("Error: Blas not found.")
	endif()
else()
	message("BLAS Disabled: libSpatialSEIR will use Eigen")	
endif()

message("${SpatialSEIR_Boost_INCLUDE_DIRS}")
message("${RCPP_INCLUDE_DIR}")
message("${R_INCLUDE_DIR}")
## Require Boost Message
include_directories("${CMAKE_SOURCE_DIR}/libSpatialSEIR/include"
                    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/FullConditionals"
                    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/FullConditionals/DataModels"
                    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/FullConditionals/DataModels/Overdispersion"
                    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/DataStructures"
                    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/clBLAS"
                    "."
					"${SpatialSEIR_Boost_INCLUDE_DIRS}"
                    "${EIGEN3_INCLUDE_DIR}"
                    "${RCPP_INCLUDE_DIR}"
                    "${R_INCLUDE_DIR}")


add_library(SpatialSEIR SHARED 
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/src/RandomNumberProvider.cpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/src/CompartmentalModelMatrix.cpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/src/CovariateMatrix.cpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/src/OCLProvider.cpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/src/OCL/Deterministic.cpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/src/OCL/FullConditionals/FC_R_Star.cpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/src/ModelContext.cpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/src/ModelContextSummaryFunctions.cpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/src/FullConditionals/FullConditional.cpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/src/FullConditionals/FC_S0.cpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/src/FullConditionals/FC_E0.cpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/src/FullConditionals/FC_I0.cpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/src/FullConditionals/FC_R0.cpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/src/FullConditionals/FC_S_star.cpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/src/FullConditionals/FC_E_star.cpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/src/FullConditionals/FC_I_star.cpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/src/FullConditionals/FC_R_star.cpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/src/FullConditionals/FC_Beta.cpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/src/FullConditionals/FC_P_SE.cpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/src/FullConditionals/FC_Beta_P_RS.cpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/src/FullConditionals/FC_Rho.cpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/src/FullConditionals/FC_Gamma_EI.cpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/src/FullConditionals/FC_Gamma_IR.cpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/src/FullConditionals/DataModels/Overdispersion/FC_I_star_overdispersion.cpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/src/FullConditionals/DataModels/Overdispersion/FC_R_star_overdispersion.cpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/src/FullConditionals/DataModels/Overdispersion/FC_Phi.cpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/src/IterationTasks/SetCompartmentSamplingIndicesTask.cpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/src/IterationTasks/PerformDecorrelationStep.cpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/src/Samplers/CompartmentSamplers/CompartmentMetropolisSampler.cpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/src/Samplers/CompartmentSamplers/CompartmentBinomialMetropolisSampler.cpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/src/Samplers/CompartmentSamplers/IndexedCompartmentBinomialMetropolisSampler.cpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/src/Samplers/CompartmentSamplers/CompartmentBinomialSliceSampler.cpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/src/Samplers/CompartmentSamplers/CompartmentBinomialMixedSampler.cpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/src/Samplers/CompartmentSamplers/CompartmentMetropolisSampler_OCL.cpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/src/Samplers/CompartmentSamplers/IndexedCompartmentMetropolisSampler.cpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/src/Samplers/CompartmentSamplers/IndexedCompartmentSliceSampler.cpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/src/Samplers/InitCompartmentSamplers/InitCompartmentMetropolisSampler.cpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/src/Samplers/InitCompartmentSamplers/InitCompartmentMetropolisSampler_OCL.cpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/src/Samplers/InitCompartmentSamplers/IndexedInitCompartmentMetropolisSampler.cpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/src/Samplers/InitCompartmentSamplers/IndexedInitCompartmentSliceSampler.cpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/src/Samplers/ParameterSamplers/ParameterSingleMetropolisSampler.cpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/src/Samplers/ParameterSamplers/ParameterJointMetropolisSampler.cpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/src/Samplers/ParameterSamplers/ParameterJointMetropolisSampler_OCL.cpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/src/Samplers/ParameterSamplers/ParameterJointSliceSampler.cpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/src/Samplers/ParameterSamplers/ParameterDecorrelationSampler.cpp"
	"${CMAKE_SOURCE_DIR}/libSpatialSEIR/src/IterationTasks/PerformHybridSE_EI_UpdateStep.cpp"
	"${CMAKE_SOURCE_DIR}/libSpatialSEIR/src/IterationTasks/PerformHybridIR_RS_UpdateStep.cpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/src/Samplers/ParameterSamplers/ParameterNullSampler.cpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/src/Samplers/ParameterSamplers/ParameterHybridSampler.cpp"
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
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/LSS_IterationTasks.hpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/FullConditionals/LSS_FullConditional.hpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/FullConditionals/LSS_FullConditionalList.hpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/FullConditionals/LSS_FC_S0.hpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/FullConditionals/LSS_FC_E0.hpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/FullConditionals/LSS_FC_I0.hpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/FullConditionals/LSS_FC_R0.hpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/FullConditionals/LSS_FC_S_star.hpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/FullConditionals/LSS_FC_E_star.hpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/FullConditionals/LSS_FC_I_star.hpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/FullConditionals/LSS_FC_R_star.hpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/FullConditionals/LSS_FC_Beta.hpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/FullConditionals/LSS_FC_P_SE.hpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/FullConditionals/LSS_FC_Beta_P_RS.hpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/FullConditionals/LSS_FC_Rho.hpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/FullConditionals/LSS_FC_Gamma_EI.hpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/FullConditionals/LSS_FC_Gamma_IR.hpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/FullConditionals/DataModels/Overdispersion/LSS_FC_I_star_overdispersion.hpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/FullConditionals/DataModels/Overdispersion/LSS_FC_R_star_overdispersion.hpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/FullConditionals/DataModels/Overdispersion/LSS_FC_Phi.hpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/IOProvider.hpp"
	"${SpatialSEIR_Boost_LIBRARIES}"
	)
	
if (OPEN_CL_SUPPORT)
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

endif()

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
            "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/LSS_IterationTasks.hpp"
            "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/IOProvider.hpp"
            "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/FullConditionals/LSS_FullConditional.hpp"
            "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/FullConditionals/LSS_FullConditionalList.hpp"
            "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/FullConditionals/LSS_FC_S0.hpp"
            "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/FullConditionals/LSS_FC_E0.hpp"
            "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/FullConditionals/LSS_FC_I0.hpp"
            "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/FullConditionals/LSS_FC_R0.hpp"
            "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/FullConditionals/LSS_FC_S_star.hpp"
            "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/FullConditionals/LSS_FC_E_star.hpp"
            "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/FullConditionals/LSS_FC_I_star.hpp"
            "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/FullConditionals/LSS_FC_R_star.hpp"
            "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/FullConditionals/LSS_FC_Beta.hpp"
            "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/FullConditionals/LSS_FC_P_SE.hpp"
            "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/FullConditionals/LSS_FC_Beta_P_RS.hpp"
            "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/FullConditionals/LSS_FC_Rho.hpp"
            "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/FullConditionals/LSS_FC_Gamma_EI.hpp"
            "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/FullConditionals/LSS_FC_Gamma_IR.hpp"
            "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/FullConditionals/DataModels/Overdispersion/LSS_FC_I_star_overdispersion.hpp"
            "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/FullConditionals/DataModels/Overdispersion/LSS_FC_R_star_overdispersion.hpp"
            "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/FullConditionals/DataModels/Overdispersion/LSS_FC_Phi.hpp"
            DESTINATION include)

add_custom_command(TARGET SpatialSEIR PRE_BUILD
                    COMMAND ${CMAKE_COMMAND} -E copy_directory
                    ${CMAKE_CURRENT_SOURCE_DIR}/libSpatialSEIR/src/kernels/
                    ${CMAKE_BINARY_DIR}/libSpatialSEIR_Kernels/
                 )




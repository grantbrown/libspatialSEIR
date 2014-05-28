cmake_minimum_required(VERSION 2.8.7)

project(libSpatialSEIR)

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
                    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/DataStructures"
                    "."
					"${SpatialSEIR_Boost_INCLUDE_DIRS}")



add_library(SpatialSEIR SHARED 
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/src/RandomNumberProvider.cpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/src/CompartmentalModelMatrix.cpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/src/CovariateMatrix.cpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/src/OCLProvider.cpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/src/ModelContext.cpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/src/ModelContextSummaryFunctions.cpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/src/FullConditional.cpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/src/DistanceMatrix.cpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/src/IOProvider.cpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/RandomNumberProvider.hpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/DataStructures/CompartmentalModelMatrix.hpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/DataStructures/CovariateMatrix.hpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/DataStructures/DistanceMatrix.hpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/OCLProvider.hpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/ModelContext.hpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/FullConditional.hpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/IOProvider.hpp"
	"${SpatialSEIR_Boost_LIBRARIES}"
	)

add_definitions( -DLSS_KERNEL_DIRECTORY="${CMAKE_BINARY_DIR}/libSpatialSEIR_Kernels/" )

install (TARGETS SpatialSEIR DESTINATION lib)
install (FILES 
            "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/RandomNumberProvider.hpp"
            "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/DataStructures/CompartmentalModelMatrix.hpp"
            "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/DataStructures/CovariateMatrix.hpp"
            "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/DataStructures/DistanceMatrix.hpp"
            "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/OCLProvider.hpp"
            "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/ModelContext.hpp"
            "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/FullConditional.hpp"
            "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/IOProvider.hpp"
            DESTINATION include)

install (FILES 
            "${CMAKE_SOURCE_DIR}/libSpatialSEIR/src/kernels/test_kernel.cl" 
            DESTINATION bin/libSpatialSEIR_Kernels)

add_custom_command(TARGET SpatialSEIR POST_BUILD
                    COMMAND ${CMAKE_COMMAND} -E copy
                    ${CMAKE_CURRENT_SOURCE_DIR}/libSpatialSEIR/src/kernels/test_kernel.cl
                    ${CMAKE_BINARY_DIR}/libSpatialSEIR_Kernels/test_kernel.cl
                 )



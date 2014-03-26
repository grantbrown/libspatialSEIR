cmake_minimum_required(VERSION 2.8.7)

project(libSpatialSEIR)

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
endif()

include_directories("${CMAKE_SOURCE_DIR}/libSpatialSEIR/include"
                    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/DataStructures"
                    "."
					"${SpatialSEIR_Boost_INCLUDE_DIRS}"
)



add_library(SpatialSEIR SHARED 
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/src/RandomNumberProvider.cpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/src/CompartmentalModelMatrix.cpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/src/CovariateMatrix.cpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/src/OCLProvider.cpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/src/FullConditional.cpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/src/ModelContext.cpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/src/DistanceMatrix.cpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/RandomNumberProvider.hpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/DataStructures/CompartmentalModelMatrix.hpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/DataStructures/CovariateMatrix.hpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/DataStructures/DistanceMatrix.hpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/OCLProvider.hpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/ModelContext.hpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/FullConditional.hpp"
	"${SpatialSEIR_Boost_LIBRARIES}"
	)


install (TARGETS SpatialSEIR DESTINATION lib)
install (FILES 
            "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/RandomNumberProvider.hpp"
            "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/DataStructures/CompartmentalModelMatrix.hpp"
            "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/DataStructures/CovariateMatrix.hpp"
            "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/DataStructures/DistanceMatrix.hpp"
            "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/OCLProvider.hpp"
            "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/ModelContext.hpp"
            "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/FullConditional.hpp"
            DESTINATION include)




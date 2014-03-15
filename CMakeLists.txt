cmake_minimum_required(VERSION 2.8.7)

project(libSpatialSEIR)


include_directories("${CMAKE_SOURCE_DIR}/libSpatialSEIR/include"
                    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/DataStructures"
                    ".")




add_library(SpatialSEIR SHARED 
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/src/CompartmentalModelMatrix.cpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/src/CovariateMatrix.cpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/src/OCLProvider.cpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/src/FullConditional.cpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/src/ModelContext.cpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/DataStructures/CompartmentalModelMatrix.hpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/DataStructures/CovariateMatrix.hpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/OCLProvider.hpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/ModelContext.hpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/FullConditional.hpp"
    )


install (TARGETS SpatialSEIR DESTINATION lib)
install (FILES 
            "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/DataStructures/CompartmentalModelMatrix.hpp"
            "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/DataStructures/CovariateMatrix.hpp"
            "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/OCLProvider.hpp"
            "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/ModelContext.hpp"
            "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/FullConditional.hpp"
            DESTINATION include)




cmake_minimum_required(VERSION 2.8.7)

project(libSpatialSEIR)


include_directories("${CMAKE_SOURCE_DIR}/libSpatialSEIR/include"
                    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/DataStructures"
                    ".")




add_library(libSpatialSEIR SHARED 
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/src/CompartmentalModelMatrix.cpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/src/CovariateMatrix.cpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/DataStructures/CompartmentalModelMatrix.hpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/DataStructures/CovariateMatrix.hpp")

install (TARGETS libSpatialSEIR DESTINATION bin)
install (FILES 
            "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/DataStructures/CompartmentalModelMatrix.hpp"
            "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/DataStructures/CovariateMatrix.hpp"
            DESTINATION include)




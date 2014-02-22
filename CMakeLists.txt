cmake_minimum_required(VERSION 2.8.7)

project(libSpatialSEIR)


include_directories("${CMAKE_SOURCE_DIR}/libSpatialSEIR/include"
                    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/DataStructures"
                    ".")

add_subdirectory(libSpatialSEIR)






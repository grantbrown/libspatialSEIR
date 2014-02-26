add_library(libSpatialSEIR SHARED 
            CompartmentalModelMatrix.cpp 
            CovariateMatrix.cpp
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/DataStructures/CompartmentalModelMatrix.hpp"
    "${CMAKE_SOURCE_DIR}/libSpatialSEIR/include/DataStructures/CovariateMatrix.hpp")

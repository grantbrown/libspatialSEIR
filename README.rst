libspatialSEIR
===============

Disclaimer
-----------
This software is under heavy development, and as such it's full of jagged edges, changing API's, and broken things.
You have been warned. 

Introduction
---------------

libspatialSEIR is a C++/OpenCL based framework for fitting Bayesian spatial SEIR and SEIRS compartmental epidemic models.
The numerical heavy lifting is all done in C++ with the optional use of OpenCL calls to multi-core CPU's and 
GPU's to accelerate computation. The primary interface is provided via the included R package, which uses the Rcpp 
library to access the lower level C++ API. 

The most up to date documentation is hosted via github pages here_.

.. _here: http://grantbrown.github.io/libspatialSEIR/

In addition, example R code is available in the /scripts directory.  

Dependencies 
-------------
While a CRAN friendly version of libspatialSEIR is tentatively planned, the project currently requires some prerequisites:

1. Git
2. CMake
3. An appropriate C++ compiler 
4. Boost random >= 1.4.9 
5. R >= 3.0
6. A BLAS capable library

In addition, if you would like to enable OpenCL support, you must install the appropriate OpenCL SDK for your hardware. This 
will also require that you tweak the CMakeLists.txt file present in the project root directory. 

Installation - Linux
-----------------------
The installation process is still very much in "development" mode. As of now, it has only been tested on Linux Mint, but similar instructions 
should apply for other platforms. After installing the appropriate pre-requisites, clone a copy of the repository. In order to make sure you 
get a copy of the embedded version of clBLAS (if you're using OpenCL) and Eigen3 (no matter what), be sure to clone the repository recursively:

::
    
    git clone https://github.com/grantbrown/libspatialSEIR.git --recursive


Once you have a copy of the code, the easiest way to install the library is to run the install script in the R/release_package folder: 


::
    
    cd R
    cd release_package
    ./install.sh


Assuming this step completed correctly, you're ready to use the library from R. Should you wish to contribute to the development of the
C++ library, it's probably better to use the "dev_package" instead:


:: 

    mkdir build
    cd build
    cmake ../
    make
    sudo make install

    cd ../R/
    R CMD INSTALL dev_package

This allows you to re-build the C++ shared library without re-compiling the R package (unless you make API changes, in which case both must be re-compiled).
This process may also be useful for troubleshooting errors encountered in the build process. 


Installation - Windows
-------------------------
I'm currently working on getting a good build process for Windows, and plan to eventually release pre-compiled binaries. In the meantime, if you happen to be a CMake wizard, you should probably be able to adapt the above steps to suit your needs.   


Installation - OSX
-------------------------
Without access to an OSX machine, I can't troubleshoot the install process. This is planned for the future, though as a unix-alike the build process should closely resemble that of linux machines. 




That should be it, though as previously mentioned this is still very much a work in progress. Assuming it worked, you could try out some of
the example scripts. 

::     

    cd ../scripts
    R
    > source("./singleLocationTest.R")



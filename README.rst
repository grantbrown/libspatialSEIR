libspatialSEIR
===============

Introduction
---------------

libspatialSEIR is a C++/OpenCL based framework for fitting Bayesian spatial SEIR and SEIRS compartmental epidemic models.
The numerical heavy lifting is all done in C++ with the optional use of OpenCL calls to multi-core CPU's and 
GPU's to accelerate computation. The primary interface is provided via the included R package, which uses the Rcpp 
library to access the lower level C++ API. 

At this point documentation is pretty limited, but a description of the class of models which can be fit
is available in /doc/models/SEIR_ALGO.pdf

In addition, example R code is available in the /scripts directory.  

Dependencies 
-------------
While a CRAN friendly version of libspatialSEIR is tentatively planned, the project currently requires some prerequisites:

1. Git
2. CMake
3. An appropriate C++ compiler 
4. Boost random >= 1.4.9 
5. An OpenCL distribution
6. (preferably) R >= 3.0



Installation
-------------
The installation process is still very much in "development" mode. As of now, it has only been tested on Linux Mint, but similar instructions 
should apply for other platforms. After installing the appropriate pre-requisites, clone a copy of the repository. In order to make sure you 
get a copy of the embedded version of clBLAS, be sure to clone the repository recursively:


::
    
    git clone https://github.com/grantbrown/libspatialSEIR.git --recursive



Once you have a copy of the code, change to the project directory and create a build folder. Use CMake to generate the 
appropriate build files. 


::
    
    cd libspatialSEIR
    mkdir build
    cd build
    cmake ../


Assuming this step completed correctly, make and install the library. At some point it should be feasible to just make the library
and link to it, but first I need to figure out how to make R play nice with shared objects being placed in non-standard directories. 


:: 
    
    make
    sudo make install

Next you'll need to build the R package itself. 

::
    
    cd ../R
    R CMD INSTALL package

That should be it, though as previously mentioned this is still very much a work in progress. Assuming it worked, you could try out some of
the example scripts. 

:: 
    
    cd ../scripts
    R
    > source("./singleLocationTest.R")



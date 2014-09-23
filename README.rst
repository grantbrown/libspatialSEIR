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

Installation
----------------

The most up to date installation information is available on the libspatialSEIR wiki_.

.. _wiki: https://github.com/grantbrown/libspatialSEIR/wiki/Installation

So far I've only set this up on linux (Ubuntu/Mint), but support for Windows and OSX is planned. 


That should be it, though as previously mentioned this is still very much a work in progress. Assuming it worked, you could try out some of
the example scripts. 

::     

    cd ../scripts
    R
    > source("./singleLocationTest.R")



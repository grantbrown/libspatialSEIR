# Obtain a copy of the development package
cp -r ../dev_package/* ./
# Create a folder for the libspatialSEIR code. 
cd src
mkdir LSS
cd include
mkdir LSS
cd ..
# Create appropriate makefile using cmake.
cp -r ../../../libSpatialSEIR/src/* ./LSS 
cp -r ../../../libSpatialSEIR/include/* ./include/LSS
#cmake -G "Unix Makefiles" ../../../../
#cmake -G "MinGW Makefiles" ../../../../
cd ../
# Replace Makevars
mv ./src/Makevars_release.win ./src/Makevars.win
# Move up one directory and install
cd ..
R CMD INSTALL release_package --library=~/Documents/R/win-library/3.1
#R CMD INSTALL release_package



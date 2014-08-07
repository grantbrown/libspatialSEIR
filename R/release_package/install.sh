# Obtain a copy of the development package
cp ../dev_package/* ./ -r
# Create a folder for the libspatialSEIR code. 
mkdir -p ./src/LSS
# Create appropriate makefile using cmake.
cd ./src/LSS
cmake ../../../../
cd ../../
# Replace Makevars
mv ./src/Makevars_release ./src/Makevars
# Move up one directory and install
cd ..
R CMD INSTALL release_package



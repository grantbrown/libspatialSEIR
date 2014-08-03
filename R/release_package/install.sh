# Obtain a copy of the development package
cp ../dev_package/* ./ -r
# Create a folder for the libspatialSEIR code. 
mkdir ./src/LSS
# Copy libspatialSEIR code 
cp ../../libSpatialSEIR/* ./src/LSS/ -r
# Replace Makevars
mv ./src/Makevars_release ./src/Makevars
# Move up one directory and install
cd ..
R CMD INSTALL release_package



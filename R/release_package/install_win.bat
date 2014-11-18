cp -r ../dev_package/* ./
cd src
mkdir LSS
cd include
mkdir LSS
cd ..
cp -r ../../../libSpatialSEIR/src/* ./LSS 
cp -r ../../../libSpatialSEIR/include/* ./include/LSS
cd ../
mv ./src/Makevars_release.win ./src/Makevars.win
cd ..
R CMD INSTALL release_package


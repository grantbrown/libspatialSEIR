#include<iostream>
#include<fstream>
#include<cstdio>

#ifndef IO_PROVIDER_INC
#define IO_PROVIDER_INC

namespace SpatialSEIR
{
    //Forward declare required classes
    class ModelContext;

    class IOProvider
    {
        public:
            //methods
            //Empty Constructor
            IOProvider();
            //Full Constructor
            IOProvider(ModelContext* context,
                       std::string* outFilePath, 
                       int* variableList,
                       int* iterationStride);
            //Initialize
            int populate(ModelContext* context,
                         std::string* outFilePath,
                         int* variableList,
                         int* iterationStride);
            int close();
            int fileInit();
            int catIter(int iteration);
            ~IOProvider();

            //Attributes
            ModelContext** context;
            int* variableList;
            int* iterationStride;
            bool* isOpen;
            std::ofstream* outFileStream;
            std::string* outFilePath;
    };

}

#endif

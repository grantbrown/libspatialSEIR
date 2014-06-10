#include<iostream>
#include<fstream>
#include<vector>
#include<cstdio>
#include<time.h>

#ifndef IO_PROVIDER_INC
#define IO_PROVIDER_INC

namespace SpatialSEIR
{
    //Forward declare required classes
    class ModelContext;

    struct TimeLocationTrace
    {
        int locationIndex;
        int timeIndex;
    };

    class IOProvider
    {
        public:
            //methods
            //Empty Constructor
            IOProvider();
            //Full Constructor
            IOProvider(ModelContext* context,
                       std::string* outFilePath, 
                       int* iterationStride);
            //Initialize
            int populate(ModelContext* context,
                         std::string* outFilePath,
                         int* iterationStride);
            //Methods
            void setTrace(int locationIndex);
            void setTrace(int locationIndex, int timeIndex);
            int close();
            int fileInit();
            int catIter(int iteration);
            ~IOProvider();

            //Attributes
            ModelContext** context;
            int* iterationStride;
            bool* isOpen;
            time_t* timer;
            time_t* startTime;
            std::vector<TimeLocationTrace*> *timeLocationTraces;
            std::ofstream* outFileStream;
            std::string* outFilePath;
    };

}

#endif

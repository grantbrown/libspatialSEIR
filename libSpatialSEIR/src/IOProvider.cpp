#include <IOProvider.hpp>
#include <ModelContext.hpp>


namespace SpatialSEIR
{

    IOProvider::IOProvider(ModelContext* _context,
                           std::string* _outFilePath,
                           int* _variableList,
                           int* _iterationStride)
    {
        context = new ModelContext*;
        iterationStride = new int;
        variableList = new int[10]; // This will need to be generalized at some point.
        outFilePath = new std::string; 
        *context = _context; 
        *outFilePath = *_outFilePath;
        *iterationStride = *_iterationStride;
        outFileStream = new std::ofstream;
        isOpen = new bool;
        *isOpen = false;
        int i;
        for (i = 0; i < 10; i++)
        {
            variableList[i] = _variableList[i];
        } 
        this -> fileInit();
    }

    int IOProvider::fileInit()
    {
        // Clear file 
        FILE* tmp = fopen(outFilePath -> c_str(), "w");
        fclose(tmp);

        // Open file as output stream
        *isOpen = true;
        outFileStream -> open(outFilePath -> c_str());

        // Write header.
        // Is there a more concise way to code this?
        // Write Iteration header
        if (variableList[0] != 0)
        {
            // Write beta header
        }
        if (variableList[1] != 0)
        {
            // Write rho header
        }
        if (variableList[2] != 0)
        {
            // Write p_se header
        }
        if (variableList[3] != 0)
        {
            // Write p_ei header
        }
        if (variableList[4] != 0)
        {
            // Write p_ir header
        }
        if (variableList[5] != 0)
        {
            // Write p_rs header
        }
        if (variableList[6] != 0)
        {
            // Write S* header
        }
        if (variableList[7] != 0)
        {
            // Write E* header
        }
        if (variableList[8] != 0)
        {
            // Write I* header
        }
        if (variableList[9] != 0)
        {
            // Write R* header
        }
        //Newline
        return(0);
    }
    int IOProvider::catIter(int iteration)
    {
        //Not Implemented 
        throw(-1);
    }
    int IOProvider::close()
    {
        outFileStream -> close();
        return(0);
    }

    IOProvider::~IOProvider()
    {
        if (*isOpen)
        {
            this -> close();
            delete isOpen;
        }
        delete context;
        delete iterationStride;
        delete[] variableList;
        delete outFilePath;
        delete outFileStream;
    }

}

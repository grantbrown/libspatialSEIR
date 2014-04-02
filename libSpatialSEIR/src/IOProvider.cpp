#include <IOProvider.hpp>
#include <ModelContext.hpp>
#include <CompartmentalModelMatrix.hpp>
#include <CovariateMatrix.hpp>


namespace SpatialSEIR
{


    IOProvider::IOProvider()
    {
        // Empty constructor, do nothing. 
    }

    IOProvider::IOProvider(ModelContext* _context,
                           std::string* _outFilePath,
                           int* _variableList,
                           int* _iterationStride)
    {
        this -> populate(&*_context,&*_outFilePath, &*_variableList, 
                &*_iterationStride);
    }


    int IOProvider::populate(ModelContext* _context,
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
        return(this -> fileInit());
    }

    int IOProvider::fileInit()
    {
        int i, j;

        int nLoc = *((*context)->S->ncol);
        int nTpt = *((*context)->S->ncol);

        // Clear file 
        FILE* tmp = fopen(outFilePath -> c_str(), "w");
        fclose(tmp);

        // Open file as output stream
        *isOpen = true;
        outFileStream -> open(outFilePath -> c_str());

        // Is there a more concise way to code this?
        if (variableList[0] != 0)
        {
            int betaLen = (*((*context) -> X -> ncol_x)) + (*((*context) -> X -> ncol_z)); 
            for (i = 0; i < betaLen; i++)
            {
                (*outFileStream) << "B" << i << ", "; 
            }    
        }
        if (variableList[1] != 0)
        {
            // Write rho header
            (*outFileStream) << "rho,";
        
        }
        if (variableList[2] != 0)
        {
            // Write p_se header
            for (i = 0; i < nLoc; i++)
            {
                for (j = 0; j < nTpt; j++)
                {
                    (*outFileStream) << "pSE_" << i << "_" << j << ",";
                }
            }
        }
        if (variableList[3] != 0)
        {
            // Write p_ei header
            (*outFileStream) << "p_ei" << ",";
        }
        if (variableList[4] != 0)
        {
            // Write p_ir header
            (*outFileStream) << "p_ir" << ",";

        }
        if (variableList[5] != 0)
        {
            // Write p_rs header
            for (j = 0; j < nTpt; j++)
            {
                (*outFileStream) << "pRS_" << j << ",";
            }

        }
        if (variableList[6] != 0)
        {
            // Write S* header
            for (i = 0; i < nLoc; i++)
            {
                for (j = 0; j < nTpt; j++)
                {
                    (*outFileStream) << "Sstar" << i << "_" << j << ",";
                }
            }
        }
        if (variableList[7] != 0)
        {
            // Write E* header
            for (i = 0; i < nLoc; i++)
            {
                for (j = 0; j < nTpt; j++)
                {
                    (*outFileStream) << "Estar" << i << "_" << j << ",";
                }
            }

        }
        if (variableList[8] != 0)
        {
            // Write I* header
            for (i = 0; i < nLoc; i++)
            {
                for (j = 0; j < nTpt; j++)
                {
                    (*outFileStream) << "Istar" << i << "_" << j << ",";
                }
            }

        }
        if (variableList[9] != 0)
        {
            // Write R* header
            for (i = 0; i < nLoc; i++)
            {
                for (j = 0; j < nTpt; j++)
                {
                    (*outFileStream) << "Rstar" << i << "_" << j << ",";
                }
            }
        }
        (*outFileStream) << "Iteration\n";
        // Write iteration number
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

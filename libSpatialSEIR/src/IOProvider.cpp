#include <IOProvider.hpp>
#include <ModelContext.hpp>
#include <CompartmentalModelMatrix.hpp>
#include <CovariateMatrix.hpp>
#include<time.h>


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

    IOProvider::~IOProvider()
    {
        if (*isOpen)
        {
            this -> close();
            delete isOpen;
        }
        delete startTime;
        delete timer;
        delete iterationStride;
        delete[] variableList;
        delete outFilePath;
        delete outFileStream;
        delete context;
    }

    int IOProvider::populate(ModelContext* _context,
                           std::string* _outFilePath,
                           int* _variableList,
                           int* _iterationStride)
    {
        context = new ModelContext*;
        iterationStride = new int;
        variableList = new int[31]; // This will need to be generalized at some point.
        outFilePath = new std::string; 
        *context = _context; 
        *outFilePath = *_outFilePath;
        *iterationStride = *_iterationStride;
        outFileStream = new std::ofstream;
        isOpen = new bool;
        *isOpen = false;
        int i;
        for (i = 0; i < 30; i++)
        {
            variableList[i] = _variableList[i];
        } 
        timer = new time_t;
        startTime = new time_t;
        time(startTime); // Set start time
        return(this -> fileInit());
    }

    int IOProvider::fileInit()
    {
        int i, j;
        int nLoc = *((*context)->S->ncol);
        int nTpt = *((*context)->S->nrow);

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
            // Write rho header
            for (j = 0; j < nTpt; j++)
            {
                (*outFileStream) << "gamma_" << j << ",";        
            }
        }
        if (variableList[3] != 0)
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
        if (variableList[4] != 0)
        {
            // Write p_ei header
            (*outFileStream) << "p_ei" << ",";
        }
        if (variableList[5] != 0)
        {
            // Write p_ir header
            (*outFileStream) << "p_ir" << ",";

        }
        if (variableList[6] != 0)
        {
            // Write p_rs header
            for (j = 0; j < nTpt; j++)
            {
                (*outFileStream) << "pRS_" << j << ",";
            }

        }
        if (variableList[7] != 0)
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
        if (variableList[8] != 0)
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
        if (variableList[9] != 0)
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
        if (variableList[10] != 0)
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
        if (variableList[11] != 0)
        {
            // Write S total Header
            (*outFileStream) << "S_Total,";
        }
        if (variableList[12] != 0)
        {
            // Write E total Header
            (*outFileStream) << "E_Total,";
        }
        if (variableList[13] != 0)
        {
            // Write I total Header
            (*outFileStream) << "I_Total,";
        }
        if (variableList[14] != 0)
        {
            // Write R total Header
            (*outFileStream) << "R_Total,";
        }
        if (variableList[15] != 0)
        {
            // Write S_star total Header
            (*outFileStream) << "S_star_Total,";
        }
        if (variableList[16] != 0)
        {
            // Write E_star total Header
            (*outFileStream) << "E_star_Total,";
        }
        if (variableList[17] != 0)
        {
            // Write I_star total Header
            (*outFileStream) << "I_star_Total,";
        }
        if (variableList[18] != 0)
        {
            // Write R_star total Header
            (*outFileStream) << "R_star_Total,";
        }
        if (variableList[19] != 0)
        {
            // Write average p_se header
            (*outFileStream) << "avgP_se,";
        }
        if (variableList[20] != 0)
        {
            // Write average p_se header
            (*outFileStream) << "avgP_rs,";
        }

    
        // Time specific

        if (variableList[21] != 0)
        {
            // Write S total Header
            for (j = 0; j < nTpt; j++)
            {
                (*outFileStream) << "S_Total_" << j << ",";
            }
        }
        if (variableList[22] != 0)
        {
            // Write E total Header
            for (j = 0; j < nTpt; j++)
            {
                (*outFileStream) << "E_Total_" << j << ",";
            }
        }
        if (variableList[23] != 0)
        {
            // Write I total Header
            for (j = 0; j < nTpt; j++)
            {
                (*outFileStream) << "I_Total_" << j << ",";
            }
        }
        if (variableList[24] != 0)
        {
            // Write R total Header
            for (j = 0; j < nTpt; j++)
            {
                (*outFileStream) << "R_Total_" << j << ",";
            }
        }
        if (variableList[25] != 0)
        {
            // Write S_star total Header
            for (j = 0; j < nTpt; j++)
            {
                (*outFileStream) << "S_star_Total_" << j << ",";
            }
        }
        if (variableList[26] != 0)
        {
            // Write E_star total Header
            for (j = 0; j < nTpt; j++)
            {
                (*outFileStream) << "E_star_Total_" << j << ",";
            }
        }
        if (variableList[27] != 0)
        {
            // Write I_star total Header            
            for (j = 0; j < nTpt; j++)
            {
                (*outFileStream) << "I_star_Total_" << j << ",";
            }
        }
        if (variableList[28] != 0)
        {
            // Write R_star total Header
            for (j = 0; j < nTpt; j++)
            {
                (*outFileStream) << "R_star_Total_" << j << ",";
            }
        }
        if (variableList[29] != 0)
        {
            // Write average p_se header
            for (j = 0; j < nTpt; j++)
            {
                (*outFileStream) << "avgP_se_" << j << ",";
            }
        }
        if (variableList[30] != 0)
        {
            // Write r_0 header. 
            (*outFileStream) << "r_0, ";
        }
        if (variableList[31] != 0)
        {
            // Write r_0_t header
            for (j = 0; j < nTpt; j++)
            {
                (*outFileStream) << "r_0_" << j << ",";
            }
        }

        (*outFileStream) << "Iteration,Time\n";
        // Write iteration number
        //Newline
        return(0);
    }
    int IOProvider::catIter(int iteration)
    {
        int i, j;
        int nLoc = *((*context)->S->ncol);
        int nTpt = *((*context)->S->ncol);

        if ((iteration % (*iterationStride)) != 0)
        {
            return(1);
        }

        // Is there a more concise way to code this?
        if (variableList[0] != 0)
        {
            int betaLen = (*((*context) -> X -> ncol_x)) + (*((*context) -> X -> ncol_z)); 
            for (i = 0; i < betaLen; i++)
            {
                (*outFileStream) << ((*context) -> beta)[i] << ", "; 
            }    
        }

        if (variableList[1] != 0)
        {
            // Write rho
            (*outFileStream) << *((*context) -> rho) << ",";        
        }

        if (variableList[2] != 0)
        {
            // Write gamma
            for (j = 0; j< nTpt; j++)
            {
                (*outFileStream) << ((*context) -> gamma)[j] << ","; 
            }
        }

        if (variableList[3] != 0)
        {
            // Write p_se
            for (i = 0; i < nLoc; i++)
            {
                for (j = 0; j < nTpt; j++)
                {
                    (*outFileStream) << ((*context)->p_se)[i*nTpt + j] <<  ",";
                }
            }
        }

        if (variableList[4] != 0)
        {
            // Write p_ei 
            (*outFileStream) << *((*context) -> p_ei) << ",";
        }

        if (variableList[5] != 0)
        {
            // Write p_ir 
            (*outFileStream) << *((*context) -> p_ir) << ",";

        }
        if (variableList[6] != 0)
        {
            // Write p_rs 
            for (j = 0; j < nTpt; j++)
            {
                (*outFileStream) << ((*context)->p_rs)[j] << ",";
            }

        }

        if (variableList[7] != 0)
        {
            // Write S* 
            for (i = 0; i < nTpt*nLoc; i++)
            {
                (*outFileStream) << ((*context) -> S_star -> data)[i] << ",";
            }
        }

        if (variableList[8] != 0)
        {
            // Write E* 
            for (i = 0; i < nTpt*nLoc; i++)
            {
                (*outFileStream) << ((*context) -> E_star -> data)[i] << ",";
            }
        }
        if (variableList[9] != 0)
        {
            // Write I* 
            for (i = 0; i < nTpt*nLoc; i++)
            {
                (*outFileStream) << ((*context) -> I_star -> data)[i] << ",";
            }
        }
        if (variableList[10] != 0)
        {
            // Write R* 
            for (i = 0; i < nTpt*nLoc; i++)
            {
                (*outFileStream) << ((*context) -> R_star -> data)[i] << ",";
            }
        }

        if (variableList[11] != 0)
        {
            // Write S total 
            (*outFileStream) << (*context) -> totalS() << ",";
        }
        if (variableList[12] != 0)
        {
            // Write E total 
            (*outFileStream) << (*context) -> totalE() <<",";
        }
        if (variableList[13] != 0)
        {
            // Write I total 
            (*outFileStream) << (*context) -> totalI() << ",";
        }
        if (variableList[14] != 0)
        {
            // Write R total 
            (*outFileStream) << (*context) -> totalR() << ",";
        }
        if (variableList[15] != 0)
        {
            // Write S_star total 
            (*outFileStream) << (*context) -> totalS_star() << ",";
        }
        if (variableList[16] != 0)
        {
            // Write E_star total 
            (*outFileStream) << (*context) -> totalE_star() << ",";
        }
        if (variableList[17] != 0)
        {
            // Write I_star total 
            (*outFileStream) << (*context) -> totalI_star() << ",";
        }
        if (variableList[18] != 0)
        {
            // Write R_star total 
            (*outFileStream) << (*context) -> totalR_star() << ",";
        }
        if (variableList[19] != 0)
        {
            // Write average p_se 
            (*outFileStream) << (*context) -> avgP_SE() << ",";
        }
        if (variableList[20] != 0)
        {
            // Write average p_se 
            (*outFileStream) << (*context) -> avgP_RS() << ",";
        }


    
        // Time specific

        if (variableList[21] != 0)
        {
            // Write S total 
            for (j = 0; j < nTpt; j++)
            {
                (*outFileStream) << (*context) -> totalS(j) << ",";
            }
        }
        if (variableList[22] != 0)
        {
            // Write E total 
            for (j = 0; j < nTpt; j++)
            {
                (*outFileStream) << (*context) -> totalE(j) << ",";
            }
        }
        if (variableList[23] != 0)
        {
            // Write I total 
            for (j = 0; j < nTpt; j++)
            {
                (*outFileStream) << (*context) -> totalI(j) << ",";
            }
        }
        if (variableList[24] != 0)
        {
            // Write R total 
            for (j = 0; j < nTpt; j++)
            {
                (*outFileStream) << (*context) -> totalR(j) << ",";
            }
        }
        if (variableList[25] != 0)
        {
            // Write S_star total 
            for (j = 0; j < nTpt; j++)
            {
                (*outFileStream) << (*context) -> totalS_star(j) << ",";
            }
        }
        if (variableList[26] != 0)
        {
            // Write E_star total 
            for (j = 0; j < nTpt; j++)
            {
                (*outFileStream) << (*context) -> totalE_star(j) << ",";
            }
        }
        if (variableList[27] != 0)
        {
            // Write I_star total             
            for (j = 0; j < nTpt; j++)
            {
                (*outFileStream) << (*context) -> totalI_star(j) << ",";
            }
        }
        if (variableList[28] != 0)
        {
            // Write R_star total 
            for (j = 0; j < nTpt; j++)
            {
                (*outFileStream) << (*context) -> totalR_star(j) << ",";
            }
        }
        if (variableList[29] != 0)
        {
            // Write average p_se 
            for (j = 0; j < nTpt; j++)
            {
                (*outFileStream) << (*context) -> avgP_SE(j) << ",";
            }
        }
        if (variableList[30] != 0)
        {
            // Write r_0 
            (*outFileStream) << (*context) -> estimateR0() << ", ";
        }
        if (variableList[31] != 0)
        {
            // Write r_0_t 
            for (j = 0; j < nTpt; j++)
            {
                (*outFileStream) << (*context) -> estimateR0(i) << ", ";
            }
        }



        (*outFileStream) << iteration << "," << difftime(time(&*timer), *startTime)  <<"\n";
        outFileStream -> flush();
        return(0);
    }
    int IOProvider::close()
    {
        outFileStream -> close();
        return(0);
    }
}

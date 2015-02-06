#include <IOProvider.hpp>
#include <ModelContext.hpp>
#include <CompartmentalModelMatrix.hpp>
#include <CovariateMatrix.hpp>
#include <LSS_FullConditionalList.hpp>
#include<time.h>


LSSCout lssCout;
namespace SpatialSEIR
{

    IOProvider::IOProvider()
    {
        // Empty constructor, do nothing. 
    }

    IOProvider::IOProvider(ModelContext* _context,
                           std::string* _outFilePath,
                           int* _iterationStride)
    {
        this -> populate(&*_context,&*_outFilePath,
                &*_iterationStride);
    }

    void IOProvider::setTrace(int locationIndex)     
    {
        if (!(*(this -> isOpen)))
        {
            lssCout << "Attempt to set trace for IOProvider without output file.\n";
            throw(-1);
        }
        if (locationIndex >= *((*context) -> S_star -> ncol))
        {
            lssCout << "Invalid Location Index\n";
            throw(-1);
        }

        std::cout << "Warning: you may be requesting a LOT of data.\n";
        int i;
        int nTpts = *((*context) -> S_star -> nrow);
        unsigned int j;
        bool alreadyTraced = false;
        for (j = 0; j < locationTraces -> size(); j++)
        {
            if (((*locationTraces)[j]) -> locationIndex == locationIndex)
            {
                alreadyTraced = true;
            }
        }
        if (!alreadyTraced)
        {
            LocationTrace* newTrace = new LocationTrace();
            (newTrace -> locationIndex) = locationIndex;
            (*locationTraces).push_back(newTrace); 
        }
        for (i = 0; i < nTpts; i++)
        {
            setTrace(locationIndex, i);
        }
    }

    void IOProvider::setTrace(int locationIndex,int timeIndex)     
    {
        if (!(*(this -> isOpen)))
        {
            std::cout << "Attempt to set trace for IOProvider without output file.\n";
            throw(-1);
        }
        if (locationIndex >= *((*context) -> S_star -> ncol))
        {
            lssCout << "Invalid Location Index\n";
            throw(-1);
        }
        if (timeIndex >= *((*context) -> S_star -> nrow))
        {
            lssCout << "Invalid Time Index\n";
            throw(-1);
        }


        TimeLocationTrace* newTrace = new TimeLocationTrace();
        (newTrace -> locationIndex) = locationIndex;
        (newTrace -> timeIndex) = timeIndex;
        timeLocationTraces -> push_back(newTrace);
        fileInit();
    }

    IOProvider::~IOProvider()
    {
        if (*isOpen)
        {
            this -> close();
            delete isOpen;
        }

        while ((timeLocationTraces -> size() != 0)){delete timeLocationTraces -> back(); timeLocationTraces -> pop_back();}
        delete timeLocationTraces;
        while ((locationTraces -> size() != 0)){delete locationTraces -> back(); locationTraces -> pop_back();}
        delete locationTraces;
        delete startTime;
        delete timer;
        delete iterationStride;
        delete outFilePath;
        delete outFileStream;
        delete context;
    }

    int IOProvider::populate(ModelContext* _context,
                           std::string* _outFilePath,
                           int* _iterationStride)
    {
        context = new ModelContext*;
        timeLocationTraces = new std::vector<TimeLocationTrace*>();
        locationTraces = new std::vector<LocationTrace*>();
        iterationStride = new int;
        outFilePath = new std::string; 
        *context = _context; 
        *outFilePath = *_outFilePath;
        *iterationStride = *_iterationStride;
        outFileStream = new std::ofstream;
        isOpen = new bool;
        *isOpen = false;
        timer = new time_t;
        startTime = new time_t;
        time(startTime); // Set start time
        return(this -> fileInit());
    }

    int IOProvider::fileInit()
    {
        unsigned int i;

        // Clear file 
        if (*isOpen)
        {
            this -> close();
            *isOpen = false;
        }
        FILE* tmp = fopen(outFilePath -> c_str(), "w");
        fclose(tmp);

        // Open file as output stream
        *isOpen = true;
        outFileStream -> open(outFilePath -> c_str());
    

        // Beta
        // Beta_P_RS
        // rho
        // gamma_ei
        // gamma_ir
        // Trace: for time j, loc i
        //      S0
        //      E0
        //      I0
        //      R0
        //      S_star
        //      E_star
        //      I_star
        //      R_star
        //      S
        //      E
        //      I
        //      R

        // New behavior: assume we want beta, beta_P_RS, rho, gamma_ei, gamma_ir, iteration, time
        // Additional data must be requested via setTrace. This should simplify the R function calls somewhat. 

        // Write Beta header
        unsigned int betaLen = (*((*context) -> X -> ncol_x)); 
        for (i = 0; i < betaLen; i++)
        {
            (*outFileStream) << "BetaP_SE_" << i << ", "; 
        }
        
        // Write Beta P_RS header
        betaLen = (*((*context) -> X_pRS -> ncol_x));

        if (((*context) -> config -> reinfectionMode) == 1)
        {
            for (i = 0; i < betaLen; i++)
            {
                (*outFileStream) << "BetaP_RS_" << i << ", "; 
            }
        }

        // Write rho header
        unsigned int nRho = ((*context) -> scaledDistMatrices -> size());
        if (*((*context) -> S_star -> ncol) > 1)
        {
            for (i = 0; i < nRho; i++)
            {
                (*outFileStream) << "rho_" << i << ",";       
            }
        }
        // Need to write overdispersion header?
        if (((*context) -> config -> dataModel) == 1)
        {
            (*outFileStream) << "phi,";
        }

        // Write gamma_ei header
        (*outFileStream) << "gamma_ei" << ",";

        // Write gamma_ir header
        (*outFileStream) << "gamma_ir" << ",";





        unsigned int nTpts = *((*context) -> S_star -> nrow);
        if (((*context) -> config -> reinfectionMode) == 1)
        {
            for (i = 0; i < nTpts; i++)
            {
                (*outFileStream)  << "P_RS_" << i << ", ";
            }
        }

        // Write header for any location and time-location traces 
        LocationTrace lTrace; 
        TimeLocationTrace tlTrace;
        for (i = 0; i < locationTraces -> size(); i++)
        {
            lTrace = (*((*locationTraces)[i]));
            (*outFileStream)  << "S0_" << lTrace.locationIndex << ", ";
            (*outFileStream)  << "E0_" << lTrace.locationIndex << ", ";
            (*outFileStream)  << "I0_" << lTrace.locationIndex << ", ";
            (*outFileStream)  << "R0_" << lTrace.locationIndex << ", ";
        }
        for (i = 0; i < timeLocationTraces -> size(); i++)
        {
           tlTrace = (*((*timeLocationTraces)[i])); 

           (*outFileStream)  << "S_" << tlTrace.locationIndex << "_" << tlTrace.timeIndex << ", ";
           (*outFileStream)  << "E_" << tlTrace.locationIndex << "_" << tlTrace.timeIndex << ", ";
           (*outFileStream)  << "I_" << tlTrace.locationIndex << "_" << tlTrace.timeIndex << ", ";
           (*outFileStream)  << "R_" << tlTrace.locationIndex << "_" << tlTrace.timeIndex << ", ";
           (*outFileStream)  << "S_star_" << tlTrace.locationIndex << "_" << tlTrace.timeIndex << ", ";
           (*outFileStream)  << "E_star_" << tlTrace.locationIndex << "_" << tlTrace.timeIndex << ", ";
           (*outFileStream)  << "I_star_" << tlTrace.locationIndex << "_" << tlTrace.timeIndex << ", ";
           (*outFileStream)  << "R_star_" << tlTrace.locationIndex << "_" << tlTrace.timeIndex << ", ";
           (*outFileStream)  << "P_SE_" << tlTrace.locationIndex << "_" << tlTrace.timeIndex << ", ";
        }

        (*outFileStream) << "Iteration,Time\n";
        // Write iteration number
        //Newline
        return(0);
    }
    int IOProvider::catIter(int iteration)
    {
        unsigned int i;
        // Don't require the user to have an output file
        if (!(*(this -> isOpen)))
        {
            return(-1);
        }
        //Check iteration stride
        if ((iteration % (*iterationStride)) != 0)
        {
            return(1);
        }

        // Write Beta
        int nTpt = *((*context) -> S_star -> nrow);
        unsigned int betaLen = (*((*context) -> X -> ncol_x)); 
        for (i = 0; i < betaLen; i++)
        {
            (*outFileStream) << ((*context) -> beta)[i] << ","; 
        }
        
        // Write Beta P_RS
        betaLen = (*((*context) -> X_pRS -> ncol_x));
        if (((*context) -> config -> reinfectionMode) == 1)
        {
            for (i = 0; i < betaLen; i++)
            {
                (*outFileStream) << ((*context) -> betaPrs)[i] << ","; 
            }
        }

        // Write rho
        unsigned int nRho = ((*context) -> scaledDistMatrices -> size());
        if (*((*context) -> S_star -> ncol) > 1)
        {
            for (i = 0; i < nRho; i++)
            {
                (*outFileStream) << ((*context) -> rho)[i] << ",";        
            }
        }
        
        // Need to write overdispersion header?
        if (((*context) -> config -> dataModel) == 1)
        {
            (*outFileStream) << *((*context) -> phi) << ",";
        }


        // Write gamma_ei
        (*outFileStream) << *((*context) -> gamma_ei) << ",";        

        // Write gamma_ir
        (*outFileStream) << *((*context) -> gamma_ir) << ",";        

        unsigned int nTpts = *((*context) -> S_star -> nrow);
        if (((*context) -> config -> reinfectionMode) == 1)
        {
            
            for (i = 0; i < nTpts; i++)
            {
               (*outFileStream)  << ((*context) -> p_rs)[i] << ", ";
            }
        }



        LocationTrace lTrace; 
        TimeLocationTrace tlTrace;

        for (i = 0; i < locationTraces -> size(); i++)
        {
            lTrace = (*(*locationTraces)[i]);

            (*outFileStream)  << ((*context) -> A0 -> S0)[lTrace.locationIndex] << ", ";
            (*outFileStream)  << ((*context) -> A0 -> E0)[lTrace.locationIndex] << ", ";
            (*outFileStream)  << ((*context) -> A0 -> I0)[lTrace.locationIndex] << ", ";
            (*outFileStream)  << ((*context) -> A0 -> R0)[lTrace.locationIndex] << ", ";
        }



        // Write any time-location traces 
       

        for (i = 0; i < timeLocationTraces -> size(); i++)
        {
           tlTrace = *((*timeLocationTraces)[i]);  
           (*outFileStream)  << ((*context) -> S -> data)[tlTrace.locationIndex*nTpt + tlTrace.timeIndex] << ", ";
           (*outFileStream)  << ((*context) -> E -> data)[tlTrace.locationIndex*nTpt + tlTrace.timeIndex] << ", ";
           (*outFileStream)  << ((*context) -> I -> data)[tlTrace.locationIndex*nTpt + tlTrace.timeIndex] << ", ";
           (*outFileStream)  << ((*context) -> R -> data)[tlTrace.locationIndex*nTpt + tlTrace.timeIndex] << ", ";
           (*outFileStream)  << ((*context) -> S_star -> data)[tlTrace.locationIndex*nTpt + tlTrace.timeIndex] << ", ";
           (*outFileStream)  << ((*context) -> E_star -> data)[tlTrace.locationIndex*nTpt + tlTrace.timeIndex] << ", ";
           (*outFileStream)  << ((*context) -> I_star -> data)[tlTrace.locationIndex*nTpt + tlTrace.timeIndex] << ", ";
           (*outFileStream)  << ((*context) -> R_star -> data)[tlTrace.locationIndex*nTpt + tlTrace.timeIndex] << ", ";
           (*outFileStream)  << ((*context) -> p_se)[tlTrace.locationIndex*nTpt + tlTrace.timeIndex] << ", ";
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

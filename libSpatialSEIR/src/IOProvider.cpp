#include <IOProvider.hpp>
#include <ModelContext.hpp>
#include <CompartmentalModelMatrix.hpp>
#include <CovariateMatrix.hpp>
#include <FullConditional.hpp>
#include<time.h>



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
            std::cout << "Attempt to set trace for IOProvider without output file.\n";
            throw(-1);
        }

        std::cout << "Warning: you may be requesting a LOT of data.\n";
        int i;
        int nTpts = *((*context) -> S_star -> nrow);
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
        delete[] timeLocationTraces;
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
        // p_ei
        // p_ir
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

        // New behavior: assume we want beta, beta_P_RS, rho, p_ei, p_ir, iteration, time
        // Additional data must be requested via setTrace. This should simplify the R function calls somewhat. 

        // Write Beta header
        unsigned int betaLen = (*((*context) -> X -> ncol_x)) + (*((*context) -> X -> ncol_z)); 
        for (i = 0; i < betaLen; i++)
        {
            (*outFileStream) << "BetaP_SE_" << i << ", "; 
        }
        
        // Write Beta P_RS header
        betaLen = (*((*context) -> X_pRS -> ncol_x));

        for (i = 0; i < betaLen; i++)
        {
            (*outFileStream) << "BetaP_RS_" << i << ", "; 
        }

        // Write rho header
        (*outFileStream) << "rho,";        

        // Write p_ei header
        (*outFileStream) << "p_ei" << ",";

        // Write p_ir header
        (*outFileStream) << "p_ir" << ",";

        // Write header for any time-location traces 
       
        TimeLocationTrace tlTrace;
        for (i = 0; i < timeLocationTraces -> size(); i++)
        {
           tlTrace = (*(*timeLocationTraces)[i]); 

           (*outFileStream)  << "S0_" << tlTrace.locationIndex << ", ";
           (*outFileStream)  << "E0_" << tlTrace.locationIndex << ", ";
           (*outFileStream)  << "I0_" << tlTrace.locationIndex << ", ";
           (*outFileStream)  << "R0_" << tlTrace.locationIndex << ", ";
           (*outFileStream)  << "S_" << tlTrace.locationIndex << "_" << tlTrace.timeIndex << ", ";
           (*outFileStream)  << "E_" << tlTrace.locationIndex << "_" << tlTrace.timeIndex << ", ";
           (*outFileStream)  << "I_" << tlTrace.locationIndex << "_" << tlTrace.timeIndex << ", ";
           (*outFileStream)  << "R_" << tlTrace.locationIndex << "_" << tlTrace.timeIndex << ", ";
           (*outFileStream)  << "S_star_" << tlTrace.locationIndex << "_" << tlTrace.timeIndex << ", ";
           (*outFileStream)  << "E_star_" << tlTrace.locationIndex << "_" << tlTrace.timeIndex << ", ";
           (*outFileStream)  << "I_star_" << tlTrace.locationIndex << "_" << tlTrace.timeIndex << ", ";
           (*outFileStream)  << "R_star_" << tlTrace.locationIndex << "_" << tlTrace.timeIndex << ", ";
           (*outFileStream)  << "P_SE_" << tlTrace.locationIndex << "_" << tlTrace.timeIndex << ", ";
           (*outFileStream)  << "P_RS_" << tlTrace.timeIndex << ", ";
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
        unsigned int betaLen = (*((*context) -> X -> ncol_x)) + (*((*context) -> X -> ncol_z)); 
        for (i = 0; i < betaLen; i++)
        {
            (*outFileStream) << ((*context) -> beta)[i] << ","; 
        }
        
        // Write Beta P_RS
        betaLen = (*((*context) -> X_pRS -> ncol_x));

        for (i = 0; i < betaLen; i++)
        {
            (*outFileStream) << ((*context) -> betaPrs)[i] << ","; 
        }

        // Write rho
        (*outFileStream) << *((*context) -> rho) << ",";        

        // Write p_ei
        (*outFileStream) << *((*context) -> p_ei) << ",";        

        // Write p_ir
        (*outFileStream) << *((*context) -> p_ir) << ",";        

        // Write any time-location traces 
       
        TimeLocationTrace tlTrace;
        for (i = 0; i < timeLocationTraces -> size(); i++)
        {
           tlTrace = *((*timeLocationTraces)[i]); 

        
           (*outFileStream)  << ((*context) -> A0 -> S0)[tlTrace.locationIndex] << ", ";
           (*outFileStream)  << ((*context) -> A0 -> E0)[tlTrace.locationIndex] << ", ";
           (*outFileStream)  << ((*context) -> A0 -> I0)[tlTrace.locationIndex] << ", ";
           (*outFileStream)  << ((*context) -> A0 -> R0)[tlTrace.locationIndex] << ", ";
           (*outFileStream)  << ((*context) -> S -> data)[tlTrace.locationIndex*nTpt + tlTrace.timeIndex] << ", ";
           (*outFileStream)  << ((*context) -> E -> data)[tlTrace.locationIndex*nTpt + tlTrace.timeIndex] << ", ";
           (*outFileStream)  << ((*context) -> I -> data)[tlTrace.locationIndex*nTpt + tlTrace.timeIndex] << ", ";
           (*outFileStream)  << ((*context) -> R -> data)[tlTrace.locationIndex*nTpt + tlTrace.timeIndex] << ", ";
           (*outFileStream)  << ((*context) -> S_star -> data)[tlTrace.locationIndex*nTpt + tlTrace.timeIndex] << ", ";
           (*outFileStream)  << ((*context) -> E_star -> data)[tlTrace.locationIndex*nTpt + tlTrace.timeIndex] << ", ";
           (*outFileStream)  << ((*context) -> I_star -> data)[tlTrace.locationIndex*nTpt + tlTrace.timeIndex] << ", ";
           (*outFileStream)  << ((*context) -> R_star -> data)[tlTrace.locationIndex*nTpt + tlTrace.timeIndex] << ", ";
           (*outFileStream)  << ((*context) -> p_se)[tlTrace.locationIndex*nTpt + tlTrace.timeIndex] << ", ";
           (*outFileStream)  << ((*context) -> p_rs)[tlTrace.timeIndex] << ", ";

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

#include <LSS_FullConditional.hpp>

#ifndef FULL_CONDITIONAL_S0_INC
#define FULL_CONDITIONAL_S0_INC

namespace SpatialSEIR
{
    class ModelContext;
    class CompartmentalModelMatrix;
    class CovariateMatrix;
    class OCLProvider;

    /**
     * FC_S0 gives the full conditional distribution for the vector of initially
     * susceptible individuals. 
     */
    class FC_S0 : public InitCompartmentFullConditional
    {
        public:

            FC_S0(ModelContext *_context,
                      CompartmentalModelMatrix *_S, 
                      CompartmentalModelMatrix *_E, 
                      CompartmentalModelMatrix *_E_star,
                      CompartmentalModelMatrix *_I_star, 
                      InitData *_A0,
                      double *_p_se,
                      double *_p_ei,
                      double sliceWidth);
            virtual ~FC_S0();
            virtual int evalCPU();
            virtual int evalCPU(int startLoc);
            virtual int evalOCL() ;
            virtual void sample(int verbose);
            virtual long double getValue();
            virtual void setValue(long double value);
            virtual int calculateRelevantCompartments();
            virtual int calculateRelevantCompartments_OCL();
            virtual int calculateRelevantCompartments(int startLoc);
            virtual void printDebugInfo(int loc);


            ModelContext** context;
            CompartmentalModelMatrix** S;
            CompartmentalModelMatrix** E;
            CompartmentalModelMatrix** E_star;
            CompartmentalModelMatrix** I_star;
            InitData** A0;
            double** p_se; 
            double** p_ei;
            long double* value;
    };
}

#endif

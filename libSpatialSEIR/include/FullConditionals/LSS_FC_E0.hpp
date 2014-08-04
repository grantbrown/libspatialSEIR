#include <LSS_FullConditional.hpp>
#ifndef FULL_CONDITIONAL_E0_INC
#define FULL_CONDITIONAL_E0_INC

namespace SpatialSEIR
{

    class ModelContext;
    class CompartmentalModelMatrix;
    class CovariateMatrix;
    class OCLProvider;

    /**
     * FC_E0 gives the full conditional distribution for the vector of initially
     * exposed individuals. 
     */
    class FC_E0 : public InitCompartmentFullConditional
    { 
        public:
            FC_E0(ModelContext *_context,
                      CompartmentalModelMatrix *_S, 
                      CompartmentalModelMatrix *_E, 
                      CompartmentalModelMatrix *_I, 
                      CompartmentalModelMatrix *_E_star,
                      CompartmentalModelMatrix *_I_star,
                      CompartmentalModelMatrix *_R_star, 
                      InitData *_A0,
                      double *_p_ir,
                      double *_p_ei,
                      double *_p_se,
                      double sliceWidth);
            virtual ~FC_E0(); 
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
            CompartmentalModelMatrix** I;
            CompartmentalModelMatrix** E_star;
            CompartmentalModelMatrix** I_star;
            CompartmentalModelMatrix** R_star;
            InitData** A0;
            double** p_se;
            double** p_ir;
            double** p_ei;
            long double* value;
    };
}

#endif

#include <LSS_FullConditional.hpp>
#ifndef FULL_CONDITIONAL_I0_INC
#define FULL_CONDITIONAL_I0_INC

namespace SpatialSEIR
{

    class ModelContext;
    class CompartmentalModelMatrix;
    class CovariateMatrix;
    class OCLProvider;

    /**
     *
     * FC_I0 gives the full conditional distribution for the vector of initially
     * infectious individuals. 
     */
    class FC_I0 : public InitCompartmentFullConditional
    {
        public:
            FC_I0(ModelContext *_context,
                      CompartmentalModelMatrix *_S, 
                      CompartmentalModelMatrix *_I, 
                      CompartmentalModelMatrix *_R, 
                      CompartmentalModelMatrix *_S_star,
                      CompartmentalModelMatrix *_E_star,
                      CompartmentalModelMatrix *_R_star,
                      InitData *_A0,
                      double *_p_ir,
                      double *_p_rs,
                      double *_p_se,
                      double sliceWidth);
            virtual ~FC_I0(); 
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
            CompartmentalModelMatrix** I;
            CompartmentalModelMatrix** R;
            CompartmentalModelMatrix** S_star;
            CompartmentalModelMatrix** E_star;
            CompartmentalModelMatrix** R_star;
            InitData** A0;
            double** p_ir;
            double** p_rs;
            double** p_se; 
            long double* value;
    };

}

#endif

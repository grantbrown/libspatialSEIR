#include <LSS_FullConditional.hpp>
#ifndef FULL_CONDITIONAL_R0_INC
#define FULL_CONDITIONAL_R0_INC

namespace SpatialSEIR
{
    using std::cout;
    using std::endl;

    class ModelContext;
    class CompartmentalModelMatrix;
    class CovariateMatrix;
    class OCLProvider;

    /**
     * FC_R0 gives the full conditional distribution for the vector of initially
     * removed/recovered individuals. 
     */
    class FC_R0 : public InitCompartmentFullConditional
    {
        public:
            FC_R0(ModelContext *_context,
                      CompartmentalModelMatrix *_R, 
                      CompartmentalModelMatrix *_S,
                      CompartmentalModelMatrix *_S_star, 
                      CompartmentalModelMatrix *_E_star, 
                      CompartmentalModelMatrix *_R_star,
                      InitData *_A0,
                      double *_p_rs,
                      double *_p_se,
                      double sliceWidth,
                      int _useOCL);
            virtual ~FC_R0(); 
            virtual int evalCPU();
            virtual int evalCPU(int startLoc);
            virtual int evalOCL() ;
            virtual void sample(int verbose);
            virtual int sampleCPU();
            virtual int sampleOCL();
            virtual long double getValue();
            virtual void setValue(long double value);
            virtual int calculateRelevantCompartments();
            virtual int calculateRelevantCompartments_OCL();
            virtual int calculateRelevantCompartments(int startLoc);
            virtual void printDebugInfo(int loc);

            ModelContext** context;
            CompartmentalModelMatrix** R;
            CompartmentalModelMatrix** S;
            CompartmentalModelMatrix** S_star;
            CompartmentalModelMatrix** E_star; 
            CompartmentalModelMatrix** R_star;
            InitData** A0;
            double** p_rs;
            double** p_se; 
            long double* value;
    };

}

#endif

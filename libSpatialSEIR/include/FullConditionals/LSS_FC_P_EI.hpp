#include <LSS_FullConditional.hpp>
#ifndef FULL_CONDITIONAL_P_EI_INC
#define FULL_CONDITIONAL_P_EI_INC

namespace SpatialSEIR
{
    using std::cout;
    using std::endl;

    class ModelContext;
    class CompartmentalModelMatrix;
    class CovariateMatrix;
    class OCLProvider;

    /**
     * FC_P_EI gives the full conditional distribution of p_ei, the 
     * probability that an exposed individual becomes infectious at a
     * given time point. 
     */
    class FC_P_EI : public ParameterFullConditional
    {
        public:
            FC_P_EI(ModelContext *_context,
                    CompartmentalModelMatrix *_I_star,
                    CompartmentalModelMatrix *_E,
                    InitData *_A0,
                    double *_p_ei,
                    double _priorAlpha,
                    double _priorBeta);
            virtual int evalCPU();
            virtual int evalOCL();
            virtual int sampleCPU();
            virtual int sampleOCL();
            virtual long double getValue();
            virtual void setValue(long double val);
            virtual int calculateRelevantCompartments();
            virtual int calculateRelevantCompartments_OCL();

            ~FC_P_EI();
            ModelContext** context;
            CompartmentalModelMatrix **I_star;
            CompartmentalModelMatrix **E;
            InitData **A0;
            double **p_ei;
            long double* value;
            double* priorAlpha;
            double* priorBeta;
    };



}

#endif

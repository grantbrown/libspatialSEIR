#include <LSS_FullConditional.hpp>
#ifndef FULL_CONDITIONAL_P_EI_INC
#define FULL_CONDITIONAL_P_EI_INC

namespace SpatialSEIR
{

    class ModelContext;
    class CompartmentalModelMatrix;
    class CovariateMatrix;
    class OCLProvider;

    /**
     * FC_Gamma_EI gives the full conditional distribution of gamma_ei, 
     * the parameter driving the probability that an exposed individual
     * becomes infectious at a given time point, p_ei. 
     */
    class FC_Gamma_EI : public ParameterFullConditional
    {
        public:
            FC_Gamma_EI(ModelContext *_context,
                    CompartmentalModelMatrix *_I_star,
                    CompartmentalModelMatrix *_E,
                    InitData *_A0,
                    double *_p_ei,
                    double *_gamma_ei,
                    double _priorAlpha,
                    double _priorBeta,
                    double sliceWidth);
            virtual int evalCPU();
            virtual int evalOCL();
            virtual void sample(int verbose);
            virtual long double getValue();
            virtual void setValue(long double val);
            virtual int calculateRelevantCompartments();
            virtual int calculateRelevantCompartments_OCL();

            ~FC_Gamma_EI();
            ModelContext** context;
            CompartmentalModelMatrix **I_star;
            CompartmentalModelMatrix **E;
            InitData **A0;
            double **p_ei;
            double **gamma_ei;
            long double* value;
            double* priorAlpha;
            double* priorBeta;
    };



}

#endif

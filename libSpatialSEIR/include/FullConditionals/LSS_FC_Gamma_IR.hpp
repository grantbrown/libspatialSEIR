#include <LSS_FullConditional.hpp>
#ifndef FULL_CONDITIONAL_P_IR_INC
#define FULL_CONDITIONAL_P_IR_INC

namespace SpatialSEIR
{

    class ModelContext;
    class CompartmentalModelMatrix;
    class CovariateMatrix;
    class OCLProvider;

    /**
     * FC_Gamma_IR gives the full conditional distribution of gamma_ir, the parameter driving the, 
     * probability that an infectious individual recovers/is removed at 
     * a given time point (p_ir).
     */
    class FC_Gamma_IR : public ParameterFullConditional
    {
        
        public:
            FC_Gamma_IR(ModelContext *_context,
                    CompartmentalModelMatrix *_R_star,
                    CompartmentalModelMatrix *_I, 
                    InitData *_A0,
                    double *_p_ir,
                    double *_gamma_ir,
                    double _priorAlpha,
                    double _priorBeta,
                    double sliceWidth);
            ~FC_Gamma_IR();
            virtual double evalPrior(); 
            virtual double* minimumValue();
            virtual double* maximumValue();
            virtual int evalCPU();
            virtual int evalOCL();
            virtual void sample(int verbose);
            virtual long double getValue();
            virtual void setValue(long double val);
            virtual int calculateRelevantCompartments();
            virtual int calculateRelevantCompartments_OCL();

            ModelContext **context;
            CompartmentalModelMatrix **R_star;
            CompartmentalModelMatrix **I;
            InitData **A0;
            double **p_ir;
            double **gamma_ir;
            long double* value;
            double* priorAlpha;
            double* priorBeta;
    };


}

#endif

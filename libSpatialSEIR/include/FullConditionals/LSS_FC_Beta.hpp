#include <LSS_FullConditional.hpp>

#ifndef FULL_CONDITIONAL_BETA_INC
#define FULL_CONDITIONAL_BETA_INC

namespace SpatialSEIR
{
    using std::cout;
    using std::endl;

    class ModelContext;
    class CompartmentalModelMatrix;
    class CovariateMatrix;
    class OCLProvider;

    /**
     * FC_Beta gives the full conditional distribution of beta, the 
     * vector of regression parameters capturing the exposure intensity 
     * process. 
     */
    class FC_Beta : public ParameterFullConditional
    {
        public:
            FC_Beta(ModelContext *_context,
                    CompartmentalModelMatrix *_E_star, 
                    CompartmentalModelMatrix *_S, 
                    InitData *_A0,
                    CovariateMatrix *_X,
                    double *_p_se, 
                    double *_beta,
                    double *_rho,
                    double sliceWidth,
                    double _priorPrecision); 
            ~FC_Beta();

            virtual int evalCPU();
            virtual int evalOCL();
            virtual void sample(int verbose);
            virtual long double getValue();
            virtual void setValue(long double val);
            virtual int calculateRelevantCompartments();
            virtual int calculateRelevantCompartments_OCL();

            ModelContext **context;
            CompartmentalModelMatrix **E_star; 
            CompartmentalModelMatrix **S; 
            InitData **A0;
            CovariateMatrix **X;
            double **p_se;
            double **beta;
            double **rho;
            long double* value;
            double* priorPrecision;
    };

}

#endif

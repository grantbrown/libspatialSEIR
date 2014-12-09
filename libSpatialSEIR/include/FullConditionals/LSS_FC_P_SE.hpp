#include <LSS_FullConditional.hpp>

#ifndef FULL_CONDITIONAL_P_SE_INC
#define FULL_CONDITIONAL_P_SE_INC

namespace SpatialSEIR
{

    class ModelContext;
    class CompartmentalModelMatrix;
    class CovariateMatrix;
    class OCLProvider;

    /**
     * FC_P_SE provides a joint full conditional for Beta and rho, the
     * vector of regression parameters for the exposure intensity 
     * and the spatial dependence parameters.  
     */
    class FC_P_SE : public ParameterFullConditional
    {
        public:
            FC_P_SE(ModelContext *_context,
                    CompartmentalModelMatrix *_E_star, 
                    CompartmentalModelMatrix *_S, 
                    InitData *_A0,
                    CovariateMatrix *_X,
                    double *_p_se, 
                    double *_beta,
                    double *_rho,
                    double sliceWidth,
                    double _priorPrecision,
                    double _priorRhoAlpha,
                    double _priorRhoBeta); 
            ~FC_P_SE();

            virtual double evalPrior();
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
            int* nBeta;
            int* nRho;
            double* combinedParams;
            double* priorRhoAlpha;
            double* priorRhoBeta;
    };

}

#endif

#include <LSS_FullConditional.hpp>

#ifndef FULL_CONDITIONAL_POSTERIOR
#define FULL_CONDITIONAL_POSTERIOR

namespace SpatialSEIR
{
    class ModelContext;
    class CompartmentalModelMatrix;
    class CovariateMatrix;
    class OCLProvider;

    /**
     * LSS_Posterior is an experimental class being used to facilitate direct posterior
     * approximation and several new MCMC sampling techniques. 
     */
    class LSS_Posterior : public FullConditional
    {
        public:
            LSS_Posterior(ModelContext * _context,
                                    int* _Y,
                                    CompartmentalModelMatrix *_S_star, 
                                    CompartmentalModelMatrix *_E_star, 
                                    CompartmentalModelMatrix *_I_star, 
                                    CompartmentalModelMatrix *_R_star, 
                                    CompartmentalModelMatrix *_S,
                                    CompartmentalModelMatrix *_E,
                                    CompartmentalModelMatrix *_I, 
                                    CompartmentalModelMatrix *_R, 
                                    double *_p_se,
                                    double *_p_ei,
                                    double *_p_ir,
                                    double *_p_rs,
                                    double *_gamma_ei,
                                    double *_gamma_ir,
                                    double *_phi,
                                    double _steadyStateConstraintPrecision);
            virtual int evalCPU();
            virtual int evalCPU(int i, int j);
            virtual int evalOCL();
            virtual void sample(int verbose);
            virtual long double getValue();
            virtual void setValue(long double val);
            virtual int calculateRelevantCompartments();
            virtual int calculateRelevantCompartments(int i, int j);
            virtual int calculateRelevantCompartments_OCL();
            virtual ~LSS_Posterior();

            ModelContext **context;
            int** Y;
            CompartmentalModelMatrix **S_star; 
            CompartmentalModelMatrix **E_star; 
            CompartmentalModelMatrix **I_star; 
            CompartmentalModelMatrix **R_star; 
            CompartmentalModelMatrix **S; 
            CompartmentalModelMatrix **E;
            CompartmentalModelMatrix **I;
            CompartmentalModelMatrix **R;
            InitData **A0;
            double **p_se;
            double **p_ei;
            double **p_ir;
            double **p_rs;
            double **phi;
            long double* value;
            double* steadyStateConstraintPrecision;
    };
}

#endif

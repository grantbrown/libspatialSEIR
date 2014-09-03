#include <LSS_FullConditional.hpp>
#ifndef FULL_CONDITIONAL_E_STAR_INC
#define FULL_CONDITIONAL_E_STAR_INC

namespace SpatialSEIR
{
    class ModelContext;
    class CompartmentalModelMatrix;
    class CovariateMatrix;
    class OCLProvider;

    /**
     * FC_E_Star gives the full conditional distribution of E_star, the 
     * transition matrix capturing individuals moving from the susceptible
     * category to the exposed category. 
     */
    class FC_E_Star : public CompartmentFullConditional
    {
        public:
            FC_E_Star(ModelContext *_context,
                      CompartmentalModelMatrix *_E_star, 
                      CompartmentalModelMatrix *_E, 
                      CompartmentalModelMatrix *_S, 
                      CompartmentalModelMatrix *_I_star,
                      CovariateMatrix *_X,
                      InitData *_A0,
                      double *_p_se,
                      double *_p_ei,
                      double *_rho,
                      double *_beta,
                      double _steadyStateConstraintPrecision,
                      double sliceWidth);
            ~FC_E_Star();

            virtual int evalCPU();
            virtual int evalCPU(int i, int j);
            virtual int evalOCL();
            virtual void sample(int verbose);
            virtual long double getValue();
            virtual void setValue(long double val);
            virtual int calculateRelevantCompartments();
            virtual int calculateRelevantCompartments(int i, int j);
            virtual int calculateRelevantCompartments_OCL();
            ModelContext **context;
            CompartmentalModelMatrix **E_star; 
            CompartmentalModelMatrix **E; 
            CompartmentalModelMatrix **S; 
            CompartmentalModelMatrix **I_star;
            CovariateMatrix **X;
            InitData **A0;
            double **p_se;
            double **p_ei;
            double **rho;
            double **beta;
            long double* value;
            double* steadyStateConstraintPrecision;
    };

}

#endif

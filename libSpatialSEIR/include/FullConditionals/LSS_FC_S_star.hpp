#include <LSS_FullConditional.hpp>

#ifndef FULL_CONDITIONAL_S_STAR_INC
#define FULL_CONDITIONAL_S_STAR_INC

namespace SpatialSEIR
{
    class ModelContext;
    class CompartmentalModelMatrix;
    class CovariateMatrix;
    class OCLProvider;

    /**
     * FC_S_Star gives the full conditional distribution of S_star, the 
     * transition matrix capturing individuals moving from the removed/recovered
     * category to the susceptible category. 
     */
    class FC_S_Star : public CompartmentFullConditional
    {
        public:
            FC_S_Star(ModelContext * _context,
                      CompartmentalModelMatrix *_S_star, 
                      CompartmentalModelMatrix *_S, 
                      CompartmentalModelMatrix *_R,
                      CompartmentalModelMatrix *_E_star,
                      CompartmentalModelMatrix *_R_star,
                      InitData *_A0,
                      CovariateMatrix *_X,
                      double *_p_se,
                      double *_p_rs,
                      double *_beta,
                      double *_rho,
                      double _steadyStateConstraintPrecision,
                      double sliceWidth);
            virtual int evalCPU();
            virtual int evalOCL();
            virtual void sample(int verbose);
            virtual long double getValue();
            virtual void setValue(long double val);
            virtual int calculateRelevantCompartments();
            virtual int calculateRelevantCompartments(int i, int j);
            virtual int calculateRelevantCompartments_OCL();
            virtual ~FC_S_Star();

            ModelContext **context;
            CompartmentalModelMatrix **S_star; 
            CompartmentalModelMatrix **S; 
            CompartmentalModelMatrix **R;
            CompartmentalModelMatrix **E_star;
            CompartmentalModelMatrix **R_star;
            InitData **A0;
            CovariateMatrix **X;
            double **p_se;
            double **p_rs;
            double **beta; 
            double **rho;
            long double* value;
            double* steadyStateConstraintPrecision;
    };
}

#endif

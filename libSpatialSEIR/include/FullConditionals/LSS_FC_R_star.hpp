#include <LSS_FullConditional.hpp>
#ifndef FULL_CONDITIONAL_R_STAR_INC
#define FULL_CONDITIONAL_R_STAR_INC

namespace SpatialSEIR
{

    class ModelContext;
    class CompartmentalModelMatrix;
    class CovariateMatrix;
    class OCLProvider;

    /**
     * FC_R_Star gives the full conditional distribution of R_star, the 
     * transition matrix capturing individuals moving from the infectious
     * category to the recovered/removed category. 
     */
    class FC_R_Star : public CompartmentFullConditional
    {
        public:
            FC_R_Star(ModelContext *_context,
                      CompartmentalModelMatrix *_R_star,
                      CompartmentalModelMatrix *_R,
                      CompartmentalModelMatrix *_I,
                      CompartmentalModelMatrix *_S_star,
                      CompartmentalModelMatrix *_E_star,
                      CompartmentalModelMatrix *_I_star,
                      CompartmentalModelMatrix *_S,
                      InitData *_A0,
                      double *_p_rs,
                      double *_p_ir,
                      double *_p_se,
                      double _steadyStateConstraintPrecision,
                      double sliceWidth);
            ~FC_R_Star();

            virtual int evalCPU();
            virtual int evalOCL();
            virtual void sample(int verbose);
            virtual long double getValue();
            virtual void setValue(long double val);
            virtual int calculateRelevantCompartments();
            virtual int calculateRelevantCompartments(int i, int j);
            virtual int calculateRelevantCompartments_OCL();
            ModelContext **context;
            CompartmentalModelMatrix **R_star;
            CompartmentalModelMatrix **R;
            CompartmentalModelMatrix **I;
            CompartmentalModelMatrix **S_star;
            CompartmentalModelMatrix **E_star;
            CompartmentalModelMatrix **I_star;
            CompartmentalModelMatrix **S;
            InitData **A0;
            double **p_rs;
            double **p_ir;
            double **p_se;
            long double* value;
            double* steadyStateConstraintPrecision;
    };

}

#endif

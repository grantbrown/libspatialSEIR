#include <LSS_FullConditional.hpp>
#ifndef FULL_CONDITIONAL_I_STAR_INC
#define FULL_CONDITIONAL_I_STAR_INC

namespace SpatialSEIR
{
    class ModelContext;
    class CompartmentalModelMatrix;
    class CovariateMatrix;
    class OCLProvider;

    /**
     * FC_I_Star gives the full conditional distribution of I_star, the 
     * transition matrix capturing individuals moving from the exposed
     * category to the infectious category. This particular class 
     * is only used for data models incorporating removal times. 
     */
    class FC_I_Star : public CompartmentFullConditional
    {
        public:
            FC_I_Star(ModelContext *_context,
                      CompartmentalModelMatrix *_E_star, 
                      CompartmentalModelMatrix *_I_star, 
                      CompartmentalModelMatrix *_R_star, 
                      CompartmentalModelMatrix *_E,
                      CompartmentalModelMatrix *_I,
                      CompartmentalModelMatrix *_S,
                      InitData *_A0,
                      double *_p_ei,
                      double *_p_ir,
                      double *_p_se,
                      double _steadyStateConstraintPrecision,
                      double sliceWidth);
            ~FC_I_Star();

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
            CompartmentalModelMatrix **I_star; 
            CompartmentalModelMatrix **R_star; 
            CompartmentalModelMatrix **E; 
            CompartmentalModelMatrix **I;
            CompartmentalModelMatrix **S;
            InitData **A0;
            double **p_ei;
            double **p_ir;
            double **p_se;
            long double* value;
            double* steadyStateConstraintPrecision;
    };

}

#endif

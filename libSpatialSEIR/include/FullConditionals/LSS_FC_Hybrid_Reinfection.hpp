#include <LSS_FullConditional.hpp>

#ifndef FULL_CONDITIONAL_HYBRID_REINFECT_INC
#define FULL_CONDITIONAL_HYBRID_REINFECT_INC

namespace SpatialSEIR
{
    using std::cout;
    using std::endl;

    class ModelContext;
    class CompartmentalModelMatrix;
    class CovariateMatrix;
    class OCLProvider;

    /**
     * FC_S_Star gives the full conditional distribution of S_star, the 
     * transition matrix capturing individuals moving from the removed/recovered
     * category to the susceptible category. 
     */
    class FC_Hybrid_Reinfection : public HybridFullConditional
    {
        public:

            FC_Hybrid_Reinfection(ModelContext * _context,
                                  CompartmentalModelMatrix *_S_star, 
                                  CompartmentalModelMatrix *_S, 
                                  CompartmentalModelMatrix *_R,
                                  CompartmentalModelMatrix *_E_star,
                                  CompartmentalModelMatrix *_R_star,
                                  CompartmentFullConditional *_S_star_FC,
                                  ParameterFullConditional *_beta_p_rs_fc,
                                  InitData *_A0,
                                  CovariateMatrix *_X,
                                  CovariateMatrix *_X_p_rs,
                                  double *_p_se,
                                  double *_p_rs,
                                  double *_beta_p_rs,
                                  double _tausq,
                                  double *_beta,
                                  double *_rho,
                                  double _steadyStateConstraintPrecision,
                                  int _useOCL);
            virtual int evalCPU();
            virtual int evalOCL();
            virtual void sample(int verbose);
            virtual int sampleCPU();
            virtual int sampleOCL();
            virtual long double getValue();
            virtual void setValue(long double val);
            virtual int calculateRelevantCompartments();
            virtual int calculateRelevantCompartments_OCL();
            virtual ~FC_Hybrid_Reinfection();

            ModelContext **context;
            CompartmentalModelMatrix **S_star; 
            CompartmentalModelMatrix **S; 
            CompartmentalModelMatrix **R;
            CompartmentalModelMatrix **E_star;
            CompartmentalModelMatrix **R_star;
            ParameterFullConditional **parameterFullConditional;
            CompartmentFullConditional **compartmentFullConditional;
            InitData **A0;
            CovariateMatrix **X;
            CovariateMatrix **X_p_rs;
            double **p_se;
            double **p_rs;
            double **beta_p_rs;
            double *tausq;
            double **beta; 
            double **rho;
            double *steadyStateConstraintPrecision;
            long double* value;
    };
}

#endif

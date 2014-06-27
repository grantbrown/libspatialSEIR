#include <LSS_FullConditional.hpp>


#ifndef FULL_CONDITIONAL_BETA_P_RS_INC
#define FULL_CONDITIONAL_BETA_P_RS_INC

namespace SpatialSEIR
{
    using std::cout;
    using std::endl;

    class ModelContext;
    class CompartmentalModelMatrix;
    class CovariateMatrix;
    class OCLProvider;

    /**
     * FC_Beta_P_RS gives the full conditional distribution of beta_p_rs, the 
     * vector of regression parameters capturing the probability that an individual
     * transitions from R to S, the removed/recovered category to the susceptible category. 
     */
    class FC_Beta_P_RS : public ParameterFullConditional
    {
        public:
            FC_Beta_P_RS(ModelContext *_context,
                    CompartmentalModelMatrix *_S_star,
                    CompartmentalModelMatrix *_R,
                    CovariateMatrix* _X,
                    InitData *_A0,
                    double *_p_rs,
                    double *_beta_p_rs,
                    double _tausq,
                    double _sliceWidth
                    );
            ~FC_Beta_P_RS();
            virtual int evalCPU();
            virtual int evalOCL();
            virtual int sampleCPU();
            virtual int sampleOCL();
            virtual long double getValue();
            virtual void setValue(long double val);
            virtual int calculateRelevantCompartments();
            virtual int calculateRelevantCompartments_OCL();

            ModelContext **context;
            CompartmentalModelMatrix **S_star;
            CompartmentalModelMatrix **R;
            CovariateMatrix **X;
            InitData **A0;
            double **beta_p_rs;
            double **p_rs;
            double* tausq;
            long double* value;
    };


}

#endif

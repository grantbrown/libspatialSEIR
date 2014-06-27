#include <LSS_FullConditional.hpp>
#ifndef FULL_CONDITIONAL_RHO_INC
#define FULL_CONDITIONAL_RHO_INC

namespace SpatialSEIR
{
    using std::cout;
    using std::endl;

    class ModelContext;
    class CompartmentalModelMatrix;
    class CovariateMatrix;
    class OCLProvider;

    /**
     * FC_Rho gives the full conditional distribution of rho, the 
     * scalar spatial dependence parameter. 
     */
    class FC_Rho : public ParameterFullConditional 
    {
        public:
            FC_Rho(ModelContext *_context,
                   CompartmentalModelMatrix *_E_star, 
                   CompartmentalModelMatrix *_S, 
                   InitData *_A0,
                   CovariateMatrix *_X,
                   double *_p_se, 
                   double *_beta, 
                   double *_rho,
                   double sliceWidth
                   );
            ~FC_Rho();
            virtual int evalCPU();
            virtual int evalOCL();
            virtual int sampleCPU();
            virtual int sampleOCL();
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
    };

}

#endif

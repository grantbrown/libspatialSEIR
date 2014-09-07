#include <LSS_FullConditional.hpp>
#ifndef FULL_CONDITIONAL_RHO_INC
#define FULL_CONDITIONAL_RHO_INC

namespace SpatialSEIR
{

    class ModelContext;
    class CompartmentalModelMatrix;
    class CovariateMatrix;
    class OCLProvider;

    /**
     * FC_Phi gives the full conditional distribution of rho, the 
     * scalar spatial dependence parameter. 
     */
    class FC_Phi : public ParameterFullConditional 
    {
        public:
            FC_Phi(ModelContext *_context,
                   CompartmentalModelMatrix *_I_star, 
                   double* _phi,
                   double _priorAlpha,
                   double _priorBeta,
                   int* _Y,
                   double sliceWidth
                   );
            ~FC_Phi();
            virtual double evalPrior();
            virtual int evalCPU();
            virtual int evalOCL();
            virtual void sample(int verbose);
            virtual long double getValue();
            virtual void setValue(long double val);
            virtual int calculateRelevantCompartments();
            virtual int calculateRelevantCompartments_OCL();

            ModelContext **context;
            CompartmentalModelMatrix **I_star; 
            double** phi;
            double* priorAlpha;
            double* priorBeta;
            int** Y;
            long double* value;
    };

}

#endif

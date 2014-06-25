#ifndef SPATIALSEIR_INCLUDEFILES
#define SPATIALSEIR_INCLUDEFILES

#include<iostream>
#include<stdio.h>
#include<cstring>
#include<vector>
#endif

#ifndef FULL_CONDITIONAL_GAMMA_INC
#define FULL_CONDITIONAL_GAMMA_INC

namespace SpatialSEIR
{
    using std::cout;
    using std::endl;

    class ModelContext;
    class CompartmentalModelMatrix;
    class CovariateMatrix;
    class OCLProvider;

    //! DEPRICATED - Prior arguments for the gamma term (external infection probability) 
    struct gammaArgs
    {
        double* priorAlpha;
        double* priorBeta;
        double* gamma;
    };

    /**
     * FC_Gamma is depricated. Gamma was used to describe a time varying external 
     * infection source, but was found to be unneccesary for most epidemic models, 
     * as well as a potential source for identifiability issues. The code remains 
     * in case this decision changes. 
     */
    class FC_Gamma : public ParameterFullConditional 
    {
        public:
            FC_Gamma(ModelContext *_context,
                   CompartmentalModelMatrix *_E_star, 
                   CompartmentalModelMatrix *_S, 
                   InitData *_A0,
                   CovariateMatrix *_X,
                   double *_p_se, 
                   double *_beta, 
                   double *_gamma,
                   double *_priorAlpha,
                   double *_priorBeta,
                   double sliceWidth
                   );
            ~FC_Gamma();
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
            double **gamma;
            double* priorAlpha;
            double* priorBeta;
            long double* value;
    };    

}

#endif

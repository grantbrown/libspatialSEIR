#include <LSS_FullConditional.hpp>

#ifndef FULL_CONDITIONAL_I_STAR_OD_INC
#define FULL_CONDITIONAL_I_STAR_OD_INC

namespace SpatialSEIR
{
    class ModelContext;
    class CompartmentalModelMatrix;
    class CovariateMatrix;
    class OCLProvider;

    /**
     * FC_I_Star_overdispersed gives the overdispersion data model for I_star. 
     */
    class FC_I_Star_overdispersed : public CompartmentFullConditional
    {
        public:
            FC_I_Star_overdispersed(ModelContext * _context,
                                    int* _Y,
                                    CompartmentalModelMatrix *_I_star, 
                                    CompartmentalModelMatrix *_S,
                                    CompartmentalModelMatrix *_I, 
                                    CompartmentalModelMatrix *_E,
                                    CompartmentalModelMatrix *_E_star,
                                    CompartmentalModelMatrix *_R_star,                                    
                                    double *_p_ei,
                                    double *_p_ir,
                                    double *_p_se,
                                    double *_phi,
                                    double _steadyStateConstraintPrecision,
                                    double _sliceWidth);
            virtual int evalCPU();
            virtual int evalCPU(int i, int j);
            virtual int evalOCL();
            virtual void sample(int verbose);
            virtual long double getValue();
            virtual void setValue(long double val);
            virtual int calculateRelevantCompartments();
            virtual int calculateRelevantCompartments(int i, int j);
            virtual int calculateRelevantCompartments_OCL();
            virtual ~FC_I_Star_overdispersed();

            ModelContext **context;
            int** Y;
            CompartmentalModelMatrix **I_star; 
            CompartmentalModelMatrix **E_star; 
            CompartmentalModelMatrix **I; 
            CompartmentalModelMatrix **E;
            CompartmentalModelMatrix **S;
            CompartmentalModelMatrix **R_star;
            InitData **A0;
            double **p_ei;
            double **p_ir;
            double **p_se;
            double **phi;
            long double* value;
            double* steadyStateConstraintPrecision;
    };
}

#endif

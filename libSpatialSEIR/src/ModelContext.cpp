#ifndef MODEL_CONTEXT_INC
#define MODEL_CONTEXT_INC
#include<ModelContext.hpp>
#endif


namespace SpatialSEIR
{
    using std::cout;
    using std::endl;

    ModelContext::ModelContext(CompartmentalModelMatrix *_S_star,
                               CompartmentalModelMatrix *_E_star,
                               CompartmentalModelMatrix *_I_star,
                               CompartmentalModelMatrix *_R_star,
                               CompartmentalModelMatrix *_S,
                               CompartmentalModelMatrix *_E,
                               CompartmentalModelMatrix *_I,
                               CompartmentalModelMatrix *_R, 
                               InitData *_A0,
                               CovariateMatrix *_X)
    {
        S_star = new CompartmentalModelMatrix*;
        E_star = new CompartmentalModelMatrix*;
        I_star = new CompartmentalModelMatrix*;
        R_star = new CompartmentalModelMatrix*;
        S = new CompartmentalModelMatrix*;
        E = new CompartmentalModelMatrix*;
        I = new CompartmentalModelMatrix*;
        R = new CompartmentalModelMatrix*;
        A0 = new InitData*;
        X = new CovariateMatrix*;
        
        *S_star = _S_star;
        *E_star = _E_star;
        *I_star = _I_star;
        *R_star = _R_star;
        *S = _S;
        *E = _E;
        *I = _I;
        *R = _R;
        *A0 = _A0;
        *X = _X;
    }
    ModelContext::~ModelContext()
    {
        delete S_star;
        delete E_star;
        delete I_star;
        delete R_star;
        delete S;
        delete E;
        delete I;
        delete R;
        delete A0;
        delete X;
    }
}



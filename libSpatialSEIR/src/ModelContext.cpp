#ifndef MODEL_CONTEXT_INC
#define MODEL_CONTEXT_INC
#include<ModelContext.hpp>
#endif


namespace SpatialSEIR
{
    using std::cout;
    using std::endl;

    ModelContext::ModelContext()
    {
        S_star = new CompartmentalModelMatrix();
        E_star = new CompartmentalModelMatrix();
        I_star = new CompartmentalModelMatrix();
        R_star = new CompartmentalModelMatrix();
        S = new CompartmentalModelMatrix();
        E = new CompartmentalModelMatrix();
        I = new CompartmentalModelMatrix();
        R = new CompartmentalModelMatrix();
        A0 = new InitData();
        X = new CovariateMatrix();
    }
    ModelContext::~ModelContext()
    {
        delete[] S_star;
        delete[] E_star;
        delete[] I_star;
        delete[] R_star;
        delete[] S;
        delete[] E;
        delete[] I;
        delete[] R;
        delete[] A0;
        delete[] X;
    }
}



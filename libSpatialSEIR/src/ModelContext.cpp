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
        rawDistMat = new DistanceMatrix();
        scaledDistMat = new DistanceMatrix();
    }

    // Method: calculateS
    // Accesses: A0, S_star, E_star
    // Updates: S
    void ModelContext::calculateS_CPU()
    {
        throw(-1);
    }
    void ModelContext::calculateS_OCL()
    {
        throw(-1);
    }

    // Method: calculateE
    // Accesses: A0, I_star, E_star, 
    // Updates: E
    void ModelContext::calculateE_CPU()
    {
        throw(-1);
    }
    void ModelContext::calculateE_OCL()
    {
        throw(-1);
    }

    // Method: calculateI
    // Accesses: A0, I_star, R_star
    // Updates: I
    void ModelContext::calculateI_CPU()
    {
        throw(-1);
    }
    void ModelContext::calculateI_OCL()
    {
        throw(-1);
    }

    // Method: calculateR
    // Accesses: A0, R_star, S_star
    // Updates: R
    void ModelContext::calculateR_CPU()
    {
        throw(-1);
    }
    void ModelContext::calculateR_OCL()
    {
        throw(-1);
    }

    void ModelContext::calculateP_SE_CPU()
    {
        throw(-1);
    }
    void ModelContext::calculateP_SE_OCL()
    {
        throw(-1);
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
        delete rawDistMat;
        delete scaledDistMat;
    }
}



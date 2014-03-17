#ifndef MODEL_CONTEXT_INC
#define MODEL_CONTEXT_INC
#include<ModelContext.hpp>
#endif

#ifndef BLAS_INC
#define BLAS_INC
#include<cblas.h> 
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
        // Load up S(t=1) from A0
        std::cout << "P_1\n";
        int i;
        int numLoc = *(A0 -> numLocations);
        std::cout << "P_2\n";
        int max2 = (*(S -> nrow))*((*(S -> ncol)));
    

        for (i = 0; i < numLoc; i++)
        {
            (S -> data)[i] = ((A0 -> S0)[i] + 
                    (A0 -> S_star0)[i] - 
                    (A0 -> E_star0)[i]);
        }
        std::cout << "P_3\n";

        for (i = numLoc; i < max2; i++)
        {
            (S -> data)[i] = ((S -> data)[i - numLoc] + 
                              (S_star -> data)[i] - 
                              (E_star -> data)[i]);
        }
        std::cout << "P_3\n";
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

    // Method: calculatePi
    // Accesses: beta, I, N, distMat, rho
    // Updates: p_se
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



#ifndef MODEL_CONTEXT_INC
#define MODEL_CONTEXT_INC
#include<ModelContext.hpp>
#endif

#ifndef BLAS_INC
#define BLAS_INC
#include<cblas.h> 
#endif

#include<cmath>


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
        N = new int; *N = -1;
        beta = new double; *beta = -1.0;
        eta = new double; *eta = -1.0;
    }

    void ModelContext::populate()
    { 
        delete N; delete beta; delete eta;
        N = new int[*(A0 -> numLocations)]; // Vector, could be changed to
                                            // a time varying matrix
        int nbeta = (*(X -> ncol_x) + (*(X -> ncol_z)));
        int neta = (*(X -> nrow_z));
        beta = new double[nbeta];
        eta = new double[neta];
        tmpContainer = new CompartmentalModelMatrix();
        tmpContainer -> createEmptyCompartment((S -> nrow), (S -> ncol));
        int i;
        for (i = 0; i < nbeta; i++)
        {
            beta[i] = 0.0;
        }
        for (i = 0; i < neta; i++)
        {
            eta[i] = 0.0;
        }
    }

    // Method: calculateS
    // Accesses: A0, S_star, E_star
    // Updates: S
    void ModelContext::calculateS_CPU()
    {
        calculateGenericCompartment_CPU(&*(this -> S), &*(this -> A0 -> S0),
                                    &*(this -> S_star), &*(this -> E_star),
                                    &*(this -> A0 -> S_star0), &*(this -> A0 -> E_star0));
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
        calculateGenericCompartment_CPU(&*(this -> E), &*(this -> A0 -> E0),
                                    &*(this -> E_star), &*(this -> I_star),
                                    &*(this -> A0 -> E_star0), &*(this -> A0 -> I_star0));
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
        calculateGenericCompartment_CPU(&*(this -> I), &*(this -> A0 -> I0),
                                    &*(this -> I_star), &*(this -> R_star),
                                    &*(this -> A0 -> I_star0), &*(this -> A0 -> R_star0));
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
        calculateGenericCompartment_CPU(&*(this -> R), &*(this -> A0 -> R0),
                                    &*(this -> R_star), &*(this -> S_star),
                                    &*(this -> A0 -> R_star0), &*(this -> A0 -> S_star0));
    }
    void ModelContext::calculateR_OCL()
    {
        throw(-1);
    }


    // Method: calculateGenericCompartment
    // Access: A0, compartments linked by compStar poitners
    // Updates: Compartment linked by comp pointer 
    void ModelContext::calculateGenericCompartment_CPU(CompartmentalModelMatrix *comp,int *comp0, 
                                                   CompartmentalModelMatrix *compStarAdd, 
                                                   CompartmentalModelMatrix *compStarSub, 
                                                   int *compStar0Add,int *compStar0Sub)
    {
        int i;
        int numLoc = *(A0 -> numLocations);
        int max2 = (*(comp -> nrow))*(*(comp -> ncol));
        for (i = 0; i < numLoc; i++)
        {
            (comp -> data)[i] = ((comp0)[i] + 
                    (compStar0Add)[i] - 
                    (compStar0Sub)[i]);
        }

        for (i = numLoc; i < max2; i++)
        {
            (comp -> data)[i] = (comp -> data)[i - numLoc] + 
                      (compStarAdd -> data)[i - numLoc] - 
                      (compStarSub -> data)[i - numLoc];
        }
       
    }
    
    void ModelContext::calculateGenericCompartment_OCL(int *comp,int *comp0, 
                                                   int *compStarAdd, int *compStarSub, 
                                                   int *compStar0Add,int *compStar0Sub)
    {
        throw(-1);
    }

    // Method: calculatePi
    // Accesses: beta, I, N, distMat, rho
    // Updates: p_se
    void ModelContext::calculateP_SE_CPU()
    {
        int i;
        //Update Eta
        this -> X -> calculate_eta_CPU(eta, beta);
        //Exponentiate
        int nrowz = *(X->nrow_z);
        for (i = 0; i < nrowz; i++)
        {
            eta[i] = std::exp(eta[i]);
        }
        // Calculate dmu: I/N * exp(eta)
        // Calculate rho*sqrt(idmat)
        // Calculate probs
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
        delete[] beta;
        delete[] eta;
    }
}



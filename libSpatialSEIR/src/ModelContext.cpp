#include <ModelContext.hpp>
#include <FullConditional.hpp>
#include <OCLProvider.hpp>
#include <CompartmentalModelMatrix.hpp>
#include <CovariateMatrix.hpp>
#include <DistanceMatrix.hpp>
#include <RandomNumberProvider.hpp>
#include <IOProvider.hpp>


#ifndef BLAS_INC
#define BLAS_INC
#include<cblas.h> 
#endif

#include<cmath>
#include<ctime>

namespace SpatialSEIR
{
    using std::cout;
    using std::endl;

    ModelContext::ModelContext()
    {
        random = new RandomNumberProvider(static_cast<unsigned int>(std::time(0)));
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
        rho = new double; *rho = 0.25;
        fileProvider = new IOProvider();
    }

    void ModelContext::populate(InitData* _A0,
                                covariateArgs* xArgs, 
                                compartmentArgs* S_starArgs,
                                compartmentArgs* E_starArgs,
                                compartmentArgs* I_starArgs,
                                compartmentArgs* R_starArgs,
                                distanceArgs* rawDistArgs,
                                scaledDistanceArgs* scaledDistArgs,
                                double* rho_, double* beta_, double* p_ei_, 
                                double* p_ir_, double* p_rs_, int* N_)
    {
        *A0 = *_A0;

        X -> genFromDataStream(xArgs -> inData_x, 
                               xArgs -> inData_z,
                               xArgs -> inRow_x,
                               xArgs -> inCol_x,
                               xArgs -> inRow_z,
                               xArgs -> inCol_z);
        S_star -> genFromDataStream(S_starArgs -> inData,
                                    S_starArgs -> inRow,
                                    S_starArgs -> inCol);
        E_star -> genFromDataStream(E_starArgs -> inData,
                                    E_starArgs -> inRow,
                                    E_starArgs -> inCol);
        I_star -> genFromDataStream(I_starArgs -> inData,
                                    I_starArgs -> inRow,
                                    I_starArgs -> inCol);
        R_star -> genFromDataStream(R_starArgs -> inData,
                                    R_starArgs -> inRow,
                                    R_starArgs -> inCol);

        S -> createEmptyCompartment(S_starArgs -> inRow,
                                    S_starArgs -> inCol);

        E -> createEmptyCompartment(S_starArgs -> inRow,
                                    S_starArgs -> inCol);

        I -> createEmptyCompartment(S_starArgs -> inRow,
                                    S_starArgs -> inCol);

        R -> createEmptyCompartment(S_starArgs -> inRow,
                                    S_starArgs -> inCol);

        rawDistMat -> genFromDataStream(rawDistArgs -> inData,
                                        rawDistArgs -> dim);
        scaledDistMat -> genFromDataStream(rawDistArgs -> inData,
                                           rawDistArgs -> dim);
        scaledDistMat -> scaledInvFunc_CPU(*(scaledDistArgs -> phi), 
                                           &*(scaledDistArgs -> inData));

        delete N; delete beta; delete eta;
        N = new int[(*(S -> nrow))*(*(S -> ncol))];                                          
        int nbeta = (*(X -> ncol_x) + (*(X -> ncol_z)));
        int neta = (*(X -> nrow_z));
        beta = new double[nbeta];
        eta = new double[neta];

        // Create empty compartment for calculation.
        tmpContainer = new CompartmentalModelMatrix();
        tmpContainer -> createEmptyCompartment((S -> nrow), (S -> ncol));

        // Initialize beta/eta
        int i;
        for (i = 0; i < nbeta; i++)
        {
            beta[i] = 0.0;
        }
        for (i = 0; i < neta; i++)
        {
            eta[i] = 0.0;
        }

        // Allocate space for the transition probabilities

        p_se = new double[*(S -> nrow)*(*(S->ncol))];
        p_se_components = new double[*(S -> nrow)*(*(S->ncol))];
        compartmentCache = new double[*(S -> nrow)*(*(S->ncol))];
        p_ei = new double;
        p_ir = new double;
        p_rs = new double[*(S->ncol)];
        rho = new double;

        // Wire up the full conditional classes
        S_star_fc = new FC_S_Star(this,
                                  S_star,
                                  S,
                                  R,
                                  E_star,
                                  A0,
                                  X,
                                  p_se,
                                  p_rs,
                                  beta,
                                  rho);
        E_star_fc = new FC_E_Star(this,
                                  E_star,
                                  E,
                                  S,
                                  I_star,
                                  X,A0,p_se,p_ei,
                                  rho,beta);
        R_star_fc = new FC_R_Star(this,
                                  R_star,
                                  R,
                                  I,
                                  S_star,
                                  A0,p_rs,p_ir);
        beta_fc = new FC_Beta(this,
                              E_star,
                              S_star,
                              A0,X,p_se,beta,rho);
        rho_fc = new FC_Rho(this,
                            E_star,
                            S,
                            A0,X,p_se,beta,rho);
        p_rs_fc = new FC_P_RS(this,S_star,R,A0,p_rs);
        p_ei_fc = new FC_P_EI(this,
                              I_star,
                              E,
                              A0,p_ei);
        p_ir_fc =  new FC_P_IR(this,
                             R_star,
                             I,
                             A0,p_ir);

        // Add initial data 

        *rho = *rho_;
        *p_ei = *p_ei_;
        *p_ir = *p_ir_;

        for (i = 0; i < (*(X -> ncol_x) + (*(X -> ncol_z))); i++)
        {
            beta[i] = beta_[i];
        }
        for (i = 0; i < *(S -> ncol); i++)
        {
            p_rs[i] = p_rs_[i];
        }
        for (i = 0; i< (*(S -> nrow))*(*(S->ncol)); i++)
        {
            N[i] = N_[i];
        } 

        // Calculate Compartments
        this -> calculateS_CPU();
        this -> calculateE_CPU();
        this -> calculateR_CPU();
        this -> calculateI_CPU();
        this -> calculateP_SE_CPU();

        // Initialize FC Values
        
        this -> beta_fc -> evalCPU();
        this -> p_rs_fc -> evalCPU();
        this -> p_ei_fc -> evalCPU();
        this -> p_ir_fc -> evalCPU();
        this -> rho_fc -> evalCPU();
        this -> S_star_fc -> evalCPU();
        this -> E_star_fc -> evalCPU();
        this -> R_star_fc -> evalCPU();
    }

    void ModelContext::simulationIter(int* useOCL, bool verbose = false)
    {
        std::cout << "Beta: " << beta_fc ->getValue() << "\n";
        std::cout << "p_rs: " << p_rs_fc ->getValue() << "\n";
        std::cout << "p_ei: " << p_ei_fc ->getValue() << "\n";
        std::cout << "p_ir: " << p_ir_fc ->getValue() << "\n";
        std::cout << "rho: " << rho_fc ->getValue() << "\n";
        std::cout << "S_star: " << S_star_fc ->getValue() << "\n";
        std::cout << "E_star: " << E_star_fc ->getValue() << "\n";
        std::cout << "R_star: " << R_star_fc ->getValue() << "\n";


        /*
        if (verbose){std::cout << "Sampling S_star\n";}
        if (useOCL[0] == 0){S_star_fc -> sampleCPU();}
        else {S_star_fc -> sampleOCL();}

        if (verbose){std::cout << "Sampling E_star\n";}
        if (useOCL[1] == 0){E_star_fc -> sampleCPU();}
        else {E_star_fc -> sampleOCL();}

        if (verbose){std::cout << "Sampling R_star\n";}
        if (useOCL[2] == 0){R_star_fc -> sampleCPU();}
        else {R_star_fc -> sampleOCL();}


        int i;
        int rowCol = (*(R->ncol))*(*(R->nrow));
        for (i = 0; i < rowCol;i++)
        {
            if ((S_star -> data)[i] > (R -> data)[i])
            {
                std::cout << "S_star too big: " << i << ", val:"<< S_star_fc -> getValue() << " \n";
                S_star_fc -> evalCPU();
                std::cout << "Value 2: " << S_star_fc -> getValue() << "\n"; 
                break;
            }
            if ((S_star -> data)[i] < 0)
            {
                std::cout << "S_star <0: " << i << ", val:"<< S_star_fc -> getValue() << " \n";
                break;
            }
        }
        for (i = 0; i < rowCol;i++)
        {
            if ((E_star -> data)[i] > (S -> data)[i])
            {
                std::cout << "E_star too big: " << i << ", val:"<< E_star_fc -> getValue() << " \n";
                break;
            }
            if ((E_star -> data)[i] < 0)
            {
                std::cout << "E_star <0: " << i << ", val:"<< E_star_fc -> getValue() << " \n";
                break;
            }
        }
        for (i = 0; i < rowCol;i++)
        {
            if ((I_star -> data)[i] > (E -> data)[i])
            {
                std::cout << "I_star too big: " << i << "\n";
                break;
            }
            if ((I_star -> data)[i] < 0)
            {
                std::cout << "I_star <0: " << i << ", val: \n";
                break;
            }
        }

        for (i = 0; i < rowCol;i++)
        {
            if ((R_star -> data)[i] > (I -> data)[i])
            {
                
                std::cout << "R_star too big: " << i << ", val:"<< R_star_fc -> getValue() << " \n";
                R_star_fc -> evalCPU();
                std::cout << "Value 2: " << R_star_fc -> getValue() << "\n"; 
                break;
            }
            if ((R_star -> data)[i] < 0)
            {
                std::cout << "R_star <0: " << i << ", val:"<< R_star_fc -> getValue() << " \n";
                break;
            }

        }
        */
        if (verbose){std::cout << "Sampling rho\n";}
        if (useOCL[7] == 0){rho_fc -> sampleCPU();}
        else {rho_fc -> sampleOCL();}
        /*
        if (verbose){std::cout << "Sampling beta\n";}
        if (useOCL[3] == 0){beta_fc -> sampleCPU();}
        else {beta_fc -> sampleOCL();}


        if (verbose){std::cout << "Sampling p_rs\n";}
        if (useOCL[4] == 0){p_rs_fc -> sampleCPU();}
        else {p_rs_fc -> sampleOCL();}
        */
        /*
        if (verbose){std::cout << "Sampling p_ei\n";}
        if (useOCL[5] == 0){p_ei_fc -> sampleCPU();}
        else {p_ei_fc -> sampleOCL();}

        if (verbose){std::cout << "Sampling p_ir\n";}
        if (useOCL[6] == 0){p_ir_fc -> sampleCPU();}
        else {p_ir_fc -> sampleOCL();}
        */

    }


    // Method: runSimulation
    // Accesses: Everything lol
    // Updates: Everything lol
    void ModelContext::runSimulation(int nIterations, int* useOCL, bool verbose = false)
    {
        int i;
        for (i = 0; i < nIterations; i++)
        {
            this -> simulationIter(&*useOCL, verbose);
            this -> fileProvider -> catIter(i);
        }

    }

    void ModelContext::runSimulation_CPU(int nIterations, bool verbose = false)
    {
        int i;
        int useOCL[8] = {0};
        for (i = 0; i < nIterations; i++)
        {
            this -> simulationIter(&*useOCL, verbose);
            this -> fileProvider -> catIter(i);
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
    void ModelContext::calculateS_CPU(int startLoc, int startTime)
    {
        calculateGenericCompartment_CPU(&*(this -> S), &*(this -> A0 -> S0),
                                    &*(this -> S_star), &*(this -> E_star),
                                    &*(this -> A0 -> S_star0), &*(this -> A0 -> E_star0),
                                    startLoc, startTime);
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

    void ModelContext::calculateE_CPU(int startLoc, int startTime)
    {
        calculateGenericCompartment_CPU(&*(this -> E), &*(this -> A0 -> E0),
                                    &*(this -> E_star), &*(this -> I_star),
                                    &*(this -> A0 -> E_star0), &*(this -> A0 -> I_star0),
                                    startLoc, startTime);
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

        int i;
        int maxItr = (*(I -> nrow))*(*(I -> ncol));
        for (i = 0; i < maxItr; i++)
        {
            (I->data)[i] = N[i] - (S->data)[i] - (E->data)[i] - (R->data)[i]; 
        }


        /*
        calculateGenericCompartment_CPU(&*(this -> I), &*(this -> A0 -> I0),
                                    &*(this -> I_star), &*(this -> R_star),
                                    &*(this -> A0 -> I_star0), &*(this -> A0 -> R_star0));
                                    */
    }

    void ModelContext::calculateI_CPU(int startLoc, int startTime)
    {
 

        int i,j,startIdx,idx;
        startIdx = startTime*(*(R->nrow)) + startLoc;

        j = 0;
        for (i = startTime; i < *(R->ncol); i++)
        {
            idx = startIdx + j*(*(R->nrow));
            (I -> data)[idx] = N[idx] - (S->data)[idx] - (E->data)[idx] - (R->data)[idx];  
            j += 1;
        }

        /*
        calculateGenericCompartment_CPU(&*(this -> I), &*(this -> A0 -> I0),
                                    &*(this -> I_star), &*(this -> R_star),
                                    &*(this -> A0 -> I_star0), &*(this -> A0 -> R_star0),
                                    startLoc, startTime);
            */
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

    void ModelContext::calculateR_CPU(int startLoc, int startTime)
    {

        calculateGenericCompartment_CPU(&*(this -> R), &*(this -> A0 -> R0),
                                    &*(this -> R_star), &*(this -> S_star),
                                    &*(this -> A0 -> R_star0), &*(this -> A0 -> S_star0),
                                    startLoc, startTime);

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

    void ModelContext::calculateGenericCompartment_CPU(CompartmentalModelMatrix *comp,int *comp0, 
                                                   CompartmentalModelMatrix *compStarAdd, 
                                                   CompartmentalModelMatrix *compStarSub, 
                                                   int *compStar0Add,int *compStar0Sub,
                                                   int startLoc, int startTime)
    {
        int i;
        int numLoc = *(A0 -> numLocations);
        int numTpts = *(comp -> ncol);
        int idx;

        if (startTime == 0)
        {
            idx = startLoc;
            (comp -> data)[idx] = ((comp0)[idx] + 
                    (compStar0Add)[idx] - 
                    (compStar0Sub)[idx]);
        }

        startTime = (startTime == 0 ? 1 : startTime);

        for (i = startTime; i < numTpts; i++)
        {
            idx = startLoc + i*numLoc;
            (comp -> data)[idx] = (comp -> data)[idx - numLoc] + 
                      (compStarAdd -> data)[idx - numLoc] - 
                      (compStarSub -> data)[idx - numLoc];
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
        this -> cacheP_SE_Calculation(); 
        int i, j, index;

        // Calculate dmu: I/N * exp(eta)
        int nLoc = *(S -> nrow);
        int nCol = *(S -> ncol);
        for (j = 0; j < nCol; j++)
        {
            for (i = 0; i < nLoc; i++) 
            {
                index = i + j*nLoc;
                p_se[index] = 0.0;
                p_se_components[index] = 
                   ((I -> data)[index] * (eta[index]))/N[index];
            }
        }

        // Calculate rho*sqrt(idmat)
        SpatialSEIR::matMult(this -> p_se, 
                scaledDistMat -> data, 
                p_se_components, 
                *(scaledDistMat -> numLocations), 
                *(scaledDistMat -> numLocations),
                *(I -> nrow),
                *(I -> ncol),false,false);
        for (j = 0; j < nCol; j++)
        {
            for (i = 0; i < nLoc; i++) 
            {
                index = i + j*nLoc;
                p_se[index] = 1-exp(-p_se_components[index] - (*rho)*p_se[index]);
            }
        }        
    }

    // To be used when beta is fixed, eta has already been exponentiated,
    // only change is to I compartment. 
    void ModelContext::calculateP_SE_CPU(int startLoc, int startTime)
    {
        int i, j, index;
        // Calculate dmu: I/N * exp(eta)
        int nLoc = *(S -> nrow);
        int nCol = *(S -> ncol);

        i = startLoc;
        for (j = startTime; j < nCol; j++)
        {
            index = i + j*nLoc;
            p_se_components[index] = 
               ((I -> data)[index] * (eta[index]))/N[i];
        }

        // Calculate rho*sqrt(idmat)
        // Full version:
        /*
        SpatialSEIR::matMult(this -> p_se, 
                scaledDistMat -> data, 
                p_se_components, 
                *(scaledDistMat -> numLocations), 
                *(scaledDistMat -> numLocations),
                *(I -> nrow),
                *(I -> ncol),false,false);
        */ 
        // Start at current time
        int numTimeUpdate = (nCol - startTime);
        SpatialSEIR::matMult(&((this -> p_se)[startTime*nLoc]), 
                scaledDistMat -> data, 
                &(p_se_components[startTime*nLoc]), 
                *(scaledDistMat -> numLocations), 
                *(scaledDistMat -> numLocations),
                *(I -> nrow),
                (nCol - startTime),false,false);

        for (j = startTime; j < nCol; j++)
        {
            for (i = 0; i < nLoc; i++) 
            {
                index = i + j*nLoc;
                p_se[index] = 1-exp(-p_se_components[index] - (*rho)*p_se[index]);
            }
        } 
    }

    void ModelContext::cacheP_SE_Calculation()
    {
        int i, j, index;
        //Update Eta
        this -> X -> calculate_eta_CPU(eta, beta);

        //Exponentiate
        int nrowz = *(X->nrow_z);
        for (i = 0; i < nrowz; i++)
        {
            eta[i] = std::exp(eta[i]);
        }
        // Calculate dmu: I/N * exp(eta)
        int nLoc = *(S -> nrow);
        int nCol = *(S -> ncol);
        for (j = 0; j < nCol; j++)
        {
            for (i = 0; i < nLoc; i++) 
            {
                index = i + j*nLoc;
                p_se_components[index] = 
                   ((I -> data)[index] * (eta[index]))/N[i];
            }
        }
        /*It would be good to cache parts of this matrix multiplication...*/
        /*
        SpatialSEIR::matMult(this -> p_se, 
                scaledDistMat -> data, 
                scratch, 
                *(scaledDistMat -> numLocations), 
                *(scaledDistMat -> numLocations),
                *(I -> nrow),
                *(I -> ncol),false,false);
        for (j = 0; j < nCol; j++)
        {
            for (i = 0; i < nLoc; i++) 
            {
                index = i + j*nLoc;
                p_se[index] = 1-exp(-scratch[index] - *rho*p_se[index]);
            }
        } 
        */
    }

    void ModelContext::calculateP_SE_OCL()
    {
        throw(-1);
    }

    ModelContext::~ModelContext()
    {
        delete fileProvider;
        delete random;
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
        delete tmpContainer;
        delete rawDistMat;
        delete scaledDistMat;
        delete[] beta;
        delete[] eta;
        delete[] p_se;
        delete[] p_se_components;
        delete[] compartmentCache;
        delete p_ei;
        delete p_ir;
        delete[] p_rs;
        delete rho;
    }
}



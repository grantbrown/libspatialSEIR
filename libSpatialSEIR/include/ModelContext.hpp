#ifndef SPATIALSEIR_INCLUDEFILES
#define SPATIALSEIR_INCLUDEFILES

#include<iostream>
#include<stdio.h>
#include<cstring>
#include<vector>
#endif

#ifndef MODEL_CONTEXT_INC
#define MODEL_CONTEXT_INC

namespace SpatialSEIR
{
    using std::cout;
    using std::endl;

    class FC_S_Star;
    class FC_E_Star;
    class FC_R_Star;
    class FC_Beta;
    class FC_P_EI;
    class FC_P_IR;
    class FC_P_RS;
    class FC_Rho;
    class FC_Gamma;
    class InitData;

    class CompartmentalModelMatrix;
    class CovariateMatrix;
    class DistanceMatrix;
    class RandomNumberProvider;
    class OCLProvider;
    class IOProvider;
    struct covariateArgs;
    struct compartmentArgs;
    struct distanceArgs;
    struct scaledDistanceArgs;
    struct gammaArgs;
    struct sliceParameters;

    class ModelContext
    {
        public:
            //Methods
            ModelContext();
            ~ModelContext(); 

            void populate(InitData* _A0,
                          covariateArgs* xArgs,
                          compartmentArgs* S_starArgs, 
                          compartmentArgs* E_starArgs,
                          compartmentArgs* I_starArgs,
                          compartmentArgs* R_starArgs,
                          distanceArgs* rawDistArgs,
                          scaledDistanceArgs* scaledDistArgs,
                          gammaArgs* gammaFCArgs,
                          double* rho, double* beta, 
                          double* p_ei, double* p_ir, double* p_rs, 
                          int* N, sliceParameters* sliceWidths);


            // Method: calculateS
            // Accesses: A0, S_star, E_star
            // Updates: S
            void calculateS_CPU();
            void calculateS_CPU(int startLoc, int startTime);
            void calculateS_givenE_CPU();
            void calculateS_givenE_CPU(int startLoc, int startTime);
            void calculateS_OCL();

            // Method: calculateE
            // Accesses: A0, I_star, E_star
            // Updates: E
            void calculateE_CPU();
            void calculateE_CPU(int startLoc, int startTime);
            void calculateE_givenI_CPU();
            void calculateE_givenI_CPU(int startLoc, int startTime);
            void calculateE_OCL();

            // Method: calculateI
            // Accesses: A0, I_star, R_star
            // Updates: I
            void calculateI_CPU();
            void calculateI_CPU(int startLoc, int startTime);
            void calculateI_givenR_CPU();
            void calculateI_givenR_CPU(int startLoc, int startTime);
            void calculateI_OCL();

            // Method: calculateR
            // Accesses: A0, R_star, S_star
            // Updates: R
            void calculateR_CPU();
            void calculateR_CPU(int startLoc, int startTime);
            void calculateR_givenS_CPU();
            void calculateR_givenS_CPU(int startLoc, int startTime);
            void calculateR_OCL();

            // Method calculateGenericCompartment
            // Accesses: A0, compartments linked by compStar pointers
            // Updates: Compartment linked by comp pointer
            void calculateGenericCompartment_CPU(CompartmentalModelMatrix *comp, int *comp0,
                                                 CompartmentalModelMatrix *compStarAdd, 
                                                 CompartmentalModelMatrix *compStarSub,
                                                 int *compStar0Add, int *compStar0Sub);

            void calculateGenericCompartment_CPU(CompartmentalModelMatrix *comp, int *comp0,
                                                 CompartmentalModelMatrix *compStarAdd, 
                                                 CompartmentalModelMatrix *compStarSub,
                                                 int *compStar0Add, int *compStar0Sub,
                                                 int startLoc, int startTime);

            void calculateGenericCompartment_OCL(int *comp, int *comp0,
                                                 int *compStarAdd, int *compStarSub,
                                                 int *compStar0Add, int *compStar0Sub);

            // Run main simulation 
            int checkCompartmentBounds();
            void printFCValues();
            void setRandomSeed(unsigned int seedValue);
            void simulationIter(int* useOCL, bool verbose, bool debug);
            void runSimulation(int nIterations, int* useOCL, bool verbose, bool debug);
            void runSimulation_CPU(int nIterations, bool verbose, bool debug);


            // Method: calculatePi
            // Accesses: beta, I, N, distMat, rho
            // Updates: p_se
            void cacheP_SE_Calculation();
            void calculateP_SE_CPU();
            void calculateP_SE_CPU(int startLoc, int startTime);
            void calculateP_SE_OCL();
        
            //Logic provider and utility classes
            IOProvider *fileProvider;
            RandomNumberProvider *random;
            OCLProvider *oclProvider; 
            FC_S_Star *S_star_fc;
            FC_E_Star *E_star_fc;
            FC_R_Star *R_star_fc;
            FC_Beta *beta_fc;
            FC_Rho *rho_fc;
            FC_Gamma *gamma_fc;
            FC_P_RS *p_rs_fc;
            FC_P_EI *p_ei_fc;
            FC_P_IR *p_ir_fc;

            //Data
            CompartmentalModelMatrix* S;
            CompartmentalModelMatrix* E;
            CompartmentalModelMatrix* I;
            CompartmentalModelMatrix* R;
            CompartmentalModelMatrix* S_star;
            CompartmentalModelMatrix* E_star;
            CompartmentalModelMatrix* I_star;
            CompartmentalModelMatrix* R_star;
            InitData* A0;
            CovariateMatrix* X;
            DistanceMatrix* rawDistMat;
            DistanceMatrix* scaledDistMat;
            CompartmentalModelMatrix* tmpContainer;

            double* beta;
            double* rho;
            double* gamma;
            double* eta;
            double* p_se;
            double* p_se_components;
            double* compartmentCache;
            double* p_ei;
            double* p_ir;
            double* p_rs;
            int* N;// Matrix?
    };
}
#endif

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

    class FC_S0;
    class FC_E0;
    class FC_I0;
    class FC_R0;
    class FC_S_Star;
    class FC_E_Star;
    class FC_R_Star;
    class FC_Beta;
    class FC_P_EI;
    class FC_P_IR;
    class FC_Beta_P_RS;
    class FC_Rho;
    class InitData;
    class TestClass;

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
    struct priorControl;
    struct modelConfiguration
    {
        int reinfectionMode;
        int samplingMode; 
    };



    class ModelContext
    {
        public:
            //Methods
            ModelContext();
            ~ModelContext(); 

            void populate(InitData* _A0,
                          covariateArgs* xArgs,
                          covariateArgs* xPrsArgs,
                          compartmentArgs* S_starArgs, 
                          compartmentArgs* E_starArgs,
                          compartmentArgs* I_starArgs,
                          compartmentArgs* R_starArgs,
                          distanceArgs* rawDistArgs,
                          scaledDistanceArgs* scaledDistArgs,
                          double* rho, double* beta, 
                          double* p_ei, double* p_ir, double* betaPrs, 
                          int* N, sliceParameters* sliceWidths,
                          priorControl* priorInformation,
                          modelConfiguration _config);

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
                                                 CompartmentalModelMatrix *compStarSub); 

            void calculateGenericCompartment_CPU(CompartmentalModelMatrix *comp, int *comp0,
                                                 CompartmentalModelMatrix *compStarAdd, 
                                                 CompartmentalModelMatrix *compStarSub,
                                                 int startLoc, int startTime);

            void calculateGenericCompartment_OCL(int *comp, int *comp0,
                                                 int *compStarAdd, int *compStarSub);

            // Run main simulation 
            int checkCompartmentBounds();
            void printFCValues();
            void setRandomSeed(unsigned int seedValue);
            void simulationIter(bool verbose, bool debug);
            void runSimulation(int nIterations,  bool verbose, bool debug);
            void updateSamplingParameters(double desiredRatio, double targetWidth, double proportionChange);


            // Method: calculatePi
            // Accesses: beta, I, N, distMat, rho
            // Updates: p_se
            void cacheP_SE_Calculation();
            void calculateP_SE_CPU();
            void calculateP_SE_CPU(int startLoc, int startTime);
            void calculateP_SE_OCL();

            // Method: calculateP_RS
            // Accesses: betaPrs, R,S,S_star,R_star
            // Updates: p_rs
            void calculateP_RS_CPU();

            // Summary Functions:
            int totalS();
            int totalE();
            int totalI();
            int totalR();

            int totalS(int tpt);
            int totalE(int tpt);
            int totalI(int tpt);
            int totalR(int tpt);

            int totalS_star();
            int totalE_star();
            int totalI_star();
            int totalR_star();

            int totalS_star(int tpt);
            int totalE_star(int tpt);
            int totalI_star(int tpt);
            int totalR_star(int tpt);

            double avgP_SE();
            double avgP_RS();

            double avgP_SE(int tpt);

            double estimateR0();
            double estimateR0(int t);
            double* calculateG(int t);


            //Logic provider and utility classes
            IOProvider *fileProvider;
            RandomNumberProvider *random;
            OCLProvider *oclProvider; 
            FC_S0 *S0_fc;
            FC_E0 *E0_fc;
            FC_I0 *I0_fc;
            FC_R0 *R0_fc;
            FC_S_Star *S_star_fc;
            FC_E_Star *E_star_fc;
            FC_R_Star *R_star_fc;
            FC_Beta *beta_fc;
            FC_Rho *rho_fc;
            FC_Beta_P_RS *betaPrs_fc;
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
            CovariateMatrix* X_pRS;
            DistanceMatrix* rawDistMat;
            DistanceMatrix* scaledDistMat;
            CompartmentalModelMatrix* tmpContainer;
            modelConfiguration* config;

            double* beta;
            double* betaPrs;
            double* rho;
            double* gamma;
            double* eta;
            double* p_se;
            double* p_se_components;
            double* compartmentCache;
            double* p_ei;
            double* p_ir;
            double* p_rs;
            int* N;
            int* isPopulated;
            int* singleLocation;
            int* numIterations;

            int* S0_OCL;
            int* I0_OCL;
            int* S_star_OCL;
            int* E_star_OCL;
            int* R_star_OCL;
            int* rho_OCL;
            int* beta_OCL;
            int* beta_P_RS_OCL;
    };
}
#endif

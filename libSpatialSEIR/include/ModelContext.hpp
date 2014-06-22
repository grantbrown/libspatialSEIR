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

    /**
     * The modelConfiguration struct stores information on the reinfection mode requested for the model
     * as well as the sampling behavior. This set of features is very much still under development.  
     */
    struct modelConfiguration
    {
        int reinfectionMode;
        int samplingMode; 
    };



    /**
     * The model context class provides the central interface of libspatialSEIR. After instantiation, the populate
     * method must be called before any meaningful work can be done. Once populated, ModelContext holds pointers to 
     * all of the required FullConditional instances as well as the data they use. In addition, ModelContext provides
     * calculation functions and holds pointers to additional utility classes such as OCLProvider, IOProvider, and 
     * RandomNumberProvider.
     */
    class ModelContext
    {
        public:
            //Methods
            ModelContext();
            ~ModelContext(); 
            /**
             * The populate method is the entry point to working with spatial SEIR models in libspatialSEIR. 
             */
            void populate(InitData* _A0, /**< A0 must be an instance of InitData, which 
                                              contains the starting values of S0, E0, I0, R0*/
                          covariateArgs* xArgs, /**< xArgs is an instance of covariateArgs which contains
                                                     the dimensions of the fixed and time varying covariate 
                                                     matrices driving the exposure process. 
                                                    */
                          covariateArgs* xPrsArgs, /**< xPrsArgs is an instance of covariateArgs which contains
                                                     the dimensions of the fixed and time varying covariate 
                                                     matrices driving the reinfection process. 
                                                    */
                          compartmentArgs* S_starArgs, /**< S_starArgs is an instance of compartmentArgs
                                                            containing the dimensions and data for S_star,
                                                            along with the steadyStateConstraintPrecision
                                                            parameter.*/ 
                          compartmentArgs* E_starArgs, /**< E_starArgs is an instance of compartmentArgs
                                                            containing the dimensions and data for E_star,
                                                            along with the steadyStateConstraintPrecision
                                                            parameter.*/ 
                          compartmentArgs* I_starArgs, /**< I_starArgs is an instance of compartmentArgs
                                                            containing the dimensions and data for I_star,
                                                            along with the steadyStateConstraintPrecision
                                                            parameter.*/ 
                          compartmentArgs* R_starArgs,/**< R_starArgs is an instance of compartmentArgs
                                                            containing the dimensions and data for R_star,
                                                            along with the steadyStateConstraintPrecision
                                                            parameter.*/ 
                          distanceArgs* rawDistArgs, /**< rawDistAgs is the distanceArgs struct 
                                                          containing the data and dimension of the unscaled
                                                          distance matrix. 
                                                        */
                          scaledDistanceArgs* scaledDistArgs, /**< scaledDistArgs is the scaledDistanceArgs struct 
                                                                   containing the data and dimension of the unscaled
                                                                   distance matrix. 
                                                               */
                          double* rho, /**< rho gives the starting value of the spatial dependence parameter*/
                          double* beta,/**< beta gives the starting vector of regression parameters driving the exposure process.*/ 
                          double* p_ei,/**< p_ei is the starting value of the exposed to infectious transition probability. */
                          double* p_ir,/**< p_ir is the starting value of the infectious to removed/recovered transition probability.*/ 
                          double* betaPrs,/**< betaPrs is the starting vector of regression parameters driving the reinfection process.*/ 
                          int* N,/**< N is the matrix of population sizes, corresponding in dimension to the CompartmentalModelMatrix instances.*/ 
                          sliceParameters* sliceWidths, /**< sliceWidths is an instance of sliceParameters which gives the slice sampling widths/
                                                            Metropolis-Hastings tuning parameters for the various FullConditional distributions.*/
                          priorControl* priorInformation, /**< The priorControl struct gives the prior precisions for the regression parameters, 
                                                             and the prior alpha and beta parameters for the p_ei and p_ir distributions.*/
                          modelConfiguration _config /**< The modelConfiguration struct gives the reinfection mode (SEIRS, SEIR, SSEIR) and 
                                                          the MCMC sampling mode. MCMC sampling modes are under very active development. */
                              );

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
            void setSamplingMode(int mode);
            int getSamplingMode();



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

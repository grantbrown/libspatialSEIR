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

            /** The zero parameter overload of calculateS_CPU calculates the S compartment from A0, S_star, and E_star.*/
            void calculateS_CPU(); 

            /** The two parameter overload of calculateS_CPU again calculates the S compartment 
                from A0, S_star and E_star, but assumes that only values of S_star and E_star 
                at column (location) startLoc and after row (time) startTime have changed since 
                S was last calculated. When applicable, this saves considerable time. */
            void calculateS_CPU(int startLoc, int startTime); 

            /** The zero parameter overload of calculate S_givenE_CPU takes advantage of the fact that N=S+E+I+R to 
                update S when E changes.*/
            void calculateS_givenE_CPU(); 

            /** The two parameter overload of calculate S_givenE_CPU takes 
                advantage of the fact that N=S+E+I+R to update S when E changes,
                while assuming that E has only changed for location startLoc and after time
                point startTime.*/
            void calculateS_givenE_CPU(int startLoc, int startTime); 

            /** calculateS_OCL works identically to calculateS_CPU using calls to oclProvider to perform 
                computations in parallel.*/
            void calculateS_OCL(); 

            /** The zero parameter overload of calculateE_CPU calculates the E compartment from A0, I_star, and E_star.*/
            void calculateE_CPU(); 

            /** The two parameter overload of calculateE_CPU again calculates the E compartment 
                from A0, I_star and E_star, but assumes that only values of I_star and E_star 
                at column (location) startLoc and after row (time) startTime have changed since 
                E was last calculated. When applicable, this saves considerable time. */
            void calculateE_CPU(int startLoc, int startTime);

            /** The zero parameter overload of calculate E_givenI_CPU takes advantage of the fact that N=S+E+I+R to 
                update E when I changes.*/
            void calculateE_givenI_CPU(); 

            /** The two parameter overload of calculate E_givenI_CPU takes 
                advantage of the fact that N=S+E+I+R to update E when I changes,
                while assuming that I has only changed for location startLoc and after time
                point startTime.*/
            void calculateE_givenI_CPU(int startLoc, int startTime);

            /** calculateE_OCL works identically to calculateS_CPU using calls to oclProvider to perform 
                computations in parallel.*/
            void calculateE_OCL();  

            /** The zero parameter overload of calculateI_CPU calculates the I compartment from A0, I_star, and R_star.*/
            void calculateI_CPU();

            /** The two parameter overload of calculateI_CPU again calculates the I compartment 
                from A0, I_star and R_star, but assumes that only values of I_star and R_star 
                at column (location) startLoc and after row (time) startTime have changed since 
                I was last calculated. When applicable, this saves considerable time. */
            void calculateI_CPU(int startLoc, int startTime);

            /** The zero parameter overload of calculate I_givenR_CPU takes advantage of the fact that N=S+E+I+R to 
                update I when R changes.*/
            void calculateI_givenR_CPU(); 

            /** The two parameter overload of calculate I_givenR_CPU takes 
                advantage of the fact that N=S+E+I+R to update I when R changes,
                while assuming that R has only changed for location startLoc and after time
                point startTime.*/
            void calculateI_givenR_CPU(int startLoc, int startTime);

            /** calculateI_OCL works identically to calculateI_CPU using calls to oclProvider to perform 
                computations in parallel.*/
            void calculateI_OCL(); 

            /** The zero parameter overload of calculateR_CPU calculates the R compartment from A0, R_star, and S_star.*/
            void calculateR_CPU();

            /** The two parameter overload of calculateR_CPU again calculates the R compartment 
                from A0, R_star and S_star, but assumes that only values of R_star and S_star 
                at column (location) startLoc and after row (time) startTime have changed since 
                I was last calculated. When applicable, this saves considerable time. */
            void calculateR_CPU(int startLoc, int startTime);

            /** The zero parameter overload of calculate R_givenS_CPU takes advantage of the fact that N=S+E+I+R to 
                update R when S changes.*/
            void calculateR_givenS_CPU();

            /** The two parameter overload of calculate R_givenS_CPU takes 
                advantage of the fact that N=S+E+I+R to update R when S changes,
                while assuming that S has only changed for location startLoc and after time
                point startTime.*/
            void calculateR_givenS_CPU(int startLoc, int startTime);

            /** calculateR_OCL works identically to calculateR_CPU using calls to oclProvider to perform 
                computations in parallel.*/
            void calculateR_OCL();

            /**< The four parameter overload of calculateGenericCompartment_CPU takes advantege of the 
                 common calculation pattern for S,E,I, and R to provide a single function mapping starting and 
                 transition values to the filan CompartmentalModelMatrix values. */
            void calculateGenericCompartment_CPU(CompartmentalModelMatrix *comp, /**< pointer to CompartmentalModelMatrix to be calculated.*/ 
                                                 int *comp0, /**< pointer to integer array containing starting values for comp.*/
                                                 CompartmentalModelMatrix *compStarAdd, /**< pointer to transition matrix into comp.*/
                                                 CompartmentalModelMatrix *compStarSub /**< pointer to transition CompartmentalModelMatrix out of comp.*/
                                                 );  

            /** The six parameter overload of calculateGenericCompartment_CPU takes advantege of the 
                common calculation pattern for S,E,I, and R to provide a single function mapping starting and 
                transition values to the filan CompartmentalModelMatrix values, assuming that only location startLoc and time points after startTime
                need be considered. */
            void calculateGenericCompartment_CPU(CompartmentalModelMatrix *comp, /**< pointer to CompartmentalModelMatrix to be calculated*/ 
                                                 int *comp0,/**< pointer to integer array containing starting values for comp*/
                                                 CompartmentalModelMatrix *compStarAdd, /**< pointer to transition matrix into comp*/ 
                                                 CompartmentalModelMatrix *compStarSub,/**< pointer to transition CompartmentalModelMatrix out of comp*/
                                                 int startLoc, /**< location to update*/
                                                 int startTime /**< time after which update is needed*/
                                                 );

            /** calculateGenericCompartment_OCL works identically to calculateGenericCompartment_CPU while using
               calls to oclProvider to perform calculations in parallel.*/
            void calculateGenericCompartment_OCL(int *comp, /**< pointer to CompartmentalModelMatrix to be calculated*/ 
                                                 int *comp0,/**< pointer to integer array containing starting values for comp*/
                                                 int *compStarAdd, /**< pointer to transition matrix into comp*/ 
                                                 int *compStarSub /**< pointer to transition CompartmentalModelMatrix out of comp*/
                                                 );

            /** checkCompartmentBounds is a debug mode function which checks for impossible compartment 
                values and prints errors accordingly. */
            int checkCompartmentBounds(); 

            /** printFCValues is a semi-depricated debug mode function which displays the current full conditional 
                likelihood values of the model parameters.*/
            void printFCValues(); 

            /** setRandomSeed sets the seed used by the pseudorandom number generator provided by
                RandomNumberProvider*/
            void setRandomSeed(unsigned int seedValue); 
              
            /** simulationIter runs a single MCMC update given the current model state. Optionally,
                verbose and debug level output is available.*/
            void simulationIter(bool verbose, bool debug); 
                                                  
            /** runSimulation runs nIterations MCMC updates given the current 
                model state. Optionally, verbose and debug level output is available.*/
            void runSimulation(int nIterations,  bool verbose, bool debug);             

            /** updateSamplingParameters moves the various slice sampling widths / Metropolis-Hastings tuning 
                parameters up or down depending on the current acceptance rate, and re-sets the acceptance 
                counters. 
             */ 
            void updateSamplingParameters(double desiredRatio, /**< desired MCMC acceptance ratio*/ 
                                          double targetWidth, /**< target ratio tolerance (how close is close enough)*/ 
                                          double proportionChange /**< proportion to change the sampling parameters by*/
                                          );

            /** setSamplingMode sets... the sampling mode. This part of the API is in flux. 
             */
            void setSamplingMode(int mode);

            /** getSamplingMode returns the current sampling mode as an integer */
            int getSamplingMode();





            /** cacheP_SE_Calculation places the pre-requisites for 
             * calculating p_se in p_se_components. This is useful 
             * when partial updates can be performed by calculateP_SE_CPU(int startLoc, int startTime)
             */
            void cacheP_SE_Calculation();

            /** The zero parameter overload of calculateP_SE_CPU calculates the exposure probability p_se. 
             */
            void calculateP_SE_CPU();

            /** The two parameter overload of calculateP_SE_CPU calculates the exposure probability starting 
             * at startLoc and startTime. cacheP_SE_Calculation must be called first.
             */
            void calculateP_SE_CPU(int startLoc, int startTime);

            /** calculateP_SE_OCL works identically to calculateP_SE_CPU, while making calls to oclProvider
                to perform calculations in parallel.*/
            void calculateP_SE_OCL();

            /** calculateP_RS_CPU calculates the reinfection probability p_rs.*/
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

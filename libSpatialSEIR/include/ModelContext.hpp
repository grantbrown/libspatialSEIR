#ifndef SPATIALSEIR_INCLUDEFILES
#define SPATIALSEIR_INCLUDEFILES

#include<cstring>
#include<vector>
#endif

#ifndef MODEL_CONTEXT_INC
#define MODEL_CONTEXT_INC

#define LSS_DEGENERATE_DATA_MODEL 1
#define LSS_OVERDISPERSED_DATA_MODEL 2

#include<Eigen/Core>
namespace SpatialSEIR
{

    typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> MatrixType;
    typedef Eigen::Map<MatrixType, Eigen::ColMajor> MatrixMapType;
    typedef Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> IntMatrixType;
    typedef Eigen::Map<IntMatrixType, Eigen::ColMajor> IntMatrixMapType;

    class FullConditional;
    class IterationTask;
    class SetCompartmentSamplingIndicesTask;
    class PerformDecorrelationStep;
    class FC_S0;
    class FC_E0;
    class FC_I0;
    class FC_R0;
    class FC_S_Star;
    class FC_E_Star;
    class FC_R_Star;
    class FC_Beta;
    class FC_Gamma_EI;
    class FC_Gamma_IR;
    class FC_Beta_P_RS;
    class FC_Rho;
    class FC_I_Star_overdispersed;
    class FC_Phi;
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
    struct priorControl;

    /**
     * The modelConfiguration struct stores information on the reinfection mode requested for the model
     * as well as the sampling behavior and data model. This set of features is very much still under development.  
     */
    struct modelConfiguration
    {
        int reinfectionMode;
        int compartmentSamplingMode; 
        int parameterSamplingMode; 
        int indexLength;
        int useDecorrelation;
        int dataModel;
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
                          double* offset, /**< Offset allows irregular spacing between time points.*/

                          int* Y, /**< Y is the actual observed data vector. It must be of the same dimension
                                      as the rest of the Compartments. In the default data model, Y is 
                                      exactly equal to I_star.*/
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
                          scaledDistanceArgs* scaledDistArgs, /**< scaledDistArgs is the scaledDistanceArgs struct 
                                                                   containing the data and dimension of the unscaled
                                                                   distance matrices. 
                                                               */
                          double* rho, /**< rho gives the starting value of the spatial dependence parameter*/
                          double* phi, /**< starting value for overdispersion precision*/
                          double* beta,/**< beta gives the starting vector of regression parameters driving the exposure process.*/ 
                          double* gamma_ei,/**< p_ei is the starting value of the parameter driving the exposed to infectious transition probability. */
                          double* gamma_ir,/**< p_ir is the starting value of the parameter driving the infectious to removed/recovered transition probability.*/ 
                          double* betaPrs,/**< betaPrs is the starting vector of regression parameters driving the reinfection process.*/ 
                          int* N,/**< N is the matrix of population sizes, corresponding in dimension to the CompartmentalModelMatrix instances.*/ 
                          sliceParameters* sliceWidths, /**< sliceWidths is an instance of sliceParameters which gives the slice sampling widths/
                                                            Metropolis-Hastings tuning parameters for the various FullConditional distributions.*/
                          priorControl* priorInformation, /**< The priorControl struct gives the prior precisions for the regression parameters, 
                                                             and the prior alpha and beta parameters for the p_ei and p_ir distributions.*/
                          modelConfiguration _config /**< The modelConfiguration struct gives the reinfection mode (SEIRS, SEIR, SSEIR) and 
                                                          the MCMC sampling mode. MCMC sampling modes are under very active development. */
                              );

            /*buildModel is called at the end of the populate function to fill the ModelContext.model vector with the required 
              full conditional distributions.*/
            void buildModel();

            /** configureIterationTasks determines whether there are any conditions which require additional logical tasks
             * to be performed during MCMC sampling.*/
            void configureIterationTasks();
            
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

            /** setCompartmentSamplingMode sets the sampling mode for the disease compartments. This part of the API is in flux. 
             */
            void setCompartmentSamplingMode(int mode);

            /** getCompartmentSamplingMode returns the current sampling mode as an integer */
            int getCompartmentSamplingMode();

            /** setParameterSamplingMode sets the sampling mode for non compartment parameters. This part of the API is in flux. 
             */
            void setParameterSamplingMode(int mode);

            /** getParameterSamplingMode returns the current sampling mode for non-compartment parameters as an integer */
            int getParameterSamplingMode();





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

            /** calculateP_EI_CPU calculates the E to I transition probability vector p_ei */
            void calculateP_EI_CPU();

            /** calculateP_IR_CPU calculates the I to R transition probability vector p_ir */
            void calculateP_IR_CPU();

            /** calculateP_RS_CPU calculates the reinfection probability p_rs.*/
            void calculateP_RS_CPU();

            /** Calculates the sum of all components of S.*/
            int totalS();
            /** Calculates the sum of all components of E.*/
            int totalE();
            /** Calculates the sum of all components of I.*/
            int totalI();
            /** Calculates the sum of all components of R.*/
            int totalR();

            /** Calculates the sum of all components of S.*/
            int totalS(int tpt);
            /** Calculates the sum of all components of E.*/
            int totalE(int tpt);
            /** Calculates the sum of all components of I.*/
            int totalI(int tpt);
            /** Calculates the sum of all components of R.*/
            int totalR(int tpt);

            /** Calculates the sum of all components of S_star.*/
            int totalS_star();
            /** Calculates the sum of all components of E_star.*/
            int totalE_star();
            /** Calculates the sum of all components of I_star.*/
            int totalI_star();
            /** Calculates the sum of all components of R_star.*/
            int totalR_star();

            /** Calculates the sum of all components of S_star.*/
            int totalS_star(int tpt);
            /** Calculates the sum of all components of E_star.*/
            int totalE_star(int tpt);
            /** Calculates the sum of all components of I_star.*/
            int totalI_star(int tpt);
            /** Calculates the sum of all components of R_star.*/
            int totalR_star(int tpt);

            /** Calculates the average of all components of p_se.*/
            double avgP_SE();
            /** Calculates the average of all components of p_rs.*/
            double avgP_RS();
            /** Calculates the average of all components of p_se for a particular time point..*/
            double avgP_SE(int tpt);

            /** Estimates the average basic reproductive number across all time periods and spatial locations.*/
            double estimateR0();
            /** Estimates the basic reproductive number for time point t*/
            double* estimateR0(int t);
            /** Estimates the effective reproductive rate across all time periods and spatial locations*/
            double estimateEffectiveR0();
            /** Estimates the effective reproductive rate at time t.*/
            double* estimateEffectiveR0(int t);
            /** Calculates the location specific components of the effective R0 quantity at time t*/
            double* calculateR0Components(int t);
            /** Calculates the location specific components of the effective R0 quantity at time t*/
            double* calculateEffectiveR0Components(int t);
            /** Calculates the next generation matrix */
            double* calculateG(int t);
            /** Calculates the next generation matrix, integrated forward in time. */
            double* calculateIntegratedG(int t);
    



            //Logic provider and utility classes
            /** pointer to an IOProvider instance which writes selected sampler output to a comma delimited text file.*/
            IOProvider *fileProvider;
            /** pointer to a RandomNumberProvider instance which provides random numbers from various distributions as well
                implementations of several probability density functions.*/
            RandomNumberProvider *random;
            /** pointer to an instance of OCLProvider which supplies an interface for selecting active OpenCL devices as well
                as parallelized computations required for various FullConditional distributions.*/
            OCLProvider *oclProvider; 
            /** Pointer to FullConditional distribution for the initial suscptible population.*/
            FC_S0 *S0_fc;
            /** Pointer to FullConditional distribution for the initial exposed population - not used*/
            FC_E0 *E0_fc;
            /** Pointer to FullConditional distribution for the initial infectious population*/
            FC_I0 *I0_fc;
            /** Pointer to FullConditional distribution for the initial recovered population - not used*/
            FC_R0 *R0_fc;
            /** Pointer to FullConditional distribution for the transition matrix S_star*/
            FC_S_Star *S_star_fc;
            /** Pointer to FullConditional distribution for the transition matrix E_star*/
            FC_E_Star *E_star_fc;
            /** Pointer to FullConditional distribution for the transition matrix R_star*/
            FC_R_Star *R_star_fc;
            /** Pointer to FullConditional distribution for the regression parameters beta*/
            FC_Beta *beta_fc;
            /** Pointer to FullConditional distribution for the spatial depenence parameter rho*/
            FC_Rho *rho_fc;
            /** Pointer to FullConditional distribution for the regression parameters betaP_RS*/
            FC_Beta_P_RS *betaPrs_fc;
            /** Pointer to FullConditional distribution for the parameter transition probability p_ei*/
            FC_Gamma_EI *gamma_ei_fc;
            /** Pointer to FullConditional distribution for the parameter controlling transition probability p_ir*/
            FC_Gamma_IR *gamma_ir_fc;
            /** Pointer to FullConditional distribution for the I_star transition matrix under overdispersion*/
            FC_I_Star_overdispersed *I_star_overdispersed_fc;
            /** Pointer to FullConditional distribution for the overdispersion parameter for I_star*/
            FC_Phi *phi_fc;

            /** Pointer to the sampling indices task.*/
            SetCompartmentSamplingIndicesTask* setSamplingIndicesTask;
            /** Pointer to the decorrelation step task*/
            PerformDecorrelationStep* decorrelationStepTask;

            //Data
            /** Pointer to the actual data vector: Y. Y is related by a data model to I_star.*/
            int* Y;
            /** Pointer to the CompartmentalModelMatrix data structure storing the S compartment.*/ 
            CompartmentalModelMatrix* S;
            /** Pointer to the CompartmentalModelMatrix data structure storing the E compartment.*/ 
            CompartmentalModelMatrix* E;
            /** Pointer to the CompartmentalModelMatrix data structure storing the I compartment.*/ 
            CompartmentalModelMatrix* I;
            /** Pointer to the CompartmentalModelMatrix data structure storing the R compartment.*/ 
            CompartmentalModelMatrix* R;
            /** Pointer to the CompartmentalModelMatrix data structure storing the S_star compartment.*/ 
            CompartmentalModelMatrix* S_star;
            /** Pointer to the CompartmentalModelMatrix data structure storing the E_star compartment.*/ 
            CompartmentalModelMatrix* E_star;
            /** Pointer to the CompartmentalModelMatrix data structure storing the I_star compartment.*/ 
            CompartmentalModelMatrix* I_star;
            /** Pointer to the CompartmentalModelMatrix data structure storing the R_star compartment.*/ 
            CompartmentalModelMatrix* R_star;
            /** Pointer to InitData instance containing the data for the initial compartment sizes.*/ 
            InitData* A0;
            /** Pointer to CovariateMatrix instance containing covariate informating driving the exposure process.*/
            CovariateMatrix* X;
            /** Pointer to CovariateMatrix instance containing covariate informating driving the reinfection process.*/
            CovariateMatrix* X_pRS;
            /** Pointer to vector of scaled DistanceMatrix  objects. */
            std::vector<DistanceMatrix*>* scaledDistMatrices;
            /** Extra compartment storage for caching integer computations*/
            CompartmentalModelMatrix* tmpContainer;
            /** Pointer to modelConfiguration instance.*/
            modelConfiguration* config;
            /** A model in libspatialSEIR is a collection of full conditional distributions which must be sampled from. This collection 
             * is assembled during ModelContext.populate based on the reinfection mode, data model (not yet implemented), and the detection
             * of special cases such as the single location case. */
            std::vector<FullConditional*>* model; 
            /** Iteration tasks are bits of logical code that need to execute periodically during MCMC sampling, but don't fall into the usual
             * FullConditional framework.*/
            std::vector<IterationTask*>* iterationTasks;
            /** Storage for the regression parameters beta*/
            double* beta;
            /** Storage for the regression parameters betaP_RS*/
            double* betaPrs;
            /** Storage for the vector of spatial dependence parameters, rho*/
            double* rho;
            /** Storage for the overdispersion parameter*/
            double* phi;
            /** Storage for the external infection parameters gamma (DEPRICATED)*/
            double* gamma;
            /** Storage for the linear predictor, eta*/
            double* eta;
            /** Storage for the vector of offsets.*/
            double* offset;
            /** Storage for the exposure probability matrix, p_se*/
            double* p_se;
            /** Storage for cacheable portions of the p_se calculation*/
            double* p_se_components;
            /** ComparmentalModelMatrix sized cache of doubles. */
            double* compartmentCache;
            /** Length of the indexList array.*/
            int* indexLength;
            /** List of indices, used by some sampling functions. */
            int* indexList;
            /** Parameter controlling transition probability vector p_ei*/
            double* gamma_ei;
            /** Transition probability vector p_ei*/
            double* p_ei;
            /** Parameter controlling transition probability vector p_ir*/
            double* gamma_ir;
            /** Transition probability vector p_ir*/
            double* p_ir;
            /** Transition probability vector p_rs*/
            double* p_rs;
            /** Matrix of population sizes */
            int* N;
            /** Indicator that ModelContext instance is ready for use.*/
            int* isPopulated;
            /** Indicator that ModelContext is operating in the single location special case.*/
            int* singleLocation;
            /** Number of MCMC samples so far.*/
            int* numIterations;
            /** Flag to indicate whether or not decorrelation should be used. */
    };
}
#endif

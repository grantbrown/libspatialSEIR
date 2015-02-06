#include <ModelContext.hpp>
#include <LSS_FullConditionalList.hpp>
#include <LSS_Samplers.hpp>
#include <LSS_IterationTasks.hpp>
#include <OCLProvider.hpp>
#include <CompartmentalModelMatrix.hpp>
#include <CovariateMatrix.hpp>
#include <DistanceMatrix.hpp>
#include <RandomNumberProvider.hpp>
#include <IOProvider.hpp>
#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#ifdef LSS_USE_BLAS
	#include <cblas.h>
#endif

#include<cmath>
#include<ctime>
				 
namespace SpatialSEIR
{
	
			
    ModelContext::ModelContext()
    {
        //random = new RandomNumberProvider(static_cast<unsigned int>(std::time(0)));
        isPopulated = new int; *isPopulated = 0;
        numIterations = new int; *numIterations = 0;
    }

    void ModelContext::setRandomSeed(unsigned int seedValue)
    {
        lssCout << "Setting seed.\n";
        delete random;
        random = new RandomNumberProvider(seedValue);
    }

    void ModelContext::setCompartmentSamplingMode(int mode)
    {
        int oldMode = (config -> compartmentSamplingMode); 
        (config -> compartmentSamplingMode) = mode;
        configureIterationTasks();
        unsigned int i;
        try
        {
            for (i = 0; i < model -> size(); i++)
            {
                if ((*model)[i] -> getFullConditionalType() == LSS_PARAMETER_FULL_CONDITIONAL_TYPE)
                {
                    // Do nothing. 
                }
                else if ((*model)[i] -> getFullConditionalType() == LSS_INIT_COMPARTMENT_FULL_CONDITIONAL_TYPE)

                {
                    switch (mode)
                    {
                        case COMPARTMENT_METROPOLIS_SAMPLER:
                            (*model)[i] -> setSamplerType(INITCOMPARTMENT_METROPOLIS_SAMPLER);
                            break;
                        case COMPARTMENT_IDX_METROPOLIS_SAMPLER:
                            (*model)[i] -> setSamplerType(INITCOMPARTMENT_METROPOLIS_SAMPLER);
                            break;
                        case COMPARTMENT_IDX_SLICE_SAMPLER:
                            (*model)[i] -> setSamplerType(INITCOMPARTMENT_METROPOLIS_SAMPLER);
                            break;
                        case COMPARTMENT_METROPOLIS_SAMPLER_OCL:
                            (*model)[i] -> setSamplerType(INITCOMPARTMENT_METROPOLIS_SAMPLER_OCL);
                            break;
                        default:
                            (*model)[i] -> setSamplerType(INITCOMPARTMENT_METROPOLIS_SAMPLER);
                            break;
                    }
                }
                else if ((*model)[i] -> getFullConditionalType() == LSS_COMPARTMENT_FULL_CONDITIONAL_TYPE)
     
                {
                    (*model)[i] -> setSamplerType(mode);
                }
            }
        }
        catch (int e){
            lssCout << "Failed to update compartment mode. Re-setting to: " << oldMode << "\n";
            (config -> parameterSamplingMode) = oldMode;
            configureIterationTasks();
            setCompartmentSamplingMode(oldMode);
        }            
    }

    int ModelContext::getCompartmentSamplingMode()
    {
        return((config -> compartmentSamplingMode));
    }

    void ModelContext::setParameterSamplingMode(int mode)
    {
        int oldMode = (config -> parameterSamplingMode); 
        (config -> parameterSamplingMode) = mode;
        configureIterationTasks();
        try{
            unsigned int i;
            for (i = 0; i < model -> size(); i++)
            {
                if ((*model)[i] -> getFullConditionalType() == LSS_PARAMETER_FULL_CONDITIONAL_TYPE)
                {
                    (*model)[i] -> setSamplerType(mode);
                } 
            }
        }
        catch(int e){
            lssCout << "Failed to set parameter sampling mode, resetting to: " << oldMode << "\n";
            (config -> parameterSamplingMode) = oldMode;
            configureIterationTasks();
            setParameterSamplingMode(oldMode); 
        }
    }

    int ModelContext::getParameterSamplingMode()
    {
        return((config -> parameterSamplingMode));
    }

    void ModelContext::populate(InitData* _A0,
                                covariateArgs* xArgs, 
                                covariateArgs* xPrsArgs,
                                double *offset_,
                                int* Y_,
                                compartmentArgs* S_starArgs,
                                compartmentArgs* E_starArgs,
                                compartmentArgs* I_starArgs,
                                compartmentArgs* R_starArgs,
                                scaledDistanceArgs* scaledDistArgs,
                                double* rho_, double* phi_, double* beta_, 
                                double* gamma_ei_, double* gamma_ir_, double* betaPrs_, 
                                int* N_, sliceParameters* sliceWidths,
                                priorControl* priorValues,
                                modelConfiguration config_)
    {
        // Allocate Objects
        random = new RandomNumberProvider(12312415);
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
        X_pRS = new CovariateMatrix();
        scaledDistMatrices = new std::vector<DistanceMatrix*>();
        rho = new double[(scaledDistArgs -> inData).size()];
        phi = new double; *phi = *phi_;
        fileProvider = new IOProvider();
        singleLocation = new int; *singleLocation = -1;
        oclProvider = new OCLProvider();
        model = new std::vector<FullConditional*>;
        iterationTasks = new std::vector<IterationTask*>;
        config = new modelConfiguration();
        *config = config_;
        p_se = new double[*(S_starArgs -> inRow)*(*(S_starArgs -> inCol))];
        p_se_components = new double[*(S_starArgs -> inRow)*(*(S_starArgs -> inCol))];
        compartmentCache = new double[*(S_starArgs -> inRow)*(*(S_starArgs -> inCol))];
        Y = new int[*(S_starArgs -> inRow)*(*(S_starArgs -> inCol))];

        // Create empty index cache
        indexLength = new int; *indexLength = config -> indexLength;
        indexList = new int[*(S_starArgs -> inRow)*(*(S_starArgs -> inCol))];
        memset(indexList, 0, (*(S_starArgs -> inRow)*(*(S_starArgs -> inCol)))*sizeof(int));

        p_ei = new double[*(S_starArgs -> inRow)];
        p_ir = new double[*(S_starArgs -> inRow)];
        p_rs = new double[*(S_starArgs -> inRow)];
        gamma_ei = new double;
        gamma_ir = new double;
        offset = new double[*(S_starArgs -> inRow)];
        N = new int[(*(S_starArgs -> inRow))*(*(S_starArgs -> inCol))];                                          
        int nbeta = (xArgs -> inCol_x);
        int neta = (xArgs -> inRow_x);
        int nBetaPrs = (xPrsArgs -> inCol_x);



        beta = new double[nbeta];
        betaPrs = new double[nBetaPrs];
        eta = new double[neta];
        //Depricated
        gamma = new double[*(S_starArgs -> inRow)];
        memset(gamma, 0, *(S_starArgs -> inRow)*sizeof(double));
        // Create empty compartment for calculation.
        tmpContainer = new CompartmentalModelMatrix();
        tmpContainer -> createEmptyCompartment((S_starArgs -> inRow), (S_starArgs -> inCol));

        *singleLocation = ((*(S_starArgs -> inCol)) > 1 ? 0 : 1);

        // Initialize Stuff
        A0 -> populate(_A0 -> S0,_A0 -> E0,_A0 -> I0,_A0 -> R0,_A0 -> numLocations);
        X -> genFromDataStream(xArgs -> inData_x, 
                               xArgs -> inRow_x,
                               xArgs -> inCol_x); 

        X_pRS -> genFromDataStream(xPrsArgs -> inData_x, 
                                   xPrsArgs -> inRow_x,
                                   xPrsArgs -> inCol_x); 

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

        int i;
        for (i = 0; i < (int) (scaledDistArgs -> inData).size(); i++)
        {
            rho[i] = rho_[i];
        }

        for (i = 0; i < (int) (scaledDistArgs -> inData).size(); i++)
        {
            DistanceMatrix* scaledTmp = new DistanceMatrix();
            scaledTmp -> genFromDataStream((scaledDistArgs -> inData)[i], scaledDistArgs -> dim);
            scaledDistMatrices -> push_back(scaledTmp);
        }

        // Initialize Data
        for (i = 0; i < *(S_star -> nrow); i++)
        {
            offset[i] = offset_[i]; 
        }

        for (i = 0; i < ((*(S_star -> nrow)) * (*(S_star -> ncol))); i++)
        {
            Y[i] = Y_[i]; 
        }

        for (i = 0; i < nbeta; i++)
        {
            beta[i] = beta_[i];
        }

        for (i = 0; i < neta; i++)
        {
            eta[i] = 0.0;
        }

        for (i = 0; i < nBetaPrs; i++)
        {
            betaPrs[i] = betaPrs_[i];
        }
        for (i = 0; i< (*(S -> nrow))*(*(S->ncol)); i++)
        {
            N[i] = N_[i];
        } 

        *gamma_ei = *gamma_ei_;
        *gamma_ir = *gamma_ir_;


        // Wire up the full conditional classes 
        S0_fc = new FC_S0(this,
                          S,
                          E,
                          E_star,
                          I_star,
                          A0,
                          p_se,
                          p_ei,
                          *(sliceWidths -> S0Width));

        E0_fc = new FC_E0(this,
                          S,
                          E,
                          I,
                          E_star,
                          I_star,
                          R_star,
                          A0,
                          p_ir,
                          p_ei,
                          p_se,
                          *(sliceWidths -> S0Width));

        I0_fc = new FC_I0(this, 
                          S,
                          I,
                          R,
                          S_star,
                          E_star,
                          R_star,
                          A0,
                          p_ir,
                          p_rs,
                          p_se,
                          *(sliceWidths -> I0Width));

        R0_fc = new FC_R0(this,
                          R,
                          S,
                          S_star,
                          E_star,
                          R_star,
                          A0,
                          p_rs,
                          p_se, 
                          *(sliceWidths -> I0Width));

        S_star_fc = new FC_S_Star(this,
                                  S_star,
                                  S,
                                  R,
                                  E_star,
                                  R_star,
                                  A0,
                                  X,
                                  p_se,
                                  p_rs,
                                  beta,
                                  rho,
                                  (S_starArgs -> steadyStateConstraintPrecision),
                                  *(sliceWidths -> S_starWidth));

        E_star_fc = new FC_E_Star(this,
                                  E_star,
                                  E,
                                  S,
                                  I_star,
                                  X,A0,p_se,p_ei,
                                  rho,beta,
                                  (E_starArgs -> steadyStateConstraintPrecision),
                                  *(sliceWidths -> E_starWidth));

        I_star_fc = new FC_I_Star(this,
                                  E_star,
                                  I_star,
                                  R_star,
                                  E,
                                  I,
                                  A0,p_ei,p_ir,
                                  (E_starArgs -> steadyStateConstraintPrecision),
                                  *(sliceWidths -> E_starWidth));

        R_star_fc = new FC_R_Star(this,
                                  R_star,
                                  R,
                                  I,
                                  S_star,
                                  E_star,
                                  I_star,
                                  S,
                                  A0,p_rs,p_ir,p_se,
                                  (R_starArgs -> steadyStateConstraintPrecision),
                                  *(sliceWidths -> R_starWidth));

        beta_fc = new FC_Beta(this,
                              E_star,
                              S,
                              A0,X,p_se,beta,rho,
                              *(sliceWidths -> betaWidth),
                              (priorValues -> betaPriorPrecision),
                              (priorValues -> betaPriorMean));

        p_se_fc = new FC_P_SE(this,
                              E_star,
                              S,
                              A0,X,p_se,beta,rho,
                              *(sliceWidths -> betaWidth),
                              priorValues -> betaPriorPrecision,
                              priorValues -> betaPriorMean,
                              scaledDistArgs -> priorAlpha_rho,
                              scaledDistArgs -> priorBeta_rho);

        rho_fc = new FC_Rho(this,
                            E_star,
                            S,
                            A0,X,p_se,beta,rho,
                            *(sliceWidths -> rhoWidth),
                            scaledDistArgs -> priorAlpha_rho,
                            scaledDistArgs -> priorBeta_rho);

        betaPrs_fc = new FC_Beta_P_RS(this,S_star,R,X_pRS,A0,p_rs,betaPrs, 
                                      *(sliceWidths -> betaPrsWidth),
                                       (priorValues->betaPrsPriorPrecision), 
                                       (priorValues->betaPrsPriorMean)
                                      );

        gamma_ei_fc = new FC_Gamma_EI(this,
                              I_star,
                              E,
                              A0,p_ei,
                              gamma_ei,
                              (priorValues -> P_EI_priorAlpha),
                              (priorValues -> P_EI_priorBeta),
                              *(sliceWidths -> gammaEiWidth)
                              );

        gamma_ir_fc =  new FC_Gamma_IR(this,
                             R_star,
                             I,
                             A0,p_ir,
                             gamma_ir,
                             (priorValues -> P_IR_priorAlpha),
                             (priorValues -> P_IR_priorBeta),
                             *(sliceWidths -> gammaIrWidth));

        phi_fc = new FC_Phi(this,
                            I_star,
                            phi,
                            (priorValues -> Phi_priorAlpha),
                            (priorValues -> Phi_priorBeta),
                            Y,
                            *(sliceWidths -> phiWidth));

        I_star_overdispersed_fc = new FC_I_Star_overdispersed(this,
                                                              Y,
                                                              I_star,
                                                              I,
                                                              E,
                                                              R_star,
                                                              p_ei,
                                                              p_ir,
                                                              phi,
                                                              (R_starArgs -> steadyStateConstraintPrecision),
                                                              *(sliceWidths -> E_starWidth));

        R_star_overdispersed_fc = new FC_R_Star_overdispersed(this,
                                                Y,
                                                R_star,
                                                R,
                                                I,
                                                S_star,
                                                E_star,
                                                I_star,
                                                S,
                                                A0,p_rs,p_ir,p_se,phi,
                                                (R_starArgs -> steadyStateConstraintPrecision),
                                                *(sliceWidths -> R_starWidth));

        
        // Set up iteration tasks
        setSamplingIndicesTask = new SetCompartmentSamplingIndicesTask(this);

        if (!(*singleLocation))
        {
            performHybridSE_EI_UpdateTask = new PerformHybridSE_EI_UpdateStep(this, gamma_ei_fc, p_se_fc,100);
            decorrelationStepTask = new PerformDecorrelationStep(this, 100, p_se_fc);
        }
        else
        {
            performHybridSE_EI_UpdateTask = new PerformHybridSE_EI_UpdateStep(this, gamma_ei_fc, beta_fc,100);
            decorrelationStepTask = new PerformDecorrelationStep(this, 100, beta_fc);
        }
        performHybridIR_RS_UpdateTask = new PerformHybridIR_RS_UpdateStep(this, gamma_ir_fc, betaPrs_fc,100);


        // Calculate Compartments
        this -> calculateS_CPU();
        this -> calculateE_CPU();
        this -> calculateI_CPU();
        this -> calculateR_CPU();
        this -> calculateP_RS_CPU();
        this -> calculateP_SE_CPU();
        this -> calculateP_EI_CPU();
        this -> calculateP_IR_CPU();

        this -> buildModel();
        *isPopulated = 1;


        this -> setCompartmentSamplingMode(COMPARTMENT_METROPOLIS_SAMPLER);
        this -> setParameterSamplingMode(PARAMETER_JOINT_METROPOLIS_SAMPLER);
    }

    void ModelContext::configureIterationTasks()
    {
        // Clear iterationTasks queue. 
        unsigned int i, sz;
        sz = iterationTasks -> size();
        for (i = 0; i < sz; i++){iterationTasks -> pop_back();}
        if  ((config -> compartmentSamplingMode) == COMPARTMENT_IDX_METROPOLIS_SAMPLER || 
             (config -> compartmentSamplingMode) == COMPARTMENT_IDX_SLICE_SAMPLER ||
             (config -> compartmentSamplingMode) == COMPARTMENT_BINOM_IDX_METROPOLIS_SAMPLER || 
             (config -> compartmentSamplingMode) == COMPARTMENT_BINOM_MIXED_SAMPLER)
        {
            iterationTasks -> push_back(setSamplingIndicesTask);
        } 

        if  ((config -> useDecorrelation) > 0) 
        {
            iterationTasks -> push_back(decorrelationStepTask);
            *(decorrelationStepTask -> iterationCount) = (config -> useDecorrelation);
            
        }
        if (config -> performHybridStep > 0)
        {
            iterationTasks -> push_back(performHybridSE_EI_UpdateTask);
             *(performHybridSE_EI_UpdateTask -> iterationCount) = (config -> performHybridStep); 

            if ((config -> reinfectionMode) == 1)
            {
                iterationTasks -> push_back(performHybridIR_RS_UpdateTask);
                 *(performHybridIR_RS_UpdateTask -> iterationCount) = (config -> performHybridStep); 
            }
        }
    }

    void ModelContext::buildModel()
    {
        // Configure iteration tasks
        configureIterationTasks();

        // clear the previous model
        unsigned int i;
        unsigned int modelSize = model -> size();
        for (i = 0; i < modelSize; i++)
        {
            model -> pop_back();
        }
        if ((model -> size()) != 0)
        {
            lssCout << "Error clearing model.\n";
            throw(-1);
        }

        CompartmentFullConditional* nonDataModel;
        if (config -> dataModelCompartment == 0) 
        {
            nonDataModel = R_star_fc;
        }
        else
        {
            nonDataModel = I_star_fc;
        }

        // build the model here. 
        if ((config -> reinfectionMode) == 1)
        {
            model -> push_back(S0_fc);
            model -> push_back(I0_fc);
            model -> push_back(S_star_fc);
            model -> push_back(E_star_fc);
            model -> push_back(nonDataModel); // I_star or R_star

            model -> push_back(betaPrs_fc);
            if (!(*singleLocation))
            {
                model -> push_back(p_se_fc);
            }
            else
            {
                model -> push_back(beta_fc);
            }
            model -> push_back(gamma_ei_fc);
            model -> push_back(gamma_ir_fc);
        }
        else if (config -> reinfectionMode == 2)
        {
            model -> push_back(S0_fc);
            model -> push_back(I0_fc);
            model -> push_back(S_star_fc);
            model -> push_back(E_star_fc);
            model -> push_back(nonDataModel); // I_star or R_star
            if (!(*singleLocation))
            {
                model -> push_back(p_se_fc);
            }
            else
            {
                model -> push_back(beta_fc);
            }
            model -> push_back(gamma_ei_fc);
            model -> push_back(gamma_ir_fc);
        }
        else if (config -> reinfectionMode  == 3)
        {
            model -> push_back(S0_fc);
            model -> push_back(I0_fc);
            model -> push_back(E_star_fc);
            model -> push_back(nonDataModel); // I_star or R_star
            if (!(*singleLocation))
            {
                model -> push_back(p_se_fc);
            }
            else
            {
                model -> push_back(beta_fc);
            }
            model -> push_back(gamma_ei_fc);
            model -> push_back(gamma_ir_fc);
        } 
        if ((config -> dataModel) == 1)
        {
            if ((config -> dataModelCompartment) == 0)
            {
                model -> push_back(I_star_overdispersed_fc);
            }
            else if ((config -> dataModelCompartment) == 1)
            {
                model -> push_back(R_star_overdispersed_fc);
            }
            else
            {
                lssCout << "Invalid Data Model Compartment Index: " << (config -> dataModelCompartment) <<"\n";
                throw(-1);
            }
            model -> push_back(phi_fc);
        }
    }

    int ModelContext::checkCompartmentBounds()
    {
        int i;
        int err = 0;
        int rowCol = (*(R->ncol))*(*(R->nrow));
        for (i = 0; i < rowCol;i++)
        {
            if ((S_star -> data)[i] > (R -> data)[i])
            {
                lssCout << "S_star too big: " << i << ", val:"<< S_star_fc -> getValue() << " \n";
                err = 1;
                break;
            }
            if ((S_star -> data)[i] < 0)
            {
                lssCout << "S_star <0: " << i << ", val:"<< S_star_fc -> getValue() << " \n";
                err = 1;
                break;
            }
            if ((S -> data)[i] < 0)
            {
                lssCout << "S <0: " << i << " \n";
                err = 1;
                break;
            }

        }
        for (i = 0; i < rowCol;i++)
        {
            if ((E_star -> data)[i] > (S -> data)[i])
            {
                lssCout << "E_star too big: " << i << ", val:"<< E_star_fc -> getValue() << " \n";
                err = 1;
                break;
            }
            if ((E_star -> data)[i] < 0)
            {
                lssCout << "E_star <0: " << i << ", val:"<< E_star_fc -> getValue() << " \n";
                err = 1;
                break;
            }
            if ((E -> data)[i] < 0)
            {
                lssCout << "E <0: " << i << " \n";
                err = 1;
                break;
            }

        }
        for (i = 0; i < rowCol;i++)
        {
            if ((I_star -> data)[i] > (E -> data)[i])
            {
                lssCout << "I_star too big: " << i << "\n";
                err = 1;
                break;
            }
            if ((I_star -> data)[i] < 0)
            {
                lssCout << "I_star <0: " << i << ", val: \n";
                err = 1;
                break;
            }

            if ((I -> data)[i] < 0)
            {
                lssCout << "I_star <0: " << i << " \n";
                err = 1;
                break;
            }

        }
        for (i = 0; i < rowCol;i++)
        {
            if ((R_star -> data)[i] > (I -> data)[i])
            {
                
                lssCout << "R_star too big: " << i << ", val:"<< R_star_fc -> getValue() << " \n";
                err = 1;
                break;
            }
            if ((R_star -> data)[i] < 0)
            {
                lssCout << "R_star <0: " << i << ", val:"<< R_star_fc -> getValue() << " \n";
                err = 1;
                break;
            }
            if ((R -> data)[i] < 0)
            {
                lssCout << "R <0: " << i << " \n";
                err = 1;
                break;
            }
        }
        if (config -> reinfectionMode <= 2)
        {
            int nTpt = *(S_star -> nrow);
            for (i = 0; i < nTpt; i++)
            {
                if (p_rs[i] <= 0 || p_rs[i] >= 1)
                {
                    lssCout << "Invalid P_RS Value at Index" << i << "\n";
                    err = 1;
                }
            }
        }

        return(err);
    }

    void ModelContext::printFCValues()
    {
        beta_fc -> evalCPU();
        if (config -> reinfectionMode == 1)
        {
            betaPrs_fc -> evalCPU();
        }
        gamma_ei_fc -> evalCPU();
        gamma_ir_fc -> evalCPU();
        rho_fc -> evalCPU();
        if (config -> reinfectionMode <= 2)
        {
            S_star_fc -> evalCPU();
        }
        E_star_fc -> evalCPU();
        R_star_fc -> evalCPU();
        lssCout << "  FC Values:\n";
        lssCout << "    Beta: " << beta_fc ->getValue() << "\n";
        if (config -> reinfectionMode == 1)
        {
            lssCout << "    betaPrs: " << betaPrs_fc -> getValue() << "\n";
        }
        lssCout << "    p_ei: " << gamma_ei_fc ->getValue() << "\n";
        lssCout << "    p_ir: " << gamma_ir_fc ->getValue() << "\n";
        if (!(*singleLocation))
        {
            lssCout << "    rho: " << rho_fc ->getValue() << "\n";
        }
        if (config -> reinfectionMode <= 2)
        {
            lssCout << "    S_star: " << S_star_fc -> getValue() << "\n";
        }
        lssCout << "    E_star: " << E_star_fc -> getValue() << "\n";
        lssCout << "    R_star: " << R_star_fc -> getValue() << "\n";
    }

    void ModelContext::updateSamplingParameters(double desiredRatio, double targetWidth, double proportionChange)
    {
        unsigned int i;
        for (i = 0; i < model -> size(); i++)
        {
            (*model)[i] -> updateSamplingParameters(desiredRatio, targetWidth, proportionChange);
        }
    }

    void ModelContext::simulationIter(bool verbose = false, bool debug = false)
    {
        if (debug)
        {
            this -> printFCValues();
            this -> checkCompartmentBounds();

            lssCout << "S: " << S -> marginSum(3,-1) << "\n";
            lssCout << "E: " << E -> marginSum(3,-1) << "\n";
            lssCout << "I: " << I -> marginSum(3,-1) << "\n";
            lssCout << "R: " << R -> marginSum(3,-1) << "\n";

            lssCout << "S_star: " << S_star -> marginSum(3,-1) << "\n";
            lssCout << "E_star: " << E_star -> marginSum(3,-1) << "\n";
            lssCout << "I_star: " << I_star -> marginSum(3,-1) << "\n";
            lssCout << "R_star: " << R_star -> marginSum(3,-1) << "\n";
        }

        unsigned int i;

        for (i = 0; i < model->size(); i++)
        {
           (*model)[i] -> sample(verbose);
        }
        for (i = 0; i < iterationTasks -> size(); i++)
        {
            (*iterationTasks)[i] -> executeTask();
        }
    }

    // Method: runSimulation
    // Accesses: Everything lol
    // Updates: Everything lol
    void ModelContext::runSimulation(int nIterations, bool verbose = false, bool debug = false)
    {
        int i;
        int itrStart = *numIterations;
        int itrMax = nIterations + (*numIterations);
        for (i = itrStart; i < itrMax; i++)
        {
            if (verbose)
            {
                lssCout << "Iteration: " << i << "\n";
            }
            this -> simulationIter(verbose, debug);
            this -> fileProvider -> catIter(i);
            (*numIterations) = (*numIterations + 1);
        }
    }

    // Method: calculateS
    // Accesses: A0, S_star, E_star
    // Updates: S
    void ModelContext::calculateS_CPU()
    {
        calculateGenericCompartment_CPU(S, A0 -> S0,
                                    S_star, E_star);
    }
    void ModelContext::calculateS_CPU(int startLoc, int startTime)
    {
        calculateGenericCompartment_CPU(S, A0 -> S0,
                                    S_star, E_star,
                                    startLoc, startTime);
    }

    void ModelContext::calculateS_givenE_CPU()
    {
        int i;
        int nLoc = *(S -> ncol);
        int nTpt = *(S -> nrow);
        int maxItr = nLoc*nTpt;
        for (i = 0; i < nLoc; i++)
        {
            (A0 -> S0)[i] = N[i*nTpt] - (A0 -> E0)[i] - (A0 -> I0)[i] - (A0 -> R0)[i];
        }

        for (i = 0; i < maxItr; i++)
        {
            (S->data)[i] = N[i] - (E->data)[i] - (I->data)[i] - (R->data)[i]; 
        }
    }

    void ModelContext::calculateS_givenE_CPU(int startLoc, int startTime)
    {
        int i,startIdx,idx;
        int nTpt = *(S -> nrow);
        startIdx = startTime + startLoc*nTpt;
        idx = startIdx;
        if (startTime == 0)
        {
            (A0 -> S0)[startLoc] = (N[startLoc*nTpt] - (A0 -> E0)[startLoc] - (A0 -> I0)[startLoc]
                                    - (A0 -> R0)[startLoc]);
        }

        for (i = startTime; i < *(S->nrow); i++)
        {
            (S -> data)[idx] = N[idx] - (R->data)[idx] - (E->data)[idx] - (I->data)[idx];  
            idx += 1;
        }
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
        calculateGenericCompartment_CPU(E, A0 -> E0,
                                        E_star, I_star);
    }

    void ModelContext::calculateE_CPU(int startLoc, int startTime)
    {
        calculateGenericCompartment_CPU(E, A0 -> E0,
                                        E_star, I_star,
                                        startLoc, startTime);
    }

    void ModelContext::calculateE_givenI_CPU()
    {
        int i;
        int nLoc = *(E -> ncol);
        int nTpt = *(E -> nrow);
        int maxItr = nLoc*nTpt;
        for (i = 0; i < nLoc; i++)
        {
            (A0 -> E0)[i] = N[i*nTpt] - (A0 -> S0)[i] - (A0 -> I0)[i] - (A0 -> R0)[i];
        }

        for (i = 0; i < maxItr; i++)
        {
            (E->data)[i] = N[i] - (S->data)[i] - (I->data)[i] - (R->data)[i]; 
        }
    }


    void ModelContext::calculateE_givenI_CPU(int startLoc, int startTime)
    {
        int i,startIdx,idx;
        int nTpt = (*(E -> nrow));
        startIdx = startTime + startLoc*nTpt;
        idx = startIdx;

        if (startTime == 0)
        {
            (A0 -> E0)[startLoc] = (N[idx] - (A0 -> S0)[startLoc] - (A0 -> I0)[startLoc]
                                    - (A0 -> R0)[startLoc]);
        }
        for (i = startTime; i < *(S->nrow); i++)
        {
            (E -> data)[idx] = N[idx] - (S->data)[idx] - (R->data)[idx] - (I->data)[idx];  
            idx += 1;
        }
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
        calculateGenericCompartment_CPU(I, A0 -> I0,
                                        I_star, R_star); 
    }

    void ModelContext::calculateI_CPU(int startLoc, int startTime)
    { 
        calculateGenericCompartment_CPU(I, A0 -> I0,
                                        I_star, R_star,
                                        startLoc, startTime);
    }

    void ModelContext::calculateI_givenR_CPU()
    {
        int i;
        int nLoc = *(I -> ncol);
        int nTpt = *(I -> nrow);
        int maxItr = nLoc*nTpt;
        for (i = 0; i < nLoc; i++)
        {
            (A0 -> I0)[i] = N[i*nTpt] - (A0 -> S0)[i] - (A0 -> E0)[i] - (A0 -> R0)[i];
        }

        for (i = 0; i < maxItr; i++)
        {
            (I->data)[i] = N[i] - (S->data)[i] - (E->data)[i] - (R->data)[i]; 
        }
    }

    void ModelContext::calculateI_givenR_CPU(int startLoc, int startTime)
    {
        int i,startIdx,idx;
        int nTpt = (*(E -> nrow));
        startIdx = startTime + startLoc*nTpt;
        idx = startIdx;
        if (startTime == 0)
        {
            (A0 -> I0)[startLoc] = (N[idx] - (A0 -> S0)[startLoc] - (A0 -> E0)[startLoc]
                                    - (A0 -> R0)[startLoc]);
        }

        for (i = startTime; i < *(S->nrow); i++)
        {
            (I -> data)[idx] = N[idx] - (S->data)[idx] - (R->data)[idx] - (E->data)[idx];  
            idx += 1;
        }
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

        calculateGenericCompartment_CPU(R, A0 -> R0,
                                        R_star, S_star);

    }

    void ModelContext::calculateR_CPU(int startLoc, int startTime)
    {

        calculateGenericCompartment_CPU(R, A0 -> R0,
                                        R_star, S_star,
                                        startLoc, startTime);

    }

    void ModelContext::calculateR_givenS_CPU()
    {

        int i;
        int nTpt = *(I -> nrow);
        int nLoc = *(I -> ncol);
        int maxItr = nLoc*nTpt;
        for (i = 0; i < nLoc; i++)
        {
            (A0 -> R0)[i] = N[i*nTpt] - (A0 -> S0)[i] - (A0 -> E0)[i] - (A0 -> I0)[i];
        }

        for (i = 0; i < maxItr; i++)
        {
            (R->data)[i] = N[i] - (S->data)[i] - (E->data)[i] - (I->data)[i]; 
        }
    }

    void ModelContext::calculateR_givenS_CPU(int startLoc, int startTime)
    {
        int i,startIdx,idx;
        int nTpt = *(I -> nrow);
        startIdx = startTime + startLoc*nTpt;
        idx = startIdx;
        if (startTime == 0)
        {
            (A0 -> R0)[startLoc] = (N[startLoc*nTpt] - (A0 -> S0)[startLoc] - (A0 -> E0)[startLoc]
                                    - (A0 -> I0)[startLoc]);
        }
        for (i = startTime; i < *(S->nrow); i++)
        {
            (R -> data)[idx] = N[idx] - (S->data)[idx] - (E->data)[idx] - (I->data)[idx];  
            idx += 1;
        }
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
                                                   CompartmentalModelMatrix *compStarSub)
    {
        int i,j,idx1,idxL1;
        int numLoc = *(comp -> ncol);
        int numTpts = *(comp -> nrow);

        for (i = 0; i < numLoc; i++)
        {

            (comp -> data)[i*numTpts] = ((comp0)[i]);
            idx1 = i*numTpts;
            for (j = 1; j < numTpts; j++)
            {
                idxL1 = idx1;
                idx1++;
                (comp -> data)[idx1] = (comp -> data)[idxL1] + 
                                       (compStarAdd -> data)[idxL1] - 
                                       (compStarSub -> data)[idxL1];
            }
        }
    }

    void ModelContext::calculateGenericCompartment_CPU(CompartmentalModelMatrix *comp,int *comp0, 
                                                   CompartmentalModelMatrix *compStarAdd, 
                                                   CompartmentalModelMatrix *compStarSub,
                                                   int startLoc, int startTime)
    {
        int j,idx1,idxL1;
        int startIdx;
        int numTpts = *(comp -> nrow);

        startIdx = startLoc*numTpts + startTime; 

        if (startTime == 0)
        {
            (comp -> data)[startIdx] = ((comp0)[startLoc]);
            idx1 = startLoc*numTpts;
            for (j = 1; j < numTpts; j++)
            {
                idxL1 = idx1;
                idx1++;
                (comp -> data)[idx1] = (comp -> data)[idxL1] + 
                                       (compStarAdd -> data)[idxL1] - 
                                       (compStarSub -> data)[idxL1];
            }
        }
        else
        {
            idx1 = startLoc*numTpts + startTime;
            idxL1 = idx1 - 1; 
            for (j = startTime; j<numTpts; j++)
            {

                
                (comp -> data)[idx1] = (comp -> data)[idxL1] + 
                                       (compStarAdd -> data)[idxL1] - 
                                       (compStarSub -> data)[idxL1];    
                idxL1 = idx1;
                idx1++;
            }
        }
    }
 
    void ModelContext::calculateGenericCompartment_OCL(int *comp,int *comp0, 
                                                   int *compStarAdd, int *compStarSub)
    {
        throw(-1);
    }


    void ModelContext::calculateP_EI_CPU()
    {
        int i;
        int nTpt = *(S -> nrow);
        double gammaVal = *gamma_ei;
        for (i = 0; i < nTpt; i++)
        {
            p_ei[i] = 1-exp(-offset[i]*gammaVal); 
        }
    }

    void ModelContext::calculateP_IR_CPU()
    {
        int i;
        int nTpt = *(S -> nrow);
        double gammaVal = *gamma_ir;
        for (i = 0; i < nTpt; i++)
        {
            p_ir[i] = 1-exp(-offset[i]*gammaVal); 
        }
    }

    // Method: calculateP_RS
    // Accesses: betaPrs, X_pRS
    // Updates: p_rs
    void ModelContext::calculateP_RS_CPU()
    {
        int i;
        int neta = *(X_pRS -> nrow_x);
        X_pRS -> calculate_eta_CPU(p_rs, betaPrs);
        for (i = 0; i < neta; i++)
        {
            p_rs[i] = 1-exp(-offset[i]*exp(p_rs[i]));
        }
    }

    // Method: calculatePi
    // Accesses: beta, I, N, distMat, rho
    // Updates: p_se
    void ModelContext::calculateP_SE_CPU()
    {
        this -> cacheP_SE_Calculation(); 
        unsigned int i, j, index;

        // Calculate dmu: I/N * exp(eta)
        unsigned int nLoc = (unsigned int) *(S -> ncol);
        unsigned int nTpt = (unsigned int) *(S -> nrow);

        memset(p_se, 0, nLoc*nTpt*sizeof(double));
        // Calculate rho*sqrt(idmat)
        
        // C := alpha * op( A ) %*% op ( B ) + beta*C
        // op( X ) is one of op( X ) = X, op( X ) == X**T
        // A, B, and C are matrices with
        // op( A ) a m by k matrix, 
        // op( B ) a k by n matrix,
        // C an m by n matrix
        double* DM;
        for (i = 0; i < (scaledDistMatrices -> size()); i++)
        {
            DM = (*scaledDistMatrices)[i] -> data;
#ifdef LSS_USE_BLAS
            cblas_dgemm(CblasColMajor,      // order
                        CblasNoTrans,       // TransA
                        CblasNoTrans,       // TransB
                        *(I->nrow),         // M
                        *(I->ncol),         // N
                        *(I->ncol),         // K
                        rho[i],             // alpha
                        p_se_components,    // A 
                        *(I->nrow),         // ldA
                        DM,                 // B 
                        *(I-> ncol),        // ldB
                        1.0,                // beta
                        p_se,               // C
                        *(I->nrow));        // ldC 
#else
			// Implement Eigen here
            MatrixMapType Amap(p_se_components, *(I->nrow), *(I->ncol));
            MatrixMapType Bmap(DM, *(I->ncol), *(I->ncol));
            MatrixMapType seMap(p_se, *(I->nrow), *(I->ncol));
            seMap += rho[i]*Amap*Bmap; 
#endif			
        }

        for (i = 0; i < nLoc; i++) 
        {
            index = i*nTpt;
            for (j = 0; j < nTpt; j++)
            {
                p_se[index] = 1-exp(-offset[j]*(p_se_components[index] + p_se[index]));
                index++;
            }
        }        
    }

    // To be used when beta is fixed, eta has already been exponentiated,
    // only change is to I compartment. 
    void ModelContext::calculateP_SE_CPU(int startLoc, int startTime)
    {

        unsigned int i, j, index;

        // Calculate dmu: I/N * exp(eta)
        unsigned int nLoc = (unsigned int) *(S -> ncol);
        unsigned int nTpt = (unsigned int) *(S -> nrow);

        for (i = 0; i < nLoc; i++)
        {
            index = i*nTpt + startTime;
            for (j = startTime; j < nTpt; j++)
            {
                p_se[index] = 0.0;
                index++;
            }
        }

        index = startLoc*nTpt + startTime;
        for (j = startTime; j < nTpt; j++)
        {
            p_se_components[index] = 
               ((I -> data)[index] * (eta[index]))/N[index];
            index++;
        }

        double* DM;
        for (i = 0; i < (scaledDistMatrices -> size()); i++)
        {
            DM = (*scaledDistMatrices)[i] -> data;
#ifdef LSS_USE_BLAS
            cblas_dgemm(CblasColMajor,         // order
                        CblasNoTrans,          // TransA
                        CblasNoTrans,          // TransB
                        (*(I->nrow)-startTime),// M
                        *(I->ncol),            // N
                        *(I->ncol),            // K
                        rho[i],                // alpha
                        &(p_se_components[startTime]),       // A 
                        *(I->nrow),            // ldA
                        DM,                    // B 
                        *(I-> ncol),           // ldB
                        1.0,                   // beta
                        &(p_se[startTime]),                  // C
                        *(I->nrow));           // ldC 

#else
            Eigen::Map<MatrixType, Eigen::ColMajor, Eigen::OuterStride<> > 
			Amap(&(p_se_components[startTime]), (*(I -> nrow) - startTime), *(I->ncol), Eigen::OuterStride<>((*(I->nrow))));
            //MatrixMapType Amap(p_se_components, *(I->nrow), *(I->ncol));
            MatrixMapType Bmap(DM, *(I->ncol), *(I->ncol));
			
            Eigen::Map<MatrixType, Eigen::ColMajor, Eigen::OuterStride<> > 
			seMap(&(p_se[startTime]), (*(I -> nrow) - startTime), *(I->ncol), Eigen::OuterStride<>((*(I->nrow))));
            //MatrixMapType seMap(p_se, *(I->nrow), *(I->ncol));
            seMap += rho[i]*Amap*Bmap; 
#endif

        }
        for (i = 0; i < nLoc; i++) 
        {
            index = i*nTpt + startTime;
            for (j = startTime; j < nTpt; j++)
            {
                p_se[index] = 1-exp(-offset[j]*(p_se_components[index] + p_se[index]));
                index++;
            }
        }        
    }

    void ModelContext::cacheP_SE_Calculation()
    {
        int i, j, index;
        //Update Eta
        this -> X -> calculate_eta_CPU(eta, beta);

        //Exponentiate
        int nrowx = *(X->nrow_x);
        for (i = 0; i < nrowx; i++)
        {
            eta[i] = std::exp(eta[i]);
        }

        // Calculate dmu: I/N * exp(eta)
        int nLoc = *(S -> ncol);
        int nTpt = *(S -> nrow);

        for (i = 0; i < nLoc; i++) 
        {
            index = i*nTpt;
            for (j = 0; j < nTpt; j++)
            {

                p_se_components[index] = 
                   ((I -> data)[index] * (eta[index]))/N[index];
                index++;
            }
        }
    }

    void ModelContext::calculateP_SE_OCL()
    {
        oclProvider -> calculateP_SE(this); 
    }

    ModelContext::~ModelContext()
    {
        delete numIterations;
        if (*isPopulated)
        {
            delete[] Y;
            delete[] offset;
            delete S0_fc;
            delete E0_fc;
            delete I0_fc;
            delete R0_fc;
            delete S_star_fc;
            delete E_star_fc;
            delete R_star_fc;
            delete beta_fc;
            delete rho_fc;
            delete phi_fc;
            delete I_star_overdispersed_fc;
            delete betaPrs_fc;
            delete gamma_ei_fc;
            delete p_se_fc;
            delete gamma_ir_fc;
            delete gamma_ei;
            delete gamma_ir;
            delete tmpContainer;
            delete[] p_se;
            delete[] p_se_components;
            delete[] compartmentCache;
            delete[] p_ei;
            delete[] p_ir;
            delete[] p_rs;
            delete config;
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
            delete X_pRS;
            delete[] betaPrs;
            delete[] N;
            delete[] beta;
            delete[] eta;
            delete[] rho;
            delete phi;
            delete[] gamma;
            delete singleLocation;   
            // FC's have already been disposed of.
            delete model;
            while (iterationTasks -> size() != 0){(iterationTasks -> pop_back());}
            delete setSamplingIndicesTask;
            delete performHybridSE_EI_UpdateTask;
            delete performHybridIR_RS_UpdateTask;
            delete decorrelationStepTask;
            while (scaledDistMatrices -> size() != 0){delete (*scaledDistMatrices).back(); (scaledDistMatrices -> pop_back());}
            delete scaledDistMatrices;
            delete iterationTasks;
            delete oclProvider;
            delete[] indexList;
            delete indexLength;

        }
        delete isPopulated;
    }
}



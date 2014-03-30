#include<iostream>
#include<stdio.h>
#include<math.h>
#include<cstring>
#include<vector>
#include<cblas.h>
#include<cmath>
#include<algorithm>
#include<FullConditional.hpp>
#include<ModelContext.hpp>
#include<OCLProvider.hpp>
#include<CompartmentalModelMatrix.hpp>
#include<CovariateMatrix.hpp>
#include<RandomNumberProvider.hpp>

namespace SpatialSEIR
{
    using std::cout;
    using std::endl;

    /*
     * Helper functions
     */
    int matMult(double* output, double * A, double * B, int Arow, int Acol, int Brow, int Bcol, bool TransA = false, bool TransB = false)
    {
        // Use BLAS to matrix multiply, assume column major, non transposing.
        // Double check this code when I have internet access. 
        cblas_dgemm(CblasColMajor,
                    TransA ? CblasTrans : CblasNoTrans,
                    TransB ? CblasTrans : CblasNoTrans,
                    Arow, Bcol, Brow,
                    1.0, 
                    A, Arow, 
                    B, Brow, 
                    0.0, output, Brow);
        return 0; 
    }

    //Log scale
    double dbeta(double x, double a, double b)
    {
        double out = (a-1)*std::log(x) + 
            (b-1)*std::log(1-x) + 
            (lgamma(a+b)) - 
            ((lgamma(a)) + (lgamma(b)));
        return(out);
    }

    /*
     *
     * Implement the data container class InitData
     *
     */    
 
 
    InitData::InitData(int *_S0, 
                       int *_E0,
                       int *_I0,
                       int *_R0,
                       int *_S_star0,
                       int *_E_star0, 
                       int *_I_star0, 
                       int *_R_star0,
                       int *nLoc)
    {
        this -> populate(*&_S0, *&_E0, *&_I0, *&_R0, 
                *&_S_star0, *&_E_star0, *&_I_star0, *&_R_star0,
                *&nLoc);
    }

    InitData::InitData()
    {
        // Do nothing
    }

    void InitData::populate(int *_S0, 
                       int *_E0,
                       int *_I0,
                       int *_R0,
                       int *_S_star0,
                       int *_E_star0, 
                       int *_I_star0, 
                       int *_R_star0,
                       int *nLoc
                       )
    {
        S0 = new int[*nLoc];
        E0 = new int[*nLoc];
        I0 = new int[*nLoc];
        R0 = new int[*nLoc];
        S_star0 = new int[*nLoc]; 
        E_star0 = new int[*nLoc];
        I_star0 = new int[*nLoc];
        R_star0 = new int[*nLoc];
        numLocations = new int;
        *numLocations = *nLoc;
        int i;
        for (i = 0; i < *nLoc; i++)
        {
            S0[i] = _S0[i]; 
            E0[i] = _E0[i];
            I0[i] = _I0[i];
            R0[i] = _R0[i];
            S_star0[i] = _S_star0[i];
            E_star0[i] = _E_star0[i];
            I_star0[i] = _I_star0[i];
            R_star0[i] = _R_star0[i];
        }
    }

    InitData::~InitData()
    {
        delete S0;
        delete E0;
        delete I0;
        delete R0;
        delete S_star0;
        delete E_star0;
        delete I_star0;
        delete R_star0;
        delete numLocations;
    }

    int FullConditional::sampleCompartment(ModelContext* context,
                                       InitData* A0,
                                       CompartmentalModelMatrix* drawCompartment,
                                       CompartmentalModelMatrix* destCompartment,
                                       CompartmentalModelMatrix* starCompartment,
                                       double width)
    {
        // Declare required variables
        int i, j, compIdx;
        int nLoc = *(A0 -> numLocations);
        int nTpts = *(drawCompartment -> ncol);
        double l,r,y,x,x0;
        
        // Update the relevant CompartmentalModelMatrix instances
        this -> calculateRelevantCompartments();

        // Allocate storage for the cached components of
        // the full conditional calculation
        double* cachedValues = new double[nLoc*nTpts];
        this -> cacheEvalCalculation(cachedValues);

        // Set the "value" attribute appropriately
        this -> evalCPU();
   
        // Main loop: 
        for (j = 0; j < nTpts; j ++)
        { 
            compIdx = j*nLoc - 1;
            for (i = 0; i < nLoc; i++)
            {
                compIdx++;
                x = (starCompartment -> data)[compIdx];
                this -> calculateRelevantCompartments(i,j); 
                this -> evalCPU(i,j,cachedValues);
                y = (this->getValue()) - (context -> random -> gamma());
                l = 0.0;
                r = l + width;

                while (y >= (this -> getValue()))
                {
                    x0 = std::floor((context -> random -> uniform())*(r));
                    (starCompartment -> data)[compIdx] = x0;
                    this -> calculateRelevantCompartments(i,j);
                    this -> evalCPU(i,j,cachedValues);
                    r = (x0 < x ? r : x0); 
                }
            }
        }

        delete[] cachedValues;
        return 0;
    }

    int FullConditional::sampleDouble(ModelContext* context,
                                       InitData* A0,
                                       double* variable, 
                                       int varLen, 
                                       double width)
    {
        // Declare required variables
        int i;
        double l,r,y,x,x0;
        
        // Update the relevant CompartmentalModelMatrix instances
        this -> calculateRelevantCompartments();

        // Set the "value" attribute appropriately
        this -> evalCPU();
   
        // Main loop: 
        for (i = 0; i < varLen; i++)
        { 
            x = variable[i];
            this -> calculateRelevantCompartments(); 
            this -> evalCPU();
            y = (this->getValue()) - (context -> random -> gamma());
            l = 0.0;
            r = l + width;

            while (y >= (this -> getValue()))
            {
                x0 = (context -> random -> uniform())*(r);
                variable[i] = x0;
                this -> calculateRelevantCompartments();
                this -> evalCPU();
                r = (x0 < x ? r : x0);  
            }
        }
        return 0;

    }


    /*
     *
     * Implement the full conditional distribution for S_star
     *
     */    

    FC_S_Star::FC_S_Star(ModelContext * _context,
                         CompartmentalModelMatrix *_S_star, 
                         CompartmentalModelMatrix *_S, 
                         CompartmentalModelMatrix *_R, 
                         InitData *_A0,
                         CovariateMatrix *_X, 
                         double *_p_se,
                         double *_p_rs,
                         double *_beta,
                         double *_rho)
    {
       context = new ModelContext*;
       S_star = new CompartmentalModelMatrix*;
       S = new CompartmentalModelMatrix*;
       R = new CompartmentalModelMatrix*;
       A0 = new InitData*;
       X = new CovariateMatrix*;
       p_se = new double*;
       p_rs = new double*;
       beta = new double*;
       rho = new double*;
       value = new double;
       *context = _context;
       *S_star = _S_star;
       *S = _S;
       *R = _R;
       *A0 = _A0;
       *X = _X;
       *p_se = _p_se;
       *p_rs = _p_rs;
       *beta = _beta;
       *rho = _rho;
       *value = -1.0;
    }    
    FC_S_Star::~FC_S_Star()
    {
        delete S_star;
        delete S;
        delete R;
        delete A0;
        delete X;
        delete p_se;
        delete p_rs;
        delete beta;
        delete rho;
        delete value;
        delete context;
    }

    // Cache the components of the FC_S_Star calculation in cachedValues
    int FC_S_Star::cacheEvalCalculation(double* cachedValues)
    {
        int i, j, tmp, compIdx;
        int nLoc = *((*A0) -> numLocations);
        int nTpts = *((*S) -> ncol);
        for (j = 0; j < nTpts; j++)     
        {
            compIdx = j*nLoc - 1;
            for (i = 0; i < nLoc; i++)    
            {
                compIdx += 1;
                tmp = ((*S_star) -> data)[compIdx];
                cachedValues[compIdx] = (std::log((*p_rs)[j])*tmp +  
                                   std::log(1-(*p_rs)[j])*(((*R) -> data)[compIdx] - tmp) +
                                   std::log(1-(*p_se)[compIdx])*(((*S) -> data)[compIdx]));
            }
        } 
        return(0);
    }
    // Evaluate the S_star FC at the current values provided by the context.
    int FC_S_Star::evalCPU()
    {
        *value = 0.0;
        int i, j, tmp, compIdx;
        int nLoc = *((*A0) -> numLocations);
        int nTpts = *((*S) -> ncol);
        double term1, term2, term3;
        term1 = 0.0; term2 = 0.0; term3 = 0.0;
        compIdx = 0;
        for (j = 0; j < nTpts; j++)     
        {
            compIdx = j*nLoc - 1;
            for (i = 0; i < nLoc; i++)    
            {
                compIdx += 1;
                tmp = ((*S_star) -> data)[compIdx];
                // We're only getting non-negative values
                //if (tmp < 0)
                //{
                //    *value = -std::log(0.0);
                //    return(-1);
                //}
                term1 += std::log((*p_rs)[j])*tmp; 
                term2 += std::log(1-(*p_rs)[j])*(((*R) -> data)[compIdx] - tmp);
                term3 += std::log(1-(*p_se)[compIdx])*(((*S) -> data)[compIdx]) ;
            }
        } 
        *value = term1 + term2 + term3;
        return(0);
    }
    int FC_S_Star::evalCPU(int startLoc, int startTime, double* cachedValues)
    {
        *value = 0.0;
        int i, j, tmp, compIdx;
        int nLoc = *((*A0) -> numLocations);
        int nTpts = *((*S) -> ncol);

        compIdx = startLoc + startTime*nLoc;
        for (j = startTime; j < nTpts; j++)
        {
            tmp = ((*S_star) -> data)[compIdx];
            cachedValues[compIdx] = (std::log((*p_rs)[j])*tmp +  
                               std::log(1-(*p_rs)[j])*(((*R) -> data)[compIdx] - tmp)+
                               std::log(1-(*p_se)[compIdx])*(((*S) -> data)[compIdx]));
            compIdx += nLoc;
        }
        for (j = 0; j < nTpts; j++)     
        {
            compIdx = j*nLoc - 1;
            for (i = 0; i < nLoc; i++)    
            {
                compIdx += 1;
                *value += cachedValues[compIdx]; 
            }
        } 
        return(0);
    }

    int FC_S_Star::evalOCL()
    {
        //NOT IMPLEMENTED
        return-1;
    }
    int FC_S_Star::calculateRelevantCompartments()
    {
        (*context) -> calculateS_CPU();
        (*context) -> calculateR_CPU();
    }
    int FC_S_Star::calculateRelevantCompartments(int startLoc, int startTime)
    {
        (*context) -> calculateS_CPU(startLoc, startTime);
        (*context) -> calculateR_CPU(startLoc, startTime);
    }

    int FC_S_Star::sampleCPU()
    {
        this -> sampleCompartment(*context,*A0,*R,*S,
                                  *S_star,10);
        return 0;
    }
    int FC_S_Star::sampleOCL()
    {
        //NOT IMPLEMENTED
        return -1;
    }
    double FC_S_Star::getValue()
    {
        return(*(this -> value));
    }



    /*
     *
     * Implement the full conditional distribution for E_Star
     *
     */    

    
    FC_E_Star::FC_E_Star(ModelContext *_context,
                         CompartmentalModelMatrix *_E_star,
                         CompartmentalModelMatrix *_E,  
                         CompartmentalModelMatrix *_S,
                         CovariateMatrix *_X,
                         InitData *_A0,
                         double *_p_se,
                         double *_p_ei,
                         double *_rho,
                         double *_beta) 
    {

        context = new ModelContext*;
        E_star = new CompartmentalModelMatrix*;
        E = new CompartmentalModelMatrix*;
        S = new CompartmentalModelMatrix*;
        X = new CovariateMatrix*;
        A0 = new InitData*;
        p_se = new double*;
        p_ei = new double*;
        rho = new double*;
        beta = new double*;
        value = new double;
        
       *context = _context;
        *E_star = _E_star;
        *E = _E;
        *S = _S;
        *X = _X;
        *A0 = _A0;
        *p_se = _p_se;
        *p_ei = _p_ei;
        *rho = _rho;
        *beta = _beta;
        *value = -1.0;
    }

    FC_E_Star::~FC_E_Star()
    {
        delete E_star;
        delete E;
        delete S;
        delete X;
        delete A0;
        delete p_se;
        delete p_ei;
        delete rho;
        delete beta;
        delete value;
        delete context;
    }

    int FC_E_Star::cacheEvalCalculation(double* cachedValues)
    {
        int i, j, tmp, compIdx;
        int nLoc = *((*A0) -> numLocations);
        int nTpts = *((*S) -> ncol);
    
        for (j = 0; j < nTpts; j++)     
        {
            compIdx = j*nLoc - 1;
            for (i = 0; i < nLoc; i++)    
            {
                compIdx += 1;
                tmp = ((*E_star) -> data)[compIdx];
                cachedValues[compIdx] = (std::log((*p_se)[compIdx])*tmp +
                                         std::log(1-(*p_se)[compIdx])*(((*S) -> data)[compIdx] - tmp) +
                                         std::log(1-(**p_ei))*(((*E) -> data)[compIdx]));
            }
        } 
        return(0);
    }

    int FC_E_Star::evalCPU()
    {
        *value = 0.0;
        int i, j, tmp, compIdx;
        int nLoc = *((*A0) -> numLocations);
        int nTpts = *((*S) -> ncol);
        double term1, term2, term3;
        term1 = 0.0; term2 = 0.0; term3 = 0.0;
        for (j = 0; j < nTpts; j++)     
        {
            for (i = 0; i < nLoc; i++)    
            {
                compIdx = i + j*nLoc;
                tmp = ((*E_star) -> data)[compIdx];
                if (tmp < 0)
                {
                    *value = -std::log(0.0);
                    return(-1);
                }
                term1 += std::log((*p_se)[compIdx])*tmp; 
                term2 += std::log(1-(*p_se)[compIdx])*(((*S) -> data)[compIdx] - tmp);
                term3 += std::log(1-(**p_ei))*(((*E) -> data)[compIdx]) ;
            }
        } 
        *value = term1 + term2 + term3;
        return(0);
    }
    int FC_E_Star::evalCPU(int startLoc, int startTime, double* cachedValues)
    {
        *value = 0.0;
        int i, j, tmp, compIdx;
        int nLoc = *((*A0) -> numLocations);
        int nTpts = *((*S) -> ncol);
     
        compIdx = startLoc + startTime*nLoc;
        for (j = startTime; j < nTpts; j++)
        {
            tmp = ((*E_star) -> data)[compIdx];
            cachedValues[compIdx] = (std::log((*p_se)[compIdx])*tmp +
                                     std::log(1-(*p_se)[compIdx])*(((*S) -> data)[compIdx] - tmp) +
                                     std::log(1-(**p_ei))*(((*E) -> data)[compIdx]));
            compIdx += nLoc;
        }
        for (j = 0; j < nTpts; j++)     
        {
            compIdx = j*nLoc - 1;
            for (i = 0; i < nLoc; i++)    
            {
                compIdx += 1;
                *value += cachedValues[compIdx]; 
            }
        } 
        return(0);
    }

    int FC_E_Star::evalOCL()
    {
        //NOT IMPLEMENTED
        return -1;
    }
    int FC_E_Star::calculateRelevantCompartments()
    {
        (*context) -> calculateS_CPU();
        (*context) -> calculateE_CPU();
    }
    int FC_E_Star::calculateRelevantCompartments(int startLoc, int startTime)
    {
        (*context) -> calculateS_CPU(startLoc, startTime);
        (*context) -> calculateE_CPU(startLoc, startTime);
    }
    int FC_E_Star::sampleCPU()
    {
        this -> sampleCompartment(*context,*A0,*S,*E,
                                  *E_star,10);
        return 0;
    }
    int FC_E_Star::sampleOCL()
    {
        //NOT IMPLEMENTED
        return -1;
    }

    double FC_E_Star::getValue()
    {
        return(*(this -> value));
    }

    /*
     *
     * Implement the full conditional distribution for R_Star
     *
     */    
    FC_R_Star::FC_R_Star(ModelContext *_context,
                         CompartmentalModelMatrix *_R_star,
                         CompartmentalModelMatrix *_R,
                         CompartmentalModelMatrix *_I,
                         InitData *_A0,
                         double *_p_rs,
                         double *_p_ir)
    {

        context = new ModelContext*;
        R_star = new CompartmentalModelMatrix*;
        R = new CompartmentalModelMatrix*;
        I = new CompartmentalModelMatrix*;
        A0 = new InitData*;
        p_rs = new double*;
        p_ir = new double*;
        value = new double;

       *context = _context;
        *R_star = _R_star;
        *R = _R;
        *I = _I;
        *A0 = _A0;
        *p_rs = _p_rs;
        *p_ir = _p_ir;
        *value = -1.0;
    }
    FC_R_Star::~FC_R_Star()
    {
        delete R_star;
        delete R;
        delete I;
        delete A0;
        delete p_rs;
        delete p_ir;
        delete value;
        delete context;
    }

    int FC_R_Star::cacheEvalCalculation(double* cachedValues)
    {
        int i, j, tmp, compIdx;
        int nLoc = *((*A0) -> numLocations);
        int nTpts = *((*R) -> ncol);
    
        for (j = 0; j < nTpts; j++)     
        {
            compIdx = j*nLoc - 1;
            for (i = 0; i < nLoc; i++)    
            {
                compIdx += 1;
                tmp = ((*R_star) -> data)[compIdx];
                cachedValues[compIdx] = (std::log((**p_ir))*tmp + 
                                std::log(1-(**p_ir))*(((*I) -> data)[compIdx] - tmp) +
                                std::log(1-((*p_rs)[j]))*(((*R) -> data)[compIdx])) ;
            }
        } 
        return(0);
    }


    int FC_R_Star::evalCPU()
    {
        *value = 0.0;
        int i, j, tmp, compIdx;
        int nLoc = *((*A0) -> numLocations);
        int nTpts = *((*I) -> ncol);
        double term1, term2, term3;
        term1 = 0.0; term2 = 0.0; term3 = 0.0;
        for (j = 0; j < nTpts; j++)     
        {
            for (i = 0; i < nLoc; i++)    
            {
                compIdx = i + j*nLoc;
                tmp = ((*R_star) -> data)[compIdx];
                if (tmp < 0)
                {
                    *value = -std::log(0.0);
                    return(-1);
                }
                term1 += std::log((**p_ir))*tmp; 
                term2 += std::log(1-(**p_ir))*(((*I) -> data)[compIdx] - tmp);
                term3 += std::log(1-((*p_rs)[j]))*(((*R) -> data)[compIdx]) ;
            }
        } 
        *value = term1 + term2 + term3;
        return(0);
    }

    int FC_R_Star::evalCPU(int startLoc, int startTime, double* cachedValues)
    {
        *value = 0.0;
        int i, j, tmp, compIdx;
        int nLoc = *((*A0) -> numLocations);
        int nTpts = *((*R) -> ncol);
     
        compIdx = startLoc + startTime*nLoc;
        for (j = startTime; j < nTpts; j++)
        {
            tmp = ((*R_star) -> data)[compIdx];
            cachedValues[compIdx] = (std::log((**p_ir))*tmp + 
                                std::log(1-(**p_ir))*(((*I) -> data)[compIdx] - tmp) +
                                std::log(1-((*p_rs)[j]))*(((*R) -> data)[compIdx]));
            compIdx += nLoc;
        }
        for (j = 0; j < nTpts; j++)     
        {
            compIdx = j*nLoc - 1;
            for (i = 0; i < nLoc; i++)    
            {
                compIdx += 1;
                *value += cachedValues[compIdx]; 
            }
        } 
        return(0);
    }

    int FC_R_Star::evalOCL()
    {
        //NOT IMPLEMENTED
        return -1;
    }
    int FC_R_Star::calculateRelevantCompartments()
    {
        (*context) -> calculateI_CPU();
        (*context) -> calculateR_CPU();
    }
    int FC_R_Star::calculateRelevantCompartments(int startLoc, int startTime)
    {
        (*context) -> calculateI_CPU(startLoc, startTime);
        (*context) -> calculateR_CPU(startLoc, startTime);
    }

    int FC_R_Star::sampleCPU()
    {
        this -> sampleCompartment(*context,*A0,*I,*R,
                                  *R_star,10);
        return(0);
    }
    int FC_R_Star::sampleOCL()
    {
        //NOT IMPLEMENTED
        return -1;
    }

    double FC_R_Star::getValue()
    {
        return(*(this -> value));
    }

     /*
     *
     * Implement the full conditional distribution for the regression
     * parameters: beta
     */


    FC_Beta::FC_Beta(ModelContext *_context,
                     CompartmentalModelMatrix *_E_star, 
                     CompartmentalModelMatrix *_S, 
                     InitData *_A0,
                     CovariateMatrix *_X,
                     double *_p_se, 
                     double *_beta, 
                     double *_rho)
    {

        context = new ModelContext*;
        E_star = new CompartmentalModelMatrix*;
        S = new CompartmentalModelMatrix*;
        A0 = new InitData*;
        X = new CovariateMatrix*;
        p_se = new double*;
        beta = new double*;
        rho = new double*;
        value = new double;

        *context = _context;
        *E_star = _E_star;
        *S = _S;
        *A0 = _A0;
        *X = _X;
        *p_se = _p_se;
        *beta = _beta;
        *rho = _rho;
        *value = -1.0;
    }

    FC_Beta::~FC_Beta()
    {
        delete E_star;
        delete S;
        delete A0;
        delete X;
        delete p_se;
        delete beta;
        delete rho;
        delete value;
        delete context;
    }
    
    int FC_Beta::cacheEvalCalculation(double* cachedValues)
    {
        //Not Implemented
        throw(-1);
    }


    int FC_Beta::evalCPU()
    {
        *value = 0.0;
        int i, j, tmp, compIdx;
        int nLoc = *((*A0) -> numLocations);
        int nTpts = *((*S) -> ncol);
        double term1, term2, term3;
        term1 = 0.0; term2 = 0.0; term3 = 0.0;
        for (j = 0; j < nTpts; j++)     
        {
            for (i = 0; i < nLoc; i++)    
            {
                compIdx = i + j*nLoc;
                tmp = ((*E_star) -> data)[compIdx];
                term1 += std::log((*p_se)[compIdx])*tmp; 
                term2 += std::log(1-(*p_se)[compIdx])*(((*S) -> data)[compIdx] - tmp);
            }
        } 
        for (i = 0; i < (*((*X) -> ncol_x) + *((*X) -> ncol_z)); i++)
        {
            term3 += pow((*beta)[i],2)/10; // Generalize to allow different prior precisions. 
        }
        *value = term1 + term2 + term3;
        return(0);
    }
    int FC_Beta::evalCPU(int startLoc, int startTime, double* cachedValues)
    {
        //NOT IMPLEMENTED
        throw(-1);
    }

    int FC_Beta::evalOCL()
    {
        //NOT IMPLEMENTED
        return -1;
    }
    int FC_Beta::calculateRelevantCompartments()
    {
        ((*context) -> calculateP_SE_CPU());
        return(0);

    }
    int FC_Beta::calculateRelevantCompartments(int startLoc, int startTime)
    {
        //NOT IMPLEMENTED
        throw(-1);
    }

    int FC_Beta::sampleCPU()
    {
        sampleDouble(*context, *A0, *beta, (*((*X) -> ncol_x) + *((*X) -> ncol_z)), 10.0); 
        return(0);
    }
    int FC_Beta::sampleOCL()
    {
        //NOT IMPLEMENTED
        return -1;
    }

    double FC_Beta::getValue()
    {
        return(*(this -> value));
    }
    /*
     *
     * Implement the full conditional for the R->S transition 
     * probabilities. 
     *
     */
    FC_P_RS::FC_P_RS(ModelContext *_context,
                     CompartmentalModelMatrix *_S_star, 
                     CompartmentalModelMatrix *_R,
                     InitData *_A0,
                     double *_p_rs)
    {

        context = new ModelContext*;
        S_star = new CompartmentalModelMatrix*;
        R = new CompartmentalModelMatrix*;
        A0 = new InitData*;
        p_rs = new double*;
        value = new double;

        *context = _context;
        *S_star = _S_star;
        *R = _R;
        *A0 = _A0;
        *p_rs = _p_rs;
        *value = -1.0;
    }
    FC_P_RS::~FC_P_RS()
    {
        delete S_star;
        delete R;
        delete A0;
        delete p_rs;
        delete value;
        delete context;
    }
    int FC_P_RS::cacheEvalCalculation(double* cachedValues)
    {
        //Not Implemented
        throw(-1);
    }


    int FC_P_RS::evalCPU()
    {
        *value = 0.0;
        int i, j, tmp, compIdx;
        int nLoc = *((*A0) -> numLocations);
        int nTpts = *((*R) -> ncol);
        double term1, term2, term3;
        term1 = 0.0; term2 = 0.0; term3 = 0.0;
        int* s_star_i_sum = new int[nTpts];
        int* r_star_i_diff = new int[nTpts];
        for (j =0; j<nTpts; j++)
        {
            s_star_i_sum[j]=0.0;
            r_star_i_diff[j]=0.0;
        }


        for (j = 0; j < nTpts; j++)     
        {
            for (i = 0; i < nLoc; i++)    
            {
               compIdx = i + j*nLoc;
                tmp = ((*S_star) -> data)[compIdx];
                s_star_i_sum[j] += tmp; 
                r_star_i_diff[j] += (((*R) -> data)[compIdx] - tmp);
            }
        } 
        *value = 0;
        for (j = 0; j < nTpts; j++)
        {
            *value += dbeta((*p_rs)[j],1.5 + s_star_i_sum[j], 1.5 + r_star_i_diff[j]);
        }
       
        delete[] s_star_i_sum;
        delete[] r_star_i_diff;
        return 0;
    }
    int FC_P_RS::evalCPU(int startLoc, int startTime, double* cachedValues)
    {
        //NOT IMPLEMENTED
        throw(-1);
    }

    int FC_P_RS::evalOCL()
    {
        //NOT IMPLEMENTED
        return -1;
    }
    int FC_P_RS::calculateRelevantCompartments()
    {
        // Not used. Do nothing
        return(0);
    }
    int FC_P_RS::calculateRelevantCompartments(int startLoc, int startTime)
    {
        //Not used. Do nothing. 
        return(0);
    }

    int FC_P_RS::sampleCPU()
    {
        sampleDouble(*context, *A0, *p_rs, *((*R)->ncol), 10.0); 
        return(0);
    }
    int FC_P_RS::sampleOCL()
    {
        //NOT IMPLEMENTED
        return -1;
    }

    double FC_P_RS::getValue()
    {
        return(*(this -> value));
    }


    FC_Rho::FC_Rho(ModelContext *_context,
                   CompartmentalModelMatrix *_E_star,  
                   CompartmentalModelMatrix *_S,
                   InitData *_A0,
                   CovariateMatrix *_X,
                   double *_p_se,
                   double *_beta,
                   double *_rho)
    {
        context = new ModelContext*;
        E_star = new CompartmentalModelMatrix*;
        S = new CompartmentalModelMatrix*;
        A0 = new InitData*;
        X = new CovariateMatrix*;
        p_se = new double*;
        beta = new double*;
        rho = new double*;
        value = new double;

        *context = _context;
        *E_star = _E_star;
        *S = _S;
        *A0 = _A0;
        *X = _X;
        *p_se = _p_se;
        *beta = _beta;
        *rho = _rho;
        *value = -1.0;
    }
    FC_Rho::~FC_Rho()
    {
        delete E_star;
        delete S;
        delete A0;
        delete X;
        delete p_se;
        delete beta;
        delete rho;
        delete value;
        delete context;
    }
    int FC_Rho::cacheEvalCalculation(double* cachedValues)
    {
        //Not Implemented
        throw(-1);
    }


    int FC_Rho::evalCPU()
    {
        *value = 0.0;
        int i, j, tmp, compIdx;
        int nLoc = *((*A0) -> numLocations);
        int nTpts = *((*S) -> ncol);
        double term1, term2, term3;
        term1 = 0.0; term2 = 0.0; term3 = 0.0;
        for (j = 0; j < nTpts; j++)     
        {
            for (i = 0; i < nLoc; i++)    
            {
                compIdx = i + j*nLoc;
                tmp = ((*E_star) -> data)[compIdx];
                term1 += std::log((*p_se)[compIdx])*tmp; 
                term2 += std::log(1-(*p_se)[compIdx])*(((*S) -> data)[compIdx] - tmp);
            }
        } 
        term3 += (**rho >0 && **rho < 1 ? 0 : -INFINITY); // Generalize to allow informative priors. 
                                                        // Prior specification in this area needs work. 
        *value = term1 + term2 + term3;
        return(0);
    }
    int FC_Rho::evalCPU(int startLoc, int startTime, double* cachedValues)
    {
        //NOT IMPLEMENTED
        throw(-1);
    }

    int FC_Rho::evalOCL()
    {
        //NOT IMPLEMENTED
        return -1;
    }
    int FC_Rho::calculateRelevantCompartments()
    {
       (*context) -> calculateP_SE_CPU();
       return(0); 
    }
    int FC_Rho::calculateRelevantCompartments(int startLoc, int startTime)
    {
        //NOT IMPLEMENTED
        throw(-1);
    }

    int FC_Rho::sampleCPU()
    {
        sampleDouble(*context, *A0, *rho, 1, 0.5); 
        return(0);
    }
    int FC_Rho::sampleOCL()
    {
        //NOT IMPLEMENTED
        return -1;
    }

    double FC_Rho::getValue()
    {
        return(*(this -> value));
    }


    FC_P_EI::FC_P_EI(ModelContext *_context,
                     CompartmentalModelMatrix *_I_star,
                     CompartmentalModelMatrix *_E,
                     InitData *_A0,
                     double *_p_ei)
    {

        context = new ModelContext*;
        I_star = new CompartmentalModelMatrix*;
        E = new CompartmentalModelMatrix*;
        A0 = new InitData*;
        p_ei = new double*;
        value = new double;

        *context = _context;
        *I_star = _I_star;
        *E = _E;
        *A0 = _A0;
        *p_ei = _p_ei;
        *value = -1.0;

    }
    FC_P_EI::~FC_P_EI()
    {
        delete I_star;
        delete E;
        delete A0;
        delete p_ei;
        delete value;
        delete context;
    }
    int FC_P_EI::cacheEvalCalculation(double* cachedValues)
    {
        //Not Implemented
        throw(-1);
    }


    int FC_P_EI::evalCPU()
    { 
        *value = 0.0;
        int i, j, tmp, compIdx;
        int nLoc = *((*A0) -> numLocations);
        int nTpts = *((*E) -> ncol);
        double term1, term2, term3;
        term1 = 0.0; term2 = 0.0; term3 = 0.0;
        int i_star_sum = 0;
        int e_sum = 0;
        for (j =0; j<(nLoc*nTpts); j++)
        {
            i_star_sum += ((*I_star)-> data)[j];
            e_sum += ((*E)-> data)[j];
        }

        *value = dbeta(**p_ei, 1.5 + i_star_sum, 1.5 - i_star_sum + e_sum); 
        return 0;
    }
    int FC_P_EI::evalCPU(int startLoc, int startTime, double* cachedValues)
    {
        //NOT IMPLEMENTED
        throw(-1);
    }

    int FC_P_EI::evalOCL()
    {
        //NOT IMPLEMENTED
        return -1;
    }
    int FC_P_EI::calculateRelevantCompartments()
    {
        // Not used, Do nothing
        return(0);
    }
    int FC_P_EI::calculateRelevantCompartments(int startLoc, int startTime)
    {
        //NOT VALID
        throw(-1);
    }

    int FC_P_EI::sampleCPU()
    {
        sampleDouble(*context, *A0, *p_ei, 1, 0.5); 
        return(0);
    }
    int FC_P_EI::sampleOCL()
    {
        //NOT IMPLEMENTED
        return -1;
    }

    double FC_P_EI::getValue()
    {
        return(*(this -> value));
    }



    FC_P_IR::FC_P_IR(ModelContext *_context,
                     CompartmentalModelMatrix *_R_star,
                     CompartmentalModelMatrix *_I,
                     InitData *_A0,
                     double *_p_ir)
    {

        context = new ModelContext*;
        R_star = new CompartmentalModelMatrix*;
        I = new CompartmentalModelMatrix*;
        A0 = new InitData*;
        p_ir = new double*;
        value = new double;

        *context = _context;
        *R_star = _R_star;
        *I = _I;
        *A0 = _A0;
        *p_ir = _p_ir;
        *value = -1.0;

    }
    FC_P_IR::~FC_P_IR()
    {
        delete R_star;
        delete I;
        delete A0;
        delete p_ir;
        delete value;
        delete context;
    }
    int FC_P_IR::cacheEvalCalculation(double* cachedValues)
    {
        //Not Implemented
        throw(-1);
    }

    int FC_P_IR::evalCPU()
    {
        *value = 0.0;
        int i, j, tmp, compIdx;
        int nLoc = *((*A0) -> numLocations);
        int nTpts = *((*I) -> ncol);
        double term1, term2, term3;
        term1 = 0.0; term2 = 0.0; term3 = 0.0;
        int r_star_sum = 0;
        int i_sum = 0;
        for (j =0; j<(nLoc*nTpts); j++)
        {
            r_star_sum += ((*R_star)-> data)[j];
            i_sum += ((*I)-> data)[j];
        }

        *value = dbeta(**p_ir, 1.5 + r_star_sum, 1.5 - r_star_sum + i_sum); 
        return 0;
    }
    int FC_P_IR::evalCPU(int startLoc, int startTime, double* cachedValues)
    {
        //NOT IMPLEMENTED
        throw(-1);
    }
    int FC_P_IR::evalOCL()
    {
        //NOT IMPLEMENTED
        return -1;
    }
    int FC_P_IR::calculateRelevantCompartments()
    {
        //NOT VALID
        throw(-1);
    }
    int FC_P_IR::calculateRelevantCompartments(int startLoc, int startTime)
    {
        //NOT VALID
        throw(-1);
    }

    int FC_P_IR::sampleCPU()
    {
        //NOT IMPLEMENTED
        return -1;
    }
    int FC_P_IR::sampleOCL()
    {
        //NOT IMPLEMENTED
        return -1;
    }

    double FC_P_IR::getValue()
    {
        return(*(this -> value));
    }
}



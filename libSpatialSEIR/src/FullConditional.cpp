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
        if (Acol != Brow)
        {
            std::cerr << "Invalid Matrix Dimensions" << std::endl;
            throw(-1);
        }
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


    int CompartmentFullConditional::sampleCompartment(ModelContext* context,
                                       InitData* A0,
                                       CompartmentalModelMatrix* starCompartment,
                                       double width,double* cachedValues)
    {
        // Declare required variables
        int i, j, compIdx;
        int nLoc = *(A0 -> numLocations);
        int nTpts = *(starCompartment -> ncol);
        double l,r,y,x,x0;
        
        // Update the relevant CompartmentalModelMatrix instances
        this -> calculateRelevantCompartments();
        
        // Cache Eval Calculation
        this -> cacheEvalCalculation(cachedValues);

        // Set the "value" attribute appropriately
        this -> evalCPU();
   
        int itrs;
        // Main loop: 
        for (j = 0; j < nTpts; j ++)
        { 
            //std::cout << j << "\n";
            compIdx = j*nLoc - 1;
            for (i = 0; i < nLoc; i++)
            {
                compIdx++;
                x = (starCompartment -> data)[compIdx];
                this -> calculateRelevantCompartments(i,j); 
                this -> evalCPU(i,j,cachedValues);
                y = (this->getValue()) - (context -> random -> gamma());
                l = std::max(0.0, (x - (context -> random -> uniform())*width));
                r = l + width;
                itrs = 0;
                do 
                {
                    itrs ++;
                    if (itrs > 1000)
                    {
                        if (!std::isfinite(this -> getValue()))
                        {
                            std::cout << "(Val, y): (" << (this -> getValue()) << ", " << y << ")" << "\n";
                            std::cout << "(x, x0): (" << x << ", "<< x0 << ")" << "\n";
                            std::cout << "l,r: " << l << ", " << r << "\n";
                            throw(-1);
                        }
                        else
                        {
                            y = -INFINITY;
                        }
                    }
                    x0 = ((context -> random -> uniform())*(r-l) + l);
                    (starCompartment -> data)[compIdx] = std::floor(x0);
                    this -> calculateRelevantCompartments(i,j);
                    this -> evalCPU(i,j,cachedValues);
                    if (x0 >= x){r = x0;}
                    else{l = x0;}
                } while (y >= (this -> getValue()));
            }
        }
        return 0;
    }

    int CompartmentFullConditional::sampleCompartmentMemoized(ModelContext* context,
                                           InitData* A0,
                                           CompartmentalModelMatrix* starCompartment,
                                           int width,double* cachedValues)
    {
        // Declare required variables
        double* memo = new double[width + 1];
        int* memoStored = new int[width + 1];

        int i, j, compIdx, memoStart, memoIdx, lastItrMemo;
        int memoIdxBytes = sizeof(int)*(width+1);
        int nLoc = *(A0 -> numLocations);
        int nTpts = *(starCompartment -> ncol);
        double l,r,y,x0;
        int x, x0_floor;
        
        memset(memo, 0, memoIdxBytes);
        memset(memoStored, 0, memoIdxBytes);

        // Update the relevant CompartmentalModelMatrix instances
        this -> calculateRelevantCompartments();
        
        // Cache Eval Calculation
        this -> cacheEvalCalculation(cachedValues);

        // Set the "value" attribute appropriately
        this -> evalCPU();
   
        int itrs;
        // Main loop: 
        for (j = 0; j < nTpts; j ++)
        { 
            //std::cout << j << "\n";
            compIdx = j*nLoc - 1;
            for (i = 0; i < nLoc; i++)
            {
                memset(memoStored, 0, memoIdxBytes);
                lastItrMemo = 0;
                compIdx++;
                x = (starCompartment -> data)[compIdx];
                this -> calculateRelevantCompartments(i,j); 
                this -> evalCPU(i,j,cachedValues);
                y = (this->getValue()) - (context -> random -> gamma());
                l = std::max(0.0, (x - (context -> random -> uniform())*width));
                r = l + width;
                memoStart = l;
                memoIdx = (x-memoStart);
                memo[memoIdx] = (this->getValue());
                memoStored[memoIdx] = 1;
                if (! std::isfinite(y))
                {
                    std::cerr << "Beginning Sampling from location with 0 probability.\n"; 
                    throw(-1);
                }
                itrs = 0;
                do 
                {
                    itrs ++;
                    if (itrs > 1000000)
                    {
                        std::cout << "Sampling Error:\n";
                        std::cout << "Current Value: " << (this -> getValue()) << "\n";
                        std::cout << "Starting Value: " << y << "\n";
                        std::cout << "(x, x0, l, r): (" << x << ", " << x0 << ", " << l  << ", " << r << ")\n";
                        if (std::isfinite(y))
                        {
                            (starCompartment -> data)[compIdx] = x;
                            (this -> setValue(memo[x-memoStart]));
                            return(0);
                        }
                        throw(-1);
                    }
                    x0 = ((context -> random -> uniform())*(r-l) + l);
                    x0_floor = std::floor(x0);
                    memoIdx = (x0_floor - memoStart);
                    if (memoStored[memoIdx])
                    {
                        lastItrMemo = 1;
                        (starCompartment -> data)[compIdx] = x0_floor;
                        (this -> setValue(memo[memoIdx]));  
                        this -> calculateRelevantCompartments(i,j);
                    }
                    else 
                    {
                        lastItrMemo = 0;
                        (starCompartment -> data)[compIdx] = x0_floor;
                        this -> calculateRelevantCompartments(i,j);
                        this -> evalCPU(i,j,cachedValues);
                        if (x0 >= x){r = x0;}
                        else{l = x0;}
                        memoStored[memoIdx] = 1;
                        memo[memoIdx] = (this->getValue());
                    }
                } while (y >= (this -> getValue()));
                if (lastItrMemo)
                {
                    this -> updateEvalCache(i,j,cachedValues);
                }
            }
        }
        return 0;

    }

    int CompartmentFullConditional::sampleCompartmentMetropolis(ModelContext* context,
                                                                InitData* A0,
                                                                CompartmentalModelMatrix* starCompartment,
                                                                double width,double* cachedValues)
    {
        // Declare required variables
        int i, j, compIdx;
        int nLoc = *(A0 -> numLocations);
        int nTpts = *(starCompartment -> ncol);
        int x0, x1;
        double initVal, newVal;
        double initProposal, newProposal;
        double criterion;
        
        // Update the relevant CompartmentalModelMatrix instances
        this -> calculateRelevantCompartments();
        
        // Cache Eval Calculation
        this -> cacheEvalCalculation(cachedValues);

        // Set the "value" attribute appropriately
        //this -> evalCPU();
        this -> evalCPU(0,0,cachedValues);
   
        // Main loop: 
        for (j = 0; j < nTpts; j ++)
        { 
            //std::cout << j << "\n";
            compIdx = j*nLoc - 1;
            for (i = 0; i < nLoc; i++)
            {
                compIdx++;
                this -> calculateRelevantCompartments(i,j); 
                this -> evalCPU(i,j,cachedValues);
                x0 = (starCompartment -> data)[compIdx];
                initVal = (this -> getValue());
  
                // Propose new value, bounded away from zero. 
                x1 = std::floor(std::max(0.0,(context -> random -> normal(x0,width))));
                (starCompartment -> data)[compIdx] = x1;
                this -> calculateRelevantCompartments(i,j);
                this -> evalCPU(i,j,cachedValues);
                newVal = (this->getValue());
                newProposal = (context -> random -> dnorm(x1, x0,width));
                initProposal = (context -> random -> dnorm(x0, x1,width));
                criterion = (newVal - initVal) + (initProposal - newProposal);
                if (std::log((context -> random -> uniform())) < criterion)
                {
                    // Accept new value
                }
                else
                {
                    // Keep Original Value
                    (starCompartment -> data)[compIdx] = x0;
                    this -> calculateRelevantCompartments(i,j);
                    this -> setValue(initVal); 
                }                
            }
        }
        return 0;
    }


    int ParameterFullConditional::sampleDouble(ModelContext* context,
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
            l = x - ((context -> random -> uniform())*width);
            r = l + width;

            if (! std::isfinite(y))
            {
                std::cerr << "Beginning Sampling from location with 0 probability.\n"; 
                throw(-1);
            }

            do
            {
                x0 = ((context -> random -> uniform())*(r-l) + l);
                variable[i] = x0;
                this -> calculateRelevantCompartments();
                this -> evalCPU();
                l = (x0 >= x ? l : x0);
                r = (x0 < x ? r : x0);  
            } while (y >= (this -> getValue()));
        }
        return 0;
    }

    int ParameterFullConditional::sampleDoubleMetropolis(ModelContext* context,
                                                         double* variable, 
                                                         int varLen, 
                                                         double width)
    {
        // Declare required variables
        int i;
        double x0,x1;
        double initVal, newVal, initProposal, newProposal;

        // Update the relevant CompartmentalModelMatrix instances
        this -> calculateRelevantCompartments();

        // Set the "value" attribute appropriately
        this -> evalCPU();
   
        // Main loop: 
        for (i = 0; i < varLen; i++)
        { 
            x0 = variable[i];
            this -> calculateRelevantCompartments(); 
            this -> evalCPU();
            initVal = (this->getValue());

            x1 = (context -> random -> normal(x0, width));
            variable[i] = x1;
            this -> calculateRelevantCompartments();
            this -> evalCPU();
            newVal = (this->getValue());
            initProposal = (context->random->dnorm(x0, x1, width));
            newProposal = (context->random->dnorm(x1, x0, width)); 

            if (std::log((context -> random -> uniform())) < ((newVal - initVal) + (initProposal - newProposal)))
            {
                // Accept the new value.
            }
            else
            {
                // Keep original value
                variable[i] = x0;
                this -> calculateRelevantCompartments();
                this -> setValue(initVal);
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
                         CompartmentalModelMatrix *_E_star,
                         CompartmentalModelMatrix *_R_star,
                         InitData *_A0,
                         CovariateMatrix *_X, 
                         double *_p_se,
                         double *_p_rs,
                         double *_beta,
                         double *_rho,
                         double _steadyStateConstraintPrecision,
                         double _sliceWidth)
    {
       context = new ModelContext*;
       S_star = new CompartmentalModelMatrix*;
       S = new CompartmentalModelMatrix*;
       R = new CompartmentalModelMatrix*;
       E_star = new CompartmentalModelMatrix*;
       R_star = new CompartmentalModelMatrix*;
       A0 = new InitData*;
       X = new CovariateMatrix*;
       p_se = new double*;
       p_rs = new double*;
       beta = new double*;
       rho = new double*;
       value = new double;
       steadyStateConstraintPrecision = new double;
       sliceWidth = new double;
       *context = _context;
       *S_star = _S_star;
       *S = _S;
       *R = _R;
       *E_star = _E_star;
       *R_star = _R_star;
       *A0 = _A0;
       *X = _X;
       *p_se = _p_se;
       *p_rs = _p_rs;
       *beta = _beta;
       *rho = _rho;
       *steadyStateConstraintPrecision = _steadyStateConstraintPrecision;
       *value = -1.0;
       *sliceWidth = _sliceWidth;
    }    
    FC_S_Star::~FC_S_Star()
    {
        delete S_star;
        delete S;
        delete R;
        delete E_star;
        delete R_star;
        delete A0;
        delete X;
        delete p_se;
        delete p_rs;
        delete beta;
        delete rho;
        delete value;
        delete steadyStateConstraintPrecision;
        delete sliceWidth;
        delete context;
    }

    // Cache the components of the FC_S_Star calculation in cachedValues
    int FC_S_Star::cacheEvalCalculation(double* cachedValues)
    {
        // Cache eval calulation could probably get away with 
        // returning upon encountering bad data, but the performance
        // benefit would be small given the relaitive infrequency of 
        // calls.
        int nLoc = *((*S)->nrow);
        int nTpts = *((*S)->ncol);

        int i,j,compIdx;
        int Sstar_val,Estar_val,S_val,R_val;
        double p_se_val, p_rs_val;
        int err = 0;

        for (j = 0; j < nTpts; j++)
        {
            compIdx = j*nLoc - 1;
            p_rs_val = (*p_rs)[j];
            for (i = 0;i<nLoc;i++);
            {
                compIdx ++; 
                Sstar_val = ((*S_star)->data)[compIdx]; 
                Estar_val = ((*E_star)->data)[compIdx];
                S_val = ((*S)->data)[compIdx];
                R_val = ((*R)->data)[compIdx];
                p_se_val = (*p_se)[compIdx];

                // Check Constraints
                if (Sstar_val < 0 || 
                    Sstar_val > R_val ||
                    Estar_val > S_val)
                {
                    cachedValues[compIdx] = -INFINITY;
                    err += 1;
                }
                else
                { 
                    cachedValues[compIdx] = (std::log(p_rs_val)*Sstar_val + 
                                   std::log(1-p_rs_val)*(R_val - Sstar_val) +
                                   std::log(1-p_se_val)*(S_val) + 
                                   ((*context) -> random -> choose(R_val, Sstar_val)) + 
                                   ((*context)->random->choosePartial(S_val, Estar_val)));
                    if (!std::isfinite(cachedValues[compIdx]))
                    {
                        cachedValues[compIdx] = -INFINITY;
                        err += 1;
                    }
                }
            }
        }
        return(err);
    }

    int FC_S_Star::updateEvalCache(int startLoc, int startTime, double* cachedValues)
    {
        // If update evalCache encounters invalid data, we're ALWAYS
        // going to end up back here at the same startLoc startTime with a new 
        // proposal. Therefore, return immediately upon encountering bad 
        // values. 
        int i,compIdx,Sstar_val,Estar_val,S_val,R_val;
        double p_se_val, p_rs_val;
        int nLoc = *((*S)->nrow);  
        int nTpts = *((*S)->ncol);
        
        compIdx = startLoc + startTime*nLoc;
        for (i = startTime; i < nTpts; i++)
        {
                Sstar_val = ((*S_star)->data)[compIdx]; 
                Estar_val = ((*E_star)->data)[compIdx];
                S_val = ((*S)->data)[compIdx];
                R_val = ((*R)->data)[compIdx];
                p_se_val = (*p_se)[compIdx];
                p_rs_val = (*p_rs)[i];

                if (Sstar_val < 0 || 
                    Sstar_val > R_val ||
                    Estar_val > S_val)
                {
                    cachedValues[compIdx] = -INFINITY;
                    return(-1);
                }
                else
                { 
                    cachedValues[compIdx] = (std::log(p_rs_val)*Sstar_val + 
                                   std::log(1-p_rs_val)*(R_val - Sstar_val) +
                                   std::log(1-p_se_val)*(S_val) + 
                                   ((*context) -> random -> choose(R_val, Sstar_val)) + 
                                   ((*context)->random->choosePartial(S_val, Estar_val)));
                    if (!std::isfinite(cachedValues[compIdx]))
                    {
                        cachedValues[compIdx] = -INFINITY;
                        return(-1);
                    }
                }
                compIdx += nLoc; 
        }
        return(0);
    }

    // Evaluate the S_star FC at the current values provided by the context.
    int FC_S_Star::evalCPU()
    {
        int nLoc = *((*S)->nrow);
        int nTpts = *((*S)->ncol);

        int i,j,compIdx;
        long double output = 0.0;
        int Sstar_val,Estar_val,S_val,R_val;
        double p_se_val, p_rs_val;
        *value = 0.0;
        long unsigned int R_star_sum;
        long unsigned int S_star_sum;
        int64_t aDiff; 

        for (j = 0; j < nTpts; j++)
        {
            compIdx = j*nLoc - 1;
            p_rs_val = (*p_rs)[j];
            for (i = 0;i<nLoc;i++);
            {
                compIdx ++; 
                Sstar_val = ((*S_star)->data)[compIdx]; 
                Estar_val = ((*E_star)->data)[compIdx];
                S_val = ((*S)->data)[compIdx];
                R_val = ((*R)->data)[compIdx];
                p_se_val = (*p_se)[compIdx];


                // Check Constraints
                if (Sstar_val < 0 || 
                    Sstar_val > R_val ||
                    Estar_val > S_val)
                {
                    *value = -INFINITY;
                    return(-1);
                }
                
                output += (std::log(p_rs_val)*Sstar_val + 
                           std::log(1-p_rs_val)*(R_val - Sstar_val) +
                           std::log(1-p_se_val)*(S_val) + 
                           ((*context) -> random -> choose(R_val, Sstar_val)) + 
                           ((*context)->random->choosePartial(S_val, Estar_val)));
            }

        }
        *value = output;

        for (i=0; i< nLoc; i++)
        {
            S_star_sum = (*S_star)->marginSum(1,i);
            R_star_sum = (*R_star)->marginSum(1,i);
            aDiff = (S_star_sum > R_star_sum ? S_star_sum - R_star_sum : R_star_sum - S_star_sum);
            *value -= (aDiff*aDiff)*(*steadyStateConstraintPrecision);
        }

        if (!std::isfinite(*value))
        {
            *value = -INFINITY;
            return(-1); 
        }
        return(0);

    }
    int FC_S_Star::evalCPU(int startLoc, int startTime, double* cachedValues)
    {
        int valid = 1;
        int nvals = (*((*S)->nrow))*(*((*S)->ncol));
        int i;
        long double output = 0.0;
        long unsigned int S_star_sum;
        long unsigned int R_star_sum;
        int64_t aDiff;

        valid = updateEvalCache(startLoc, startTime, cachedValues);
        if (valid < 0)
        {
            *value = -INFINITY;
            return(-1);
        }

        for (i = 0; i < nvals; i++)
        {
            output += cachedValues[i];
        }
        *value = output;

        S_star_sum = (*S_star)->marginSum(1,startLoc);
        R_star_sum = (*R_star)->marginSum(1,startLoc);
        aDiff = (S_star_sum > R_star_sum ? S_star_sum - R_star_sum : R_star_sum - S_star_sum);
        *value -= (aDiff*aDiff)*(*steadyStateConstraintPrecision);

        if (!std::isfinite(*value))
        {
            *value = -INFINITY;
            return(-1);
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
        (*context) -> calculateR_givenS_CPU();
        return(0);

    }
    int FC_S_Star::calculateRelevantCompartments(int startLoc, int startTime)
    {
        (*context) -> calculateS_CPU(startLoc, startTime);
        (*context) -> calculateR_givenS_CPU(startLoc, startTime);
        return(0);
    }

    int FC_S_Star::sampleCPU()
    {
        this -> sampleCompartmentMetropolis(*context,*A0,
                                  *S_star,*sliceWidth,(*context) -> compartmentCache);
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
    void FC_S_Star::setValue(double val)
    {
        *(this -> value) = val;
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
                         CompartmentalModelMatrix *_I_star,
                         CovariateMatrix *_X,
                         InitData *_A0,
                         double *_p_se,
                         double *_p_ei,
                         double *_rho,
                         double *_beta,
                         double _steadyStateConstraintPrecision,
                         double _sliceWidth) 
    {

        context = new ModelContext*;
        E_star = new CompartmentalModelMatrix*;
        E = new CompartmentalModelMatrix*;
        S = new CompartmentalModelMatrix*;
        I_star = new CompartmentalModelMatrix*;
        X = new CovariateMatrix*;
        A0 = new InitData*;
        p_se = new double*;
        p_ei = new double*;
        rho = new double*;
        beta = new double*;
        sliceWidth = new double;
        steadyStateConstraintPrecision = new double;
        value = new double;
        
        *context = _context;
        *E_star = _E_star;
        *E = _E;
        *S = _S;
        *I_star = _I_star;
        *X = _X;
        *A0 = _A0;
        *p_se = _p_se;
        *p_ei = _p_ei;
        *rho = _rho;
        *beta = _beta;
        *sliceWidth = _sliceWidth;
        *steadyStateConstraintPrecision = _steadyStateConstraintPrecision;
        *value = -1.0;
    }

    FC_E_Star::~FC_E_Star()
    {
        delete E_star;
        delete E;
        delete S;
        delete I_star;
        delete X;
        delete A0;
        delete p_se;
        delete p_ei;
        delete rho;
        delete beta;
        delete value;
        delete sliceWidth;
        delete steadyStateConstraintPrecision; 
        delete context;
    }

    int FC_E_Star::cacheEvalCalculation(double* cachedValues)
    {
        int i, j, compIdx;
        int nLoc = *((*A0) -> numLocations);
        int nTpts = *((*S) -> ncol);
        double ln_1m_p_ei = std::log(1-(**p_ei));     
        double p_se_val;
        int S_val, E_val, Estar_val, Istar_val;
        for (j = 0; j < nTpts; j++)     
        {
            compIdx = j*nLoc - 1;
            for (i = 0; i < nLoc; i++)    
            {
                compIdx++;
                Estar_val = ((*E_star) -> data)[compIdx];
                S_val = ((*S) -> data)[compIdx];
                E_val = ((*E) -> data)[compIdx];
                Istar_val = ((*I_star) -> data)[compIdx];
                p_se_val = (*p_se)[compIdx];

                // Check conditions
                if (Estar_val < 0 || Estar_val > S_val || 
                        Istar_val > E_val)
                {
                    cachedValues[compIdx] = (-INFINITY);
                }
                else
                {
                    cachedValues[compIdx] = (std::log(p_se_val)*Estar_val +
                                             std::log(1-p_se_val)*(S_val - Estar_val) +
                                             ln_1m_p_ei*E_val + 
                                             ((*context) -> random -> choosePartial(E_val, Istar_val)) + 
                                             ((*context)->random->choose(S_val, Estar_val)));
                }
            }
        } 
        return(0);
    }

    int FC_E_Star::updateEvalCache(int startLoc, int startTime, double* cachedValues)
    {
        int j, compIdx;
        int nLoc = *((*A0) -> numLocations);
        int nTpts = *((*S) -> ncol);
        double ln_1m_p_ei = std::log(1-(**p_ei));    
        double p_se_val;
        int S_val, E_val, Estar_val, Istar_val;

        compIdx = startLoc + startTime*nLoc;
        for (j = startTime; j < nTpts; j++)
        {
            Estar_val = ((*E_star) -> data)[compIdx];
            S_val = ((*S) -> data)[compIdx];
            E_val = ((*E) -> data)[compIdx];
            Istar_val = ((*I_star) -> data)[compIdx];
            p_se_val = (*p_se)[compIdx];

            if (Estar_val < 0 || Estar_val > S_val || 
                    Istar_val > E_val)

            {
                cachedValues[compIdx] = -INFINITY;
                return(-1);
            }
            else
            {
                cachedValues[compIdx] = (std::log(p_se_val)*Estar_val +
                                         std::log(1-p_se_val)*(S_val - Estar_val) +
                                         ln_1m_p_ei*E_val + 
                                         ((*context) -> random -> choosePartial(E_val, Istar_val)) + 
                                         ((*context)->random->choose(S_val, Estar_val)));
                if (!std::isfinite(cachedValues[compIdx]))
                {
                    cachedValues[compIdx] = -INFINITY;
                    return(-1);
                }
            }
            compIdx += nLoc;
        }
       return 0;
    }

    int FC_E_Star::evalCPU()
    {

        int i, j, compIdx;
        int nLoc = *((*A0) -> numLocations);
        int nTpts = *((*S) -> ncol);
        double ln_1m_p_ei = std::log(1-(**p_ei));     
        double p_se_val;
        int S_val, E_val, Estar_val, Istar_val;
        long double output = 0.0;
        long unsigned int E_star_sum;
        long unsigned int I_star_sum;
        int64_t aDiff; 

        for (j = 0; j < nTpts; j++)     
        {
            compIdx = j*nLoc - 1;
            for (i = 0; i < nLoc; i++)    
            {
                compIdx++;
                Estar_val = ((*E_star) -> data)[compIdx];
                S_val = ((*S) -> data)[compIdx];
                E_val = ((*E) -> data)[compIdx];
                Istar_val = ((*I_star) -> data)[compIdx];
                p_se_val = (*p_se)[compIdx];

                // Check conditions
                if (Estar_val < 0 || Estar_val > S_val || 
                        Istar_val > E_val)
                {
                    *value = (-INFINITY);
                    return(-1);
                }
                else
                {
                    output += (std::log(p_se_val)*Estar_val +
                              std::log(1-p_se_val)*(S_val - Estar_val) +
                              ln_1m_p_ei*E_val + 
                              ((*context) -> random -> choosePartial(E_val, Istar_val)) + 
                              ((*context)->random->choose(S_val, Estar_val)));
                }
            }
        } 
        

        for (i=0; i< nLoc; i++)
        {
            E_star_sum = (*E_star)->marginSum(1,i);
            I_star_sum = (*I_star)->marginSum(1,i);
            aDiff = (E_star_sum > I_star_sum ? E_star_sum - I_star_sum : I_star_sum - E_star_sum);
            output -= (aDiff*aDiff)*(*steadyStateConstraintPrecision);
        }
        if (!std::isfinite(output))
        {
            *value = -INFINITY;
            return(-1);
        }
        *value = output;
        return(0);



    }
    int FC_E_Star::evalCPU(int startLoc, int startTime, double* cachedValues)
    {
        int err;
        long unsigned int E_star_sum;
        long unsigned int I_star_sum;
        int64_t aDiff; 

        err = updateEvalCache(startLoc, startTime, cachedValues);
        if (err < 0)
        {
            *value = -INFINITY; 
            return(-1);
        }
        int maxVal = (*((*E)-> ncol))*(*((*E)-> nrow));
        int i;
        long double output = 0.0;
        for (i = 0; i < maxVal; i++)
        {
            output += cachedValues[i];
        }
        if (!std::isfinite(output))
        {
            *value = -INFINITY;
            return(-1);
        }
        E_star_sum = (*E_star)->marginSum(1,startLoc);
        I_star_sum = (*I_star)->marginSum(1,startLoc);
        aDiff = (E_star_sum > I_star_sum ? E_star_sum - I_star_sum : I_star_sum - E_star_sum);
        output -= (aDiff*aDiff)*(*steadyStateConstraintPrecision);
        *value = output;
        return(0);
    }

    int FC_E_Star::evalOCL()
    {
        //NOT IMPLEMENTED
        return -1;
    }
    int FC_E_Star::calculateRelevantCompartments()
    {
        (*context) -> calculateE_CPU();
        (*context) -> calculateS_givenE_CPU();
        return(0);
    }
    int FC_E_Star::calculateRelevantCompartments(int startLoc, int startTime)
    {
        (*context) -> calculateE_CPU(startLoc, startTime);
        (*context) -> calculateS_givenE_CPU(startLoc, startTime);
        return(0);
    }

    int FC_E_Star::sampleCPU()
    {
        this -> sampleCompartmentMetropolis(*context,*A0,
                                  *E_star,*sliceWidth,(*context) -> compartmentCache);
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
    void FC_E_Star::setValue(double val)
    {
        *(this -> value) = val;
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
                         CompartmentalModelMatrix *_S_star,
                         CompartmentalModelMatrix *_E_star,
                         CompartmentalModelMatrix *_I_star,
                         CompartmentalModelMatrix *_S,
                         InitData *_A0,
                         double *_p_rs,
                         double *_p_ir,
                         double *_p_se,
                         double _steadyStateConstraintPrecision,
                         double _sliceWidth)
    {

        context = new ModelContext*;
        R_star = new CompartmentalModelMatrix*;
        R = new CompartmentalModelMatrix*;
        I = new CompartmentalModelMatrix*;
        S_star = new CompartmentalModelMatrix*;
        E_star = new CompartmentalModelMatrix*;
        I_star = new CompartmentalModelMatrix*;
        S = new CompartmentalModelMatrix*;
        A0 = new InitData*;
        p_rs = new double*;
        p_ir = new double*;
        p_se = new double*;
        sliceWidth = new double;
        steadyStateConstraintPrecision = new double;
        value = new double;

        *context = _context;
        *R_star = _R_star;
        *R = _R;
        *I = _I;
        *S_star = _S_star;
        *E_star = _E_star;
        *I_star = _I_star;
        *S = _S;
        *A0 = _A0;
        *p_rs = _p_rs;
        *p_ir = _p_ir;
        *p_se = _p_se;
        *sliceWidth = _sliceWidth;
        *steadyStateConstraintPrecision = _steadyStateConstraintPrecision;
        *value = -1.0;
    }
    FC_R_Star::~FC_R_Star()
    {
        delete R_star;
        delete R;
        delete I;
        delete S_star;
        delete E_star;
        delete I_star;
        delete S;
        delete A0;
        delete p_rs;
        delete p_ir;
        delete p_se;
        delete value;
        delete sliceWidth;
        delete steadyStateConstraintPrecision;
        delete context;
    }

    int FC_R_Star::cacheEvalCalculation(double* cachedValues)
    {
        int i, j, compIdx;
        int nLoc = *((*A0) -> numLocations);
        int nTpts = *((*R) -> ncol);
        double p_se_val;
        double p_rs_val;
        double ln_p_ir = std::log(**p_ir);
        double ln_1m_p_ir = std::log(1-(**p_ir));
        int Rstar_val, Sstar_val, Estar_val, R_val, I_val, S_val;

    
        for (j = 0; j < nTpts; j++)     
        {
            compIdx = j*nLoc - 1;
            p_rs_val = (*p_rs)[j];          
            for (i = 0; i < nLoc; i++)    
            {
                compIdx++;
                Rstar_val = ((*R_star) -> data)[compIdx];
                Sstar_val = ((*S_star)->data)[compIdx];
                Estar_val = ((*E_star) -> data)[compIdx];
                R_val = ((*R) ->data)[compIdx];
                I_val = ((*I) ->data)[compIdx];
                S_val = ((*S)->data)[compIdx];
                p_se_val = (*p_se)[compIdx];


                if (Rstar_val < 0 || Rstar_val > I_val || 
                        Sstar_val > R_val ||
                         p_se_val >= 1 || p_se_val <= 0)

                {
                    cachedValues[compIdx] = -INFINITY;
                }
                else
                {
                    cachedValues[compIdx] = (ln_p_ir*Rstar_val + 
                                    ln_1m_p_ir*(I_val -Rstar_val) +
                                    std::log(1-p_rs_val)*(R_val) +
                                    ((*context) -> random -> choosePartial(R_val, Sstar_val)) +
                                    ((*context) -> random -> choose(I_val, Rstar_val)) +
                                    // p_se components
                                    std::log(p_se_val)*Estar_val +
                                    std::log(1-p_se_val)*(S_val - Estar_val));

                }
            }
        } 
        return(0);
    }

    int FC_R_Star::updateEvalCache(int startLoc, int startTime, double* cachedValues)
    {
        int j, compIdx;
        int nLoc = *((*A0) -> numLocations);
        int nTpts = *((*R) -> ncol);
        
        double p_se_val;
        double p_rs_val;
        double ln_p_ir = std::log(**p_ir);
        double ln_1m_p_ir = std::log(1-(**p_ir));
        int Rstar_val, Sstar_val, Estar_val, R_val, I_val, S_val;   

        compIdx = startLoc + startTime*nLoc;
        for (j = startTime; j < nTpts; j++)
        {

            Rstar_val = ((*R_star) -> data)[compIdx];
            Sstar_val = ((*S_star)->data)[compIdx];
            Estar_val = ((*E_star) -> data)[compIdx];
            R_val = ((*R) ->data)[compIdx];
            I_val = ((*I) ->data)[compIdx];
            S_val = ((*S)->data)[compIdx];
            p_se_val = (*p_se)[compIdx];
            p_rs_val = (*p_rs)[j];


            if (Rstar_val < 0 || Rstar_val > I_val || 
                    Sstar_val > R_val || 
                    p_se_val >= 1 || p_se_val <= 0)
            {
                cachedValues[compIdx] = -INFINITY;
                return(-1);
            }
            else
            {
                cachedValues[compIdx] = (ln_p_ir*Rstar_val +
                                ln_1m_p_ir*(I_val - Rstar_val) +
                                std::log(1-p_rs_val)*(R_val) +
                                ((*context) -> random -> choosePartial(R_val, Sstar_val)) +
                                ((*context) -> random -> choose(I_val, Rstar_val)) +
                                std::log(p_se_val)*Estar_val +
                                std::log(1-p_se_val)*(S_val - Estar_val));
                if (!std::isfinite(cachedValues[compIdx]))
                {
                    cachedValues[compIdx] = -INFINITY;
                    return(-1);
                }
            }
            compIdx += nLoc;
        } 
        return 0;
    }


    int FC_R_Star::evalCPU()
    {
        int i, j, compIdx;
        int nLoc = *((*A0) -> numLocations);
        int nTpts = *((*R) -> ncol);
        double p_se_val;
        double p_rs_val;
        double ln_p_ir = std::log(**p_ir);
        double ln_1m_p_ir = std::log(1-(**p_ir));
        int Rstar_val, Sstar_val, Estar_val, R_val, I_val,S_val;
        long double output = 0.0;
        long unsigned int I_star_sum;
        long unsigned int R_star_sum;
        int64_t aDiff; 

    
        for (j = 0; j < nTpts; j++)     
        {
            compIdx = j*nLoc - 1;
            p_rs_val = (*p_rs)[j];          
            for (i = 0; i < nLoc; i++)    
            {
                compIdx++;
                Rstar_val = ((*R_star) -> data)[compIdx];
                Sstar_val = ((*S_star)->data)[compIdx];
                Estar_val = ((*E_star) -> data)[compIdx];
                R_val = ((*R) ->data)[compIdx];
                I_val = ((*I) ->data)[compIdx];
                S_val = ((*S)->data)[compIdx];
                p_se_val = (*p_se)[compIdx];

                if (Rstar_val < 0 || Rstar_val > I_val || 
                        Sstar_val > R_val ||
                         p_se_val >= 1 || p_se_val <= 0)

                {
                    *value = -INFINITY;
                    return(-1);
                }
                else
                {
                    output += (ln_p_ir*Rstar_val + 
                                ln_1m_p_ir*(I_val -Rstar_val) +
                                std::log(1-p_rs_val)*(R_val) +
                                ((*context) -> random -> choosePartial(R_val, Sstar_val)) +
                                ((*context) -> random -> choose(I_val, Rstar_val)) +

                                // p_se components
                                std::log(p_se_val)*Estar_val +
                                std::log(1-p_se_val)*(S_val - Estar_val));
                }
            }
        } 

        for (i=0; i< nLoc; i++)
        {
            I_star_sum = (*I_star)->marginSum(1,i);
            R_star_sum = (*R_star)->marginSum(1,i);
            aDiff = (I_star_sum > R_star_sum ? I_star_sum - R_star_sum : R_star_sum - I_star_sum);
            output -= (aDiff*aDiff)*(*steadyStateConstraintPrecision);
        }

        if (!std::isfinite(*value))
        {
            *value = -INFINITY;
            return(-1);
        }
        *value = output;
        return(0);
    }

    int FC_R_Star::evalCPU(int startLoc, int startTime, double* cachedValues)
    {
        int err;
        int maxVal = (*((*R) -> ncol))*(*((*R)->nrow));
        int i;
        long double output = 0.0;
        long unsigned int I_star_sum;
        long unsigned int R_star_sum;
        int64_t aDiff; 

        err = updateEvalCache(startLoc, startTime, cachedValues);
        if (err < 0)
        {
            *value = -INFINITY;
            return(-1);
        }
        for (i = 0; i < maxVal; i++)
        {
            output += cachedValues[i];
        }


        I_star_sum = (*I_star)->marginSum(1,startLoc);
        R_star_sum = (*R_star)->marginSum(1,startLoc);
        aDiff = (I_star_sum > R_star_sum ? I_star_sum - R_star_sum : R_star_sum - I_star_sum);
        output -= (aDiff*aDiff)*(*steadyStateConstraintPrecision);

        if (!std::isfinite(output))
        {
            *value = -INFINITY;
            return(-1);
        }
        *value = output;
        return(0);
    }


    int FC_R_Star::evalOCL()
    {
        //NOT IMPLEMENTED
        return -1;
    }
    int FC_R_Star::calculateRelevantCompartments()
    {
        (*context) -> calculateR_CPU();
        (*context) -> calculateI_givenR_CPU();
        ((*context) -> calculateP_SE_CPU());
        return(0);
    }
    int FC_R_Star::calculateRelevantCompartments(int startLoc, int startTime)
    {
        
        (*context) -> calculateR_CPU(startLoc, startTime);
        (*context) -> calculateI_givenR_CPU(startLoc, startTime);
        ((*context) -> calculateP_SE_CPU(startLoc, startTime));
        return(0);
    }

    int FC_R_Star::sampleCPU()
    {
        this -> sampleCompartmentMetropolis(*context,*A0,
                                  *R_star,*sliceWidth,(*context) -> compartmentCache);
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
    void FC_R_Star::setValue(double val)
    {
        *(this -> value) = val;
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
                     double *_rho,
                     double _sliceWidth,
                     double _priorPrecision)
    {

        context = new ModelContext*;
        E_star = new CompartmentalModelMatrix*;
        S = new CompartmentalModelMatrix*;
        A0 = new InitData*;
        X = new CovariateMatrix*;
        p_se = new double*;
        beta = new double*;
        rho = new double*;
        sliceWidth = new double;
        priorPrecision = new double;
        value = new double;

        *context = _context;
        *E_star = _E_star;
        *S = _S;
        *A0 = _A0;
        *X = _X;
        *p_se = _p_se;
        *beta = _beta;
        *rho = _rho;
        *sliceWidth = _sliceWidth;
        *priorPrecision = _priorPrecision;
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
        delete sliceWidth;
        delete priorPrecision;
        delete context;
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
            compIdx = j*nLoc - 1;
            for (i = 0; i < nLoc; i++)    
            {
                compIdx++;
                tmp = ((*E_star) -> data)[compIdx];
                term1 += std::log((*p_se)[compIdx])*tmp; 
                term2 += std::log(1-(*p_se)[compIdx])*(((*S) -> data)[compIdx] - tmp);
            }
        } 

        for (i = 0; i < (*((*X) -> ncol_x) + *((*X) -> ncol_z)); i++)
        {
            term3 -= pow((*beta)[i],2)*(*priorPrecision); // Generalize to allow different prior precisions. 
        }

        *value = term1 + term2 + term3;
        // Catch invalid values, nans etc. 
        if (!std::isfinite(*value))
        {
            *value = -INFINITY;
        }
        return(0);
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

    int FC_Beta::sampleCPU()
    {
        sampleDoubleMetropolis(*context, *beta, (*((*X) -> ncol_x) + *((*X) -> ncol_z)), *sliceWidth); 
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
    void FC_Beta::setValue(double val)
    {
        *(this -> value) = val;
    }


    /*
     *
     * Implement the full conditional for the R->S transition 
     * probabilities. 
     *
     */
    FC_Beta_P_RS::FC_Beta_P_RS(ModelContext *_context,
                     CompartmentalModelMatrix *_S_star, 
                     CompartmentalModelMatrix *_R,
                     CovariateMatrix* _X,
                     InitData *_A0,
                     double *_p_rs,
                     double *_beta_p_rs,
                     double _tausq,
                     double _sliceWidth)
    {

        context = new ModelContext*;
        S_star = new CompartmentalModelMatrix*;
        R = new CompartmentalModelMatrix*;
        X = new CovariateMatrix*;
        A0 = new InitData*;
        p_rs = new double*;
        beta_p_rs = new double*;
        tausq = new double;
        sliceWidth = new double;
        value = new double;

        *context = _context;
        *S_star = _S_star;
        *X = _X;
        *R = _R;
        *A0 = _A0;
        *p_rs = _p_rs;
        *beta_p_rs = _beta_p_rs;
        *tausq = _tausq;
        *sliceWidth = _sliceWidth;
        *value = -1.0;
    }
    FC_Beta_P_RS::~FC_Beta_P_RS()
    {
        delete S_star;
        delete R;
        delete X;
        delete beta_p_rs;
        delete tausq;
        delete A0;
        delete p_rs;
        delete value;
        delete sliceWidth;
        delete context;
    }

    int FC_Beta_P_RS::evalCPU()
    {
        int j;
        long double a,b;
        int nbeta = *((*X) -> ncol_x);
        int nTpts = *((*R) -> ncol);
        double tmp;
        long double term1 = 0.0;
        double term2 = 0.0;

        for (j = 0; j < nTpts; j++)
        {
            tmp = (*p_rs)[j];
            if (tmp <= 0 || tmp >= 1)
            {
                *value = -INFINITY;
                return(-1);
            }
            a = ((*S_star)-> marginSum(2,j));
            b = ((*R) -> marginSum(2,j)); 
            term1 += std::log(tmp)*a; 
            term1 += std::log(1 - tmp)*(b - a);
        }
        for (j = 0; j < nbeta; j++)
        {
            term2 -= ((*tausq)/2)*pow((*beta_p_rs)[j],2);
        }
        *value = term1 + term2;
        if (!std::isfinite(*value))
        {
            *value = -INFINITY;
        }
        return(0);
    }

    int FC_Beta_P_RS::evalOCL()
    {
        //NOT IMPLEMENTED
        return -1;
    }
    int FC_Beta_P_RS::calculateRelevantCompartments()
    {
         ((*context) -> calculateP_RS_CPU());      
         return(0);
    }

    int FC_Beta_P_RS::sampleCPU()
    {
        int nbeta = *((*X) -> ncol_x);
        sampleDoubleMetropolis(*context, *beta_p_rs, nbeta, *sliceWidth); 
        return(0);
    }

    int FC_Beta_P_RS::sampleOCL()
    {
        //NOT IMPLEMENTED
        return -1;
    }

    double FC_Beta_P_RS::getValue()
    {
        return(*(this -> value));
    }
    void FC_Beta_P_RS::setValue(double val)
    {
        *(this -> value) = val;
    }

    FC_Rho::FC_Rho(ModelContext *_context,
                   CompartmentalModelMatrix *_E_star,  
                   CompartmentalModelMatrix *_S,
                   InitData *_A0,
                   CovariateMatrix *_X,
                   double *_p_se,
                   double *_beta,
                   double *_rho,
                   double _sliceWidth)
    {
        context = new ModelContext*;
        E_star = new CompartmentalModelMatrix*;
        S = new CompartmentalModelMatrix*;
        A0 = new InitData*;
        X = new CovariateMatrix*;
        p_se = new double*;
        beta = new double*;
        rho = new double*;
        sliceWidth = new double;
        value = new double;

        *context = _context;
        *E_star = _E_star;
        *S = _S;
        *A0 = _A0;
        *X = _X;
        *p_se = _p_se;
        *beta = _beta;
        *rho = _rho;
        *sliceWidth = _sliceWidth;
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
        delete sliceWidth;
        delete context;
    }

    int FC_Rho::evalCPU()
    {
        *value = 0.0;
        int i, j, Es, compIdx;
        double pse;
        int nLoc = *((*A0) -> numLocations);
        int nTpts = *((*S) -> ncol);
        double term1, term2, term3;
        term1 = 0.0; term2 = 0.0; term3 = 0.0;
        for (j = 0; j < nTpts; j++)     
        {
            compIdx = j*nLoc - 1;
            for (i = 0; i < nLoc; i++)    
            {
                compIdx++;
                Es = ((*E_star) -> data)[compIdx];
                pse = (*p_se)[compIdx];
                term1 += std::log(pse)*Es; 
                term2 += std::log(1-pse)*(((*S) -> data)[compIdx] - Es);
            }
        } 
        term3 += (**rho > 0 && **rho < 1 ? 0 : -INFINITY); // Generalize to allow informative priors. 
                                                        // Prior specification in this area needs work. 
        *value = term1 + term2 + term3;
        // Catch invalid values, nans etc. 
        if (!std::isfinite(*value))
        {
            *value = -INFINITY;
        }

        return(0);
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

    int FC_Rho::sampleCPU()
    {
        sampleDoubleMetropolis(*context, *rho, 1, *sliceWidth); 
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
    void FC_Rho::setValue(double val)
    {
        *(this -> value) = val;
    }

    FC_Gamma::FC_Gamma(ModelContext *_context,
                   CompartmentalModelMatrix *_E_star,  
                   CompartmentalModelMatrix *_S,
                   InitData *_A0,
                   CovariateMatrix *_X,
                   double *_p_se,
                   double *_beta,
                   double *_gamma,
                   double *_priorAlpha,
                   double *_priorBeta,
                   double _sliceWidth)
    {
        context = new ModelContext*;
        E_star = new CompartmentalModelMatrix*;
        S = new CompartmentalModelMatrix*;
        A0 = new InitData*;
        X = new CovariateMatrix*;
        p_se = new double*;
        beta = new double*;
        gamma = new double*;
        priorAlpha = new double;
        priorBeta = new double;
        sliceWidth = new double;
        value = new double;

        *context = _context;
        *E_star = _E_star;
        *S = _S;
        *A0 = _A0;
        *X = _X;
        *p_se = _p_se;
        *beta = _beta;
        *gamma = _gamma;
        *priorAlpha = *_priorAlpha;
        *priorBeta = *_priorBeta;
        *sliceWidth = _sliceWidth;
        *value = -1.0;
    }
    FC_Gamma::~FC_Gamma()
    {
        delete E_star;
        delete S;
        delete A0;
        delete X;
        delete p_se;
        delete beta;
        delete gamma;
        delete priorAlpha;
        delete priorBeta;
        delete sliceWidth;
        delete value;
        delete context;
    }

    int FC_Gamma::evalCPU()
    {
        *value = 0.0;
        int i, j, Es, compIdx;
        double pse;
        int nLoc = *((*A0) -> numLocations);
        int nTpts = *((*S) -> ncol);
        double term1, term2, term3;
        term1 = 0.0; term2 = 0.0; term3 = 0.0;
        for (j = 0; j < nTpts; j++)     
        {
            compIdx = j*nLoc - 1;
            for (i = 0; i < nLoc; i++)    
            {
                compIdx++;
                Es = ((*E_star) -> data)[compIdx];
                pse = (*p_se)[compIdx];
                term1 += std::log(pse)*Es; 
                term2 += std::log(1-pse)*(((*S) -> data)[compIdx] - Es);
            }
            term3 += ((*priorAlpha-1)*std::log((*gamma)[j]) - ((*gamma)[j])/(*priorBeta)); 
        } 
        *value = term1 + term2 + term3;
        // Catch invalid values, nans etc. 
        if (!std::isfinite(*value))
        {
            *value = -INFINITY;
        }

        return(0);
    }

    int FC_Gamma::evalOCL()
    {
        //NOT IMPLEMENTED
        return -1;
    }
    int FC_Gamma::calculateRelevantCompartments()
    {
       (*context) -> calculateP_SE_CPU();
       return(0); 
    }

    int FC_Gamma::sampleCPU()
    {
        sampleDoubleMetropolis(*context, *gamma, *((*A0) -> numLocations), *sliceWidth); 
        return(0);
    }
    int FC_Gamma::sampleOCL()
    {
        //NOT IMPLEMENTED
        return -1;
    }

    double FC_Gamma::getValue()
    {
        return(*(this -> value));
    }
    void FC_Gamma::setValue(double val)
    {
        *(this -> value) = val;
    }

    FC_P_EI::FC_P_EI(ModelContext *_context,
                     CompartmentalModelMatrix *_I_star,
                     CompartmentalModelMatrix *_E,
                     InitData *_A0,
                     double *_p_ei,
                     double _priorAlpha,
                     double _priorBeta)
    {

        context = new ModelContext*;
        I_star = new CompartmentalModelMatrix*;
        E = new CompartmentalModelMatrix*;
        A0 = new InitData*;
        p_ei = new double*;
        priorAlpha = new double;
        priorBeta = new double;
        value = new double;

        *context = _context;
        *I_star = _I_star;
        *E = _E;
        *A0 = _A0;
        *p_ei = _p_ei;
        *priorAlpha = _priorAlpha + 1;
        *priorBeta = _priorBeta + 1;
        *value = -1.0;

    }
    FC_P_EI::~FC_P_EI()
    {
        delete I_star;
        delete E;
        delete A0;
        delete p_ei;
        delete value;
        delete priorAlpha;
        delete priorBeta;
        delete context;
    }

    int FC_P_EI::evalCPU()
    { 
        *value = 0.0;
        int j;
        int nLoc = *((*A0) -> numLocations);
        int nTpts = *((*E) -> ncol);
        int i_star_sum = 0;
        int e_sum = 0;
        for (j =0; j<(nLoc*nTpts); j++)
        {
            i_star_sum += ((*I_star)-> data)[j];
            e_sum += ((*E)-> data)[j];
        }

        *value = dbeta(**p_ei, *priorAlpha + i_star_sum, *priorBeta - i_star_sum + e_sum); 
        return 0;
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

    int FC_P_EI::sampleCPU()
    {
        double a, b;
        a = ((*I_star) -> marginSum(3, -1));
        b = ((*E) -> marginSum(3, -1)) - a;
        //std::cout << "(a,b): (" << a << "," << b <<")\n";
        //std::cout << "(a,b): (" << a + *priorAlpha << "," << b+*priorBeta <<")\n";
        (**p_ei) = ((*context) -> random -> beta(a+*priorAlpha, b+*priorBeta));
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
    void FC_P_EI::setValue(double val)
    {
        *(this -> value) = val;
    }


    FC_P_IR::FC_P_IR(ModelContext *_context,
                     CompartmentalModelMatrix *_R_star,
                     CompartmentalModelMatrix *_I,
                     InitData *_A0,
                     double *_p_ir,
                     double _priorAlpha,
                     double _priorBeta)
    {

        context = new ModelContext*;
        R_star = new CompartmentalModelMatrix*;
        I = new CompartmentalModelMatrix*;
        A0 = new InitData*;
        p_ir = new double*;
        priorAlpha = new double;
        priorBeta = new double;
        value = new double;

        *context = _context;
        *R_star = _R_star;
        *I = _I;
        *A0 = _A0;
        *p_ir = _p_ir;
        *priorAlpha = _priorAlpha + 1;
        *priorBeta = _priorBeta + 1;
        *value = -1.0;

    }

    FC_P_IR::~FC_P_IR()
    {
        delete R_star;
        delete I;
        delete A0;
        delete p_ir;
        delete value;
        delete priorAlpha;
        delete priorBeta;
        delete context;
    }

    int FC_P_IR::evalCPU()
    {
        *value = 0.0;
        int j;
        int nLoc = *((*A0) -> numLocations);
        int nTpts = *((*I) -> ncol);
        int r_star_sum = 0;
        int i_sum = 0;
        for (j =0; j<(nLoc*nTpts); j++)
        {
            r_star_sum += ((*R_star)-> data)[j];
            i_sum += ((*I)-> data)[j];
        }

        *value = dbeta(**p_ir, *priorAlpha + r_star_sum, *priorBeta - r_star_sum + i_sum); 
        return 0;
    }

    int FC_P_IR::evalOCL()
    {
        //NOT IMPLEMENTED
        return -1;
    }
    int FC_P_IR::calculateRelevantCompartments()
    {
        // Not used, do nothing. 
        return(0);
    }

    int FC_P_IR::sampleCPU()
    {
        double a,b;
        a = (*R_star) -> marginSum(3,-1);
        b = ((*I) -> marginSum(3,-1)) - a;
        (**p_ir) = ((*context)->random->beta(a+(*priorAlpha), b+(*priorBeta)));
        return(0);
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
    void FC_P_IR::setValue(double val)
    {
        *(this -> value) = val;
    }


}



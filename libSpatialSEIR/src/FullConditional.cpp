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
    int matMult(double* output, double * A, double * B, int Arow, int Acol, 
            int Brow, int Bcol, bool TransA, bool TransB, int ldA, int ldB, int ldC)
    {
        // Use BLAS to matrix multiply, assume column major, non transposing.
        // Double check this code when I have internet access. 

        if ((TransA ? Arow : Acol) != (TransB ? Bcol : Brow))
        {
            std::cerr << "Invalid Matrix Dimensions: " << std::endl;
            std::cerr << "A: " <<  Arow << " x " << Acol << std::endl;
            std::cerr << "B: " <<  Brow << " x " << Bcol << std::endl;
            std::cerr << "Transpose: " << TransA << ", " << TransB << std::endl;
            throw(-1);
        }
        cblas_dgemm(CblasColMajor,
                    TransA ? CblasTrans : CblasNoTrans,
                    TransB ? CblasTrans : CblasNoTrans,
                    Arow, Bcol, Brow,
                    1.0, 
                    A, ldA,  
                    B, ldB,  
                    0.0, output, ldC);
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
                       int *nLoc)
    {
        this -> populate(_S0, _E0, _I0, _R0, 
                nLoc);
    }

    InitData::InitData()
    {
        // Do nothing
    }

    void InitData::populate(int *_S0, 
                       int *_E0,
                       int *_I0,
                       int *_R0,
                       int *nLoc
                       )
    {
        S0 = new int[*nLoc];
        E0 = new int[*nLoc];
        I0 = new int[*nLoc];
        R0 = new int[*nLoc];
        numLocations = new int;
        *numLocations = *nLoc;
        int i;
        for (i = 0; i < *nLoc; i++)
        {
            S0[i] = _S0[i]; 
            E0[i] = _E0[i];
            I0[i] = _I0[i];
            R0[i] = _R0[i];
        }
    }

    InitData::~InitData()
    {
        delete S0;
        delete E0;
        delete I0;
        delete R0;
        delete numLocations;
    }


    int CompartmentFullConditional::sampleCompartment_CPU(ModelContext* context,
                                                       CompartmentalModelMatrix* starCompartment,
                                                       double width)
    {
        int i;
        int nLoc = *(starCompartment -> ncol);        
        // Main loop: 
        for (i = 0; i < nLoc; i++)
        {

            sampleCompartmentLocation(i, context, starCompartment, width);
            //std::cout << "(i,val): (" << i << ", " << this->getValue() << ")\n";
        }

        return(0);
    }

    /*
    int CompartmentFullConditional::sampleCompartment_CPU(ModelContext* context,
                                                       CompartmentalModelMatrix* starCompartment,
                                                       double width)
    {

        int i;
        int nLoc = *(starCompartment -> ncol);        
        int nTpts = *(starCompartment -> nrow);
        int batchSize = 5;
        int numBatches = nTpts % batchSize; 
        int* batchCache = new int[batchSize];
        // Main loop: 
        for (i = 0; i < nLoc; i++)
        {
            jointSampleCompartmentLocation(i,batchSize,numBatches,batchCache,context,starCompartment,width);
            //sampleCompartmentLocation(i, context, starCompartment, width);
            //std::cout << "(i,val): (" << i << ", " << this->getValue() << ")\n";
        }
        delete batchCache;
        return(0);
    }
    */

    int CompartmentFullConditional::sampleCompartmentLocation(int i, ModelContext* context,
                                                       CompartmentalModelMatrix* starCompartment,
                                                       double width)
    {
        int j, compIdx;
        int nTpts = *(starCompartment -> nrow);
        int x0, x1;
        double initVal, newVal;
        double initProposal, newProposal;
        double criterion;
 
        compIdx = i*nTpts;
        for (j = 0; j < nTpts; j ++)
        { 
            //std::cout << j << "\n";
            this -> calculateRelevantCompartments(i,j); 
            this -> evalCPU(i,j);
            x0 = (starCompartment -> data)[compIdx];
            initVal = (this -> getValue());

            // Propose new value, bounded away from zero. 
            x1 = std::floor(std::max(0.0,(context -> random -> normal(x0,width))));
            (starCompartment -> data)[compIdx] = x1;
            this -> calculateRelevantCompartments(i,j);
            this -> evalCPU(i,j);
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


            if (!std::isfinite(this -> getValue()))
            {
                std::cout << "Impossible value selected:\n";
                std::cout << "(i,j): (" << i << "," << j << ")\n";
                std::cout << "Data value: " << (starCompartment -> data)[compIdx] << "\n";
                this -> printDebugInfo(i,j);
                throw(-1);
            }
            compIdx ++;
        }

        return(0);
    }

    int CompartmentFullConditional::jointSampleCompartmentLocation(int i, int batchSize, int numBatches, int* batchCache, ModelContext* context,
                                                                   CompartmentalModelMatrix* starCompartment,
                                                                   double width)
    {
        int j,k,l,compIdx;
        int nTpts = *(starCompartment -> nrow);
        int x0, x1;
        double initVal, newVal;
        double initProposal, newProposal;
        double criterion;
        //int numBatches = nTpts % batchSize; 
 
        compIdx = i*nTpts;
        j = 0;
        for (k = 0; k < numBatches; k++)
        { 
            //init
            newProposal = 0.0;
            initProposal = 0.0;
            criterion = 0.0;

            this -> calculateRelevantCompartments(i,j); 
            this -> evalCPU(i,j);

            initVal = (this -> getValue());

            // Propose new values, bounded away from zero. 
            for (l = 0; l < batchSize; l++)
            { 
                x0 = (starCompartment -> data)[compIdx + l];
                batchCache[l] = x0;
                x1 = std::floor(std::max(0.0,(context -> random -> normal(x0,width))));
                (starCompartment -> data)[compIdx + l] = x1;
                newProposal += (context -> random -> dnorm(x1, x0,width));
                initProposal += (context -> random -> dnorm(x0, x1,width));
            }
            this -> calculateRelevantCompartments(i,j);
            this -> evalCPU(i,j);
            newVal = (this->getValue());
            criterion = (newVal - initVal) + (initProposal - newProposal);
            if (std::log((context -> random -> uniform())) < criterion)
            {
                // Accept new value
            }
            else
            {
                // Keep Original Value
                for (l = 0; l < batchSize; l++)
                {
                    (starCompartment -> data)[compIdx + l] = batchCache[l];
                }
                this -> calculateRelevantCompartments(i,j);
                this -> setValue(initVal); 
            }                


            if (!std::isfinite(this -> getValue()))
            {
                std::cout << "Impossible value selected:\n";
                std::cout << "(i,j): (" << i << "," << j << ")\n";
                std::cout << "Data value: " << (starCompartment -> data)[compIdx] << "\n";
                this -> printDebugInfo(i,j);
                throw(-1);
            }
            compIdx += batchSize;
            j += batchSize;
        }
        
        //Remainder
        for (l = j; l < nTpts; l++)
        {
            this -> calculateRelevantCompartments(i,l); 
            this -> evalCPU(i,l);
            x0 = (starCompartment -> data)[compIdx];
            initVal = (this -> getValue());

            // Propose new value, bounded away from zero. 
            x1 = std::floor(std::max(0.0,(context -> random -> normal(x0,width))));
            (starCompartment -> data)[compIdx] = x1;
            this -> calculateRelevantCompartments(i,l);
            this -> evalCPU(i,l);
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
                this -> calculateRelevantCompartments(i,l);
                this -> setValue(initVal); 
            }                


            if (!std::isfinite(this -> getValue()))
            {
                std::cout << "Impossible value selected:\n";
                std::cout << "(i,l): (" << i << "," << l << ")\n";
                std::cout << "Data value: " << (starCompartment -> data)[compIdx] << "\n";
                this -> printDebugInfo(i,l);
                throw(-1);
            }
            compIdx ++;
        } 
        return(0);
    }


    int CompartmentFullConditional::sliceSampleCompartmentLocation(int i, ModelContext* context,
                                                       CompartmentalModelMatrix* starCompartment,
                                                       double width)
    {
        int j, compIdx;
        int nTpts = *(starCompartment -> nrow);
        double l,r,y,x,x0;
 
        compIdx = i*nTpts;
        for (j = 0; j < nTpts; j ++)
        { 
            //std::cout << j << "\n";
            this -> calculateRelevantCompartments(i,j); 
            this -> evalCPU(i,j);
            x = (starCompartment -> data)[compIdx];
            y = (this -> getValue()) - (context -> random -> gamma());
            l = std::max(0.0, (x-(context -> random -> uniform())*width));
            r = l + width;

            do
            {
                x0 = ((context -> random -> uniform())*(r-l) + l);
                (starCompartment -> data)[compIdx] = std::floor(x0);
                this -> calculateRelevantCompartments(i,j); 
                this -> evalCPU(i,j);
                if (x0 >= x){r=x0;}
                else{l=x0;}
            } while (y >= (this -> getValue()));
            compIdx ++;
        }
        return(0);
    }



    int InitCompartmentFullConditional::sampleCompartment_CPU(ModelContext* context, 
                                                              int* initCompartment,
                                                              double width)
    {

        int i;
        int nLoc = *(context -> S -> ncol);        
        // Main loop: 
        for (i = 0; i < nLoc; i++)
        {
            sampleCompartmentLocation(i, context, initCompartment, width);
            //std::cout << "(i,val): (" << i << ", " << this->getValue() << ")\n";
        }
        return(0);
    }


    int InitCompartmentFullConditional::sampleCompartmentLocation(int i, ModelContext* context, 
                                                                  int* initCompartment,
                                                                  double width)
    {
        int x0, x1;
        double initVal, newVal;
        double initProposal, newProposal;
        double criterion;
     
        this -> calculateRelevantCompartments(i); 
        this -> evalCPU(i);
        x0 = initCompartment[i];
        initVal = (this -> getValue());

        // Propose new value, bounded away from zero. 
        x1 = std::floor(std::max(0.0,(context -> random -> normal(x0,width))));
        initCompartment[i] = x1;
        this -> calculateRelevantCompartments(i);
        this -> evalCPU(i);
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
            initCompartment[i] = x0;
            this -> calculateRelevantCompartments(i);
            this -> setValue(initVal); 
        }                

        if (!std::isfinite(this -> getValue()))
        {
            std::cout << "Impossible value selected:\n";
            std::cout << "i: (" << i << "," << ")\n";
            this -> printDebugInfo(i);
            throw(-1);
        }
        return(0);

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
     * Implement Full Conditional for S0
     *
     */

    FC_S0::FC_S0(ModelContext* _context,
                 CompartmentalModelMatrix *_S,
                 CompartmentalModelMatrix *_E,
                 CompartmentalModelMatrix *_E_star,
                 CompartmentalModelMatrix *_I_star,
                 InitData *_A0,
                 double *_p_se,
                 double *_p_ei,
                 double _sliceWidth)
    {

        context = new ModelContext*;
        S = new CompartmentalModelMatrix*;
        E = new CompartmentalModelMatrix*;
        E_star = new CompartmentalModelMatrix*;
        I_star = new CompartmentalModelMatrix*;
        A0 = new InitData*;
        p_se = new double*;
        p_ei = new double*;
        sliceWidth = new double;
        value = new long double;

        *context = _context;
        *S = _S;
        *E = _E;
        *E_star = _E_star;
        *I_star = _I_star;
        *A0 = _A0;
        *p_se = _p_se;
        *p_ei = _p_ei;
        *sliceWidth = _sliceWidth;

    }
    FC_S0::~FC_S0()
    {

        delete context;
        delete S;
        delete E;
        delete I_star;
        delete E_star;
        delete p_ei;
        delete p_se;
        delete sliceWidth;
        delete value;

    }
    

    int FC_S0::evalCPU(int startLoc)
    {

        int j, compIdx;
        int nTpts = *((*S) -> nrow); 
        double p_se_val, p_ei_val;
        int S_val, E_val, Istar_val, Estar_val;
        long double output = 0.0;
        
        if (((*A0) -> S0)[startLoc] < 0 || 
            ((*A0) -> E0)[startLoc] < 0)
        {
            *value = -INFINITY;
            return(-1);
        }

        compIdx = startLoc*nTpts;
        p_ei_val = **p_ei;
        for (j = 0; j < nTpts; j++)
        {
            S_val = ((*S) -> data)[compIdx];
            E_val = ((*E) -> data)[compIdx];
            Istar_val = ((*I_star)-> data)[compIdx];
            Estar_val = ((*E_star) -> data)[compIdx];
            p_se_val = (*p_se)[compIdx];
            if (Estar_val > S_val || Istar_val > E_val)
            {
                *value = -INFINITY;
                return(-1);
            }
            else
            {
                output += (((*context) -> random -> dbinom(Estar_val, S_val, p_se_val)) +    
                            ((*context) -> random -> dbinom(Istar_val, E_val, p_ei_val)));

            }
            compIdx++;
        }

        if (!std::isfinite(output))
        {
            *value = -INFINITY;
            return(-1);
        }
        else
        {
            *value = output;
        }
    
       return 0;
    }

    int FC_S0::evalOCL()
    {
        // Not Implemented
        return(-1);
    }
    int FC_S0::sampleCPU()
    {
        sampleCompartment_CPU(*context, (*A0)->S0, *sliceWidth);
        return(0);
    
    }
    int FC_S0::sampleOCL()
    {
        // Not Implemented
        return(-1);
    }

    long double FC_S0::getValue()
    {
        return(*value);
    }

    void FC_S0::setValue(long double val)
    {
        *(this -> value) = val;
    }

    int FC_S0::calculateRelevantCompartments()
    {
        (*context) -> calculateS_CPU();
        (*context) -> calculateE_givenI_CPU(); 
        return(0);
    }

    int FC_S0::calculateRelevantCompartments(int startLoc)
    {
        (*context) -> calculateS_CPU(startLoc,0);
        (*context) -> calculateE_givenI_CPU(startLoc,0);
        return(0);
    }

    void FC_S0::printDebugInfo(int startLoc)
    {
        std::cout << "Error Sampling S0, location: " << startLoc << ", value: " << ((*A0) -> S0)[startLoc] << "\n";
        int j, compIdx;
        int nTpts = *((*S) -> nrow); 
        double p_se_val, p_ei_val;
        int S_val, E_val, Istar_val, Estar_val;
        long double output = 0.0;
        
        if (((*A0) -> S0)[startLoc] < 0 || 
            ((*A0) -> E0)[startLoc] < 0)
        {
            std::cout << "Invalid Value.\n";
            return;
        }

        compIdx = startLoc*nTpts;
        p_ei_val = **p_ei;
        for (j = 0; j < nTpts; j++)
        {
            S_val = ((*S) -> data)[compIdx];
            E_val = ((*E) -> data)[compIdx];
            Istar_val = ((*I_star)-> data)[compIdx];
            Estar_val = ((*E_star) -> data)[compIdx];
            p_se_val = (*p_se)[compIdx];
            if (Estar_val > S_val || Istar_val > E_val)
            {
                std::cout << "Bounds error detected, time point: " << j << "\n";
                std::cout << "S: " << S_val << "\n";
                std::cout << "E: " << E_val << "\n";
                std::cout << "E_star: " << Estar_val << "\n";
                std::cout << "I_star: " << Istar_val << "\n";
                return;
            }
            else
            {
                output += (((*context) -> random -> dbinom(Estar_val, S_val, p_se_val)) +    
                            ((*context) -> random -> dbinom(Istar_val, E_val, p_ei_val)));

                if (!std::isfinite(output))
                {
                    std::cout << "Calculation Error Detected, time point: " << j << "\n";
                    return;
                }
            }
            compIdx++;
        } 
       return;
    }

    /*
     *
     * Implement full conditional for E0
     *
     */
    

    FC_E0::FC_E0(ModelContext* _context, 
                 CompartmentalModelMatrix *_S,
                 CompartmentalModelMatrix *_E,
                 CompartmentalModelMatrix *_I,
                 CompartmentalModelMatrix *_E_star,
                 CompartmentalModelMatrix *_I_star,
                 CompartmentalModelMatrix *_R_star,
                 InitData *_A0,
                 double *_p_ir,
                 double *_p_ei,
                 double *_p_se,
                 double _sliceWidth)
    {
        context = new ModelContext*;
        S = new CompartmentalModelMatrix*;
        E = new CompartmentalModelMatrix*;
        I = new CompartmentalModelMatrix*;
        E_star = new CompartmentalModelMatrix*;
        I_star = new CompartmentalModelMatrix*;
        R_star = new CompartmentalModelMatrix*;
        A0 = new InitData*;
        p_se = new double*;
        p_ir = new double*;
        p_ei = new double*;
        sliceWidth = new double;
        value = new long double;

        *context = _context;
        *S = _S;
        *E = _E;
        *I = _I;
        *E_star = _E_star;
        *I_star = _I_star;
        *R_star = _R_star;
        *p_se = _p_se;
        *p_ir = _p_ir;
        *p_ei = _p_ei;
        *A0 = _A0;
        *sliceWidth = _sliceWidth;
    }
    FC_E0::~FC_E0()
    {
        delete context;
        delete S;
        delete E;
        delete I;
        delete E_star;
        delete I_star;
        delete R_star;
        delete p_se;
        delete p_ir;
        delete p_ei;
        delete A0;
        delete sliceWidth;
        delete value;
    }
    
    int FC_E0::evalCPU(int startLoc)
    {
        int i,j,compIdx,S_val,E_val,I_val,Istar_val,Estar_val,Rstar_val;
        double p_ei_val, p_ir_val, p_se_val;
        int nLoc = *((*E)->ncol);
        int nTpts = *((*E)->nrow);
        long double output = 0.0;

        if (((*A0) -> E0)[startLoc] < 0 || 
            ((*A0) -> I0)[startLoc] < 0)
        {
            *value = -INFINITY;
            return(-1);
        }

        compIdx = startLoc*nTpts;
        p_ei_val = **p_ei;
        p_ir_val = **p_ir;
        for (i = 0; i < nTpts; i++)
        {
                Rstar_val = ((*R_star)->data)[compIdx]; 
                Istar_val = ((*I_star)->data)[compIdx];
                E_val = ((*E)->data)[compIdx];
                I_val = ((*I)->data)[compIdx];
                if (Istar_val > E_val ||
                    Rstar_val > I_val)
                {
                    *value = -INFINITY;
                    return(-1);
                }
                else
                { 
                    output += (((*context) -> random -> dbinom(Rstar_val, I_val, p_ir_val)) + 
                               ((*context) -> random -> dbinom(Istar_val, E_val, p_ei_val)));
                }
                compIdx ++; 
        }

        // p_se changes, so need to look at p_se component for all locations and 
        // time points after 0
        for (i = 0; i < nLoc; i++)
        {
            compIdx = i*nTpts;
            for (j = 0; j< nTpts; j++)
            {

                p_se_val = (*p_se)[compIdx];
                Estar_val = ((*E_star) -> data)[compIdx];
                S_val = ((*S)->data)[compIdx];
                if (p_se_val > 1 || p_se_val < 0)
                {
                    *value = -INFINITY;
                    return(-1);
                }

                output += (*context) -> random -> dbinom(Estar_val,S_val, p_se_val);
                compIdx ++; 
            }
        }


        if (!std::isfinite(output))
        {
            *value = -INFINITY;
            return(-1);
        }
        else
        {
            *value = output;
        }
        return(0);
    }

    int FC_E0::evalOCL()
    {
        // Not Implemented
        return(-1);
    }
    int FC_E0::sampleCPU()
    {
        sampleCompartment_CPU(*context, (*A0) -> E0, *sliceWidth);
        return(0);
    }
    int FC_E0::sampleOCL()
    {
        // Not Implemented
        return(-1);
    }

    long double FC_E0::getValue()
    {
        return(*value);
    }

    void FC_E0::setValue(long double val)
    {
        *(this -> value) = val;
    }

    int FC_E0::calculateRelevantCompartments()
    {
        (*context) -> calculateE_CPU();
        (*context) -> calculateI_givenR_CPU();
        (*context) -> calculateP_SE_CPU();
        return(0);
    }

    int FC_E0::calculateRelevantCompartments(int startLoc)
    {
        (*context) -> calculateE_CPU(startLoc, 0);
        (*context) -> calculateI_givenR_CPU(startLoc,0);
        (*context) -> calculateP_SE_CPU(startLoc,0);
        return(0);
    }

    void FC_E0::printDebugInfo(int loc)
    {
        std::cout << "Error Sampling E0, location: " << loc << ", value: " << ((*A0) -> E0)[loc] << "\n";
    }

    /*
     *
     * Implement full conditional for I0
     *
     */
    
    FC_I0::FC_I0(ModelContext* _context, 
                 CompartmentalModelMatrix *_S,
                 CompartmentalModelMatrix *_I,
                 CompartmentalModelMatrix *_R,
                 CompartmentalModelMatrix *_S_star,
                 CompartmentalModelMatrix *_E_star,
                 CompartmentalModelMatrix *_R_star,
                 InitData *_A0,
                 double *_p_ir,
                 double *_p_rs,
                 double *_p_se,
                 double _sliceWidth)
    {
        context = new ModelContext*;
        S = new CompartmentalModelMatrix*;
        I = new CompartmentalModelMatrix*;
        R = new CompartmentalModelMatrix*;
        S_star = new CompartmentalModelMatrix*;
        E_star = new CompartmentalModelMatrix*;
        R_star = new CompartmentalModelMatrix*;
        A0 = new InitData*;
        p_ir = new double*;
        p_rs = new double*;
        p_se = new double*;
        sliceWidth = new double;
        value = new long double;

        *context = _context;
        *S = _S;
        *I = _I;
        *R = _R;
        *S_star = _S_star;
        *E_star = _E_star;
        *R_star = _R_star;
        *A0 = _A0;
        *p_ir = _p_ir;
        *p_se = _p_se;
        *p_rs = _p_rs;
        *sliceWidth = _sliceWidth;
    }
    FC_I0::~FC_I0()
    {
        delete context;
        delete S;
        delete I;
        delete R;
        delete S_star;
        delete E_star; 
        delete R_star;
        delete A0;
        delete sliceWidth;
        delete p_ir;
        delete p_se;
        delete p_rs;
        delete value;
    }
    

    int FC_I0::evalCPU(int startLoc)
    {

        int i,j, compIdx;
        int nTpts = *((*R) -> nrow);
        int nLoc = *((*R) -> ncol);

        long double output = 0.0;
        
        double p_se_val;
        double p_rs_val;
        double ln_1m_p_ir = std::log(1-(**p_ir));
        int Rstar_val, Sstar_val, Estar_val, R_val, I_val, S_val;   

        if (((*A0) -> I0)[startLoc] < 0 || 
            ((*A0) -> R0)[startLoc] < 0)
        {
            *value = -INFINITY;
            return(-1);
        }
        compIdx = startLoc*nTpts;

        // Is p_rs meaningful?
        if ((*context) -> config -> reinfectionMode <= 2)
        {

            for (j = 0; j < nTpts; j++)
            {
                Rstar_val = ((*R_star) -> data)[compIdx];
                Sstar_val = ((*S_star)->data)[compIdx];
                R_val = ((*R) ->data)[compIdx];
                I_val = ((*I) ->data)[compIdx];
                p_rs_val = (*p_rs)[j];

                if (Rstar_val > I_val || 
                        Sstar_val > R_val)
                {
                    *value = -INFINITY;
                    return(-1);
                }
                else
                {
                    output +=  (ln_1m_p_ir*(I_val) +
                                std::log(1-p_rs_val)*(R_val) +
                                ((*context) -> random -> choosePartial(R_val, Sstar_val)) +
                                ((*context) -> random -> choose(I_val, Rstar_val)));
                }
                compIdx++;
            } 
        }
        else
        {
            for (j = 0; j < nTpts; j++)
            {
                Rstar_val = ((*R_star) -> data)[compIdx];
                R_val = ((*R) ->data)[compIdx];
                I_val = ((*I) ->data)[compIdx];

                if (Rstar_val > I_val)
                {
                    *value = -INFINITY;
                    return(-1);
                }
                else
                {
                    output +=  (ln_1m_p_ir*(I_val) +
                                ((*context) -> random -> choose(I_val, Rstar_val)));
                }
                compIdx++;
            } 

        }

        // p_se changes, so need to look at p_se component for all locations and 
        // time points after 0
        for (i = 0; i < nLoc; i++)
        {
            compIdx = i*nTpts;
            for (j = 0; j< nTpts; j++)
            {

                p_se_val = (*p_se)[compIdx];
                Estar_val = ((*E_star) -> data)[compIdx];
                S_val = ((*S)->data)[compIdx];
                if (p_se_val > 1 || p_se_val < 0)
                {
                    *value = -INFINITY;
                    return(-1);
                }

                output += (*context) -> random -> dbinom(Estar_val,S_val, p_se_val);
                compIdx ++; 
            }
        }

        if (!std::isfinite(output))
        {
            *value = -INFINITY;   
            return(-1);
        }
        else
        {
            *value = output;
        }

        return 0;
    }

    int FC_I0::evalOCL()
    {
        // Not Implemented
        return(-1);
    }
    int FC_I0::sampleCPU()
    {
        sampleCompartment_CPU(*context, (*A0) -> I0, *sliceWidth);
        return(0);
    }
    int FC_I0::sampleOCL()
    {
        // Not Implemented
        return(-1);
    }

    long double FC_I0::getValue()
    {
        return(*value);
    }

    void FC_I0::setValue(long double val)
    {
        *(this -> value) = val;
    }

    int FC_I0::calculateRelevantCompartments()
    {
        (*context) -> calculateI_CPU();
        (*context) -> calculateR_givenS_CPU();
        (*context) -> calculateP_SE_CPU();
        return(0);
    }

    int FC_I0::calculateRelevantCompartments(int startLoc)
    {
        (*context) -> calculateI_CPU(startLoc, 0);
        (*context) -> calculateR_givenS_CPU(startLoc, 0);
        (*context) -> calculateP_SE_CPU(startLoc, 0);
        return(0);
    }

    void FC_I0::printDebugInfo(int loc)
    {
        std::cout << "Error Sampling I0, location: " << loc << ", value: " << ((*A0) -> I0)[loc] << "\n";
    }

    /*
     *
     *
     * Implement FC for R0
     * 
     *
     */

    FC_R0::FC_R0(ModelContext* _context,   
                 CompartmentalModelMatrix *_R,
                 CompartmentalModelMatrix *_S,
                 CompartmentalModelMatrix *_S_star,
                 CompartmentalModelMatrix *_E_star,
                 CompartmentalModelMatrix *_R_star,
                 InitData *_A0,
                 double *_p_rs,
                 double *_p_se,
                 double _sliceWidth)
    {
        context = new ModelContext*;
        S = new CompartmentalModelMatrix*;
        R = new CompartmentalModelMatrix*;
        S_star = new CompartmentalModelMatrix*;
        E_star = new CompartmentalModelMatrix*;
        R_star = new CompartmentalModelMatrix*;
        A0 = new InitData*;
        p_rs = new double*;
        p_se = new double*;
        sliceWidth = new double;
        value = new long double;

        *context = _context;
        *S = _S;
        *R = _R;
        *S_star = _S_star;
        *E_star = _E_star;
        *R_star = _R_star;
        *A0 = _A0;
        *p_se = _p_se;
        *p_rs = _p_rs;
        *sliceWidth = _sliceWidth;
    }
    FC_R0::~FC_R0()
    {
        delete context;
        delete S;
        delete R;
        delete S_star;
        delete E_star; 
        delete R_star;
        delete A0;
        delete sliceWidth;
        delete p_se;
        delete p_rs;
        delete value;
    }
    

    int FC_R0::evalCPU(int startLoc)
    {

        int j, compIdx;
        int nTpts = *((*R) -> nrow);

        long double output = 0.0;
        
        double p_se_val;
        double p_rs_val;
        int Sstar_val, Estar_val, R_val, S_val;   

        if (((*A0) -> S0)[startLoc] < 0 || 
            ((*A0) -> R0)[startLoc] < 0)
        {
            *value = -INFINITY;
            return(-1);
        }

        compIdx = startLoc*nTpts;

        // Is p_rs meaningful?
        if ((*context) -> config -> reinfectionMode <= 2)
        {
            for (j = 0; j < nTpts; j++)
            {
                Estar_val = ((*E_star) -> data)[compIdx];
                Sstar_val = ((*S_star)->data)[compIdx];
                R_val = ((*R) ->data)[compIdx];
                S_val = ((*S) ->data)[compIdx];
                p_rs_val = (*p_rs)[j];

                if (Estar_val > S_val || 
                        Sstar_val > R_val)
                {
                    *value = -INFINITY;
                    return(-1);
                }
                else
                {
                    output += (((*context) -> random -> dbinom(Estar_val, S_val, p_se_val)) + 
                               ((*context) -> random -> dbinom(Sstar_val, R_val, p_rs_val)));

                }
                compIdx++;
            } 
        }
        else 
        {
            for (j = 0; j < nTpts; j++)
            {
                Estar_val = ((*E_star) -> data)[compIdx];
                R_val = ((*R) ->data)[compIdx];
                S_val = ((*S) ->data)[compIdx];
                p_rs_val = (*p_rs)[j];

                if (Estar_val > S_val)
                {
                    *value = -INFINITY;
                    return(-1);
                }
                else
                {
                    output += (((*context) -> random -> dbinom(Estar_val, S_val, p_se_val)));

                }
                compIdx++;
            } 

        }

        if (!std::isfinite(output))
        {
            *value = -INFINITY;   
            return(-1);
        }
        else
        {
            *value = output;
        }

        return 0;
    }

    int FC_R0::evalOCL()
    {
        // Not Implemented
        return(-1);
    }
    int FC_R0::sampleCPU()
    {
        sampleCompartment_CPU(*context, (*A0) -> R0, *sliceWidth);
        return(0);
    }
    int FC_R0::sampleOCL()
    {
        // Not Implemented
        return(-1);
    }

    long double FC_R0::getValue()
    {
        return(*value);
    }

    void FC_R0::setValue(long double val)
    {
        *(this -> value) = val;
    }

    int FC_R0::calculateRelevantCompartments()
    {
        (*context) -> calculateR_CPU();
        (*context) -> calculateS_givenE_CPU();
        return(0);
    }

    int FC_R0::calculateRelevantCompartments(int startLoc)
    {
        (*context) -> calculateR_CPU(startLoc, 0);
        (*context) -> calculateS_givenE_CPU(startLoc,0);
        return(0);
    }

    void FC_R0::printDebugInfo(int loc)
    {
        std::cout << "Error Sampling R0, location: " << loc << ", value: " << ((*A0) -> R0)[loc] << "\n";
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
       value = new long double;
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

    int FC_S_Star::evalCPU(int startLoc, int startTime)
    {
        int i,compIdx,Sstar_val,Estar_val,S_val,R_val;
        double p_se_val, p_rs_val;
        int nTpts = *((*S)->nrow);
        long double output = 0.0;
        long unsigned int S_star_sum;
        long unsigned int R_star_sum;
        int64_t aDiff;

        compIdx = startLoc*nTpts + startTime;
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
                    *value = -INFINITY;
                    return(-1);
                }
                else
                { 
                    output += (((*context) -> random -> dbinom(Sstar_val, R_val, p_rs_val)) + 
                               ((*context) -> random -> dbinom(Estar_val, S_val, p_se_val)));
                }
                compIdx ++; 
        }

        S_star_sum = (*S_star)->marginSum(2,startLoc);
        R_star_sum = (*R_star)->marginSum(2,startLoc);
        aDiff = (S_star_sum > R_star_sum ? S_star_sum - R_star_sum : R_star_sum - S_star_sum)/nTpts;
        output -= (aDiff*aDiff)*(*steadyStateConstraintPrecision);

        if (!std::isfinite(output))
        {
            *value = -INFINITY;
            return(-1);
        }
        else
        {
            *value = output;
        }
        return(0);
    }

    void FC_S_Star::printDebugInfo(int startLoc, int startTime)
    {
        std::cout << "S_star debug info, location " << startLoc << ", time " << startTime << "\n";     

        int i,compIdx,Sstar_val,Estar_val,S_val,R_val;
        double p_se_val, p_rs_val;
        int nTpts = *((*S)->nrow);
        long double output = 0.0;
        long unsigned int S_star_sum;
        long unsigned int R_star_sum;
        int64_t aDiff;

        compIdx = startLoc*nTpts + startTime;
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
                    std::cout << "Bounds Error Detected at time " << i << "\n";
                    std::cout << "S_star: " << Sstar_val << "\n";
                    std::cout << "E_star: " << Estar_val << "\n";
                    std::cout << "R: " << R_val << "\n";
                    std::cout << "S: " << S_val << "\n";
                    return;
                }
                else
                { 
                    output += (((*context) -> random -> dbinom(Sstar_val, R_val, p_rs_val)) + 
                               ((*context) -> random -> dbinom(Estar_val, S_val, p_se_val)));
                }
                if (! std::isfinite(output))
                {
                    std::cout << "Calculation Error Detected at time " << i << "\n";
                    std::cout << "S_star: " << Sstar_val << "\n";
                    std::cout << "E_star: " << Estar_val << "\n";
                    std::cout << "R: " << R_val << "\n";
                    std::cout << "S: " << S_val << "\n";
                    std::cout << "p_se: " << p_se_val << "\n";
                    std::cout << "p_rs: " << p_rs_val << "\n"; 
                    return;
                }
                compIdx ++; 
        }

        S_star_sum = (*S_star)->marginSum(2,startLoc);
        R_star_sum = (*R_star)->marginSum(2,startLoc);
        aDiff = (S_star_sum > R_star_sum ? S_star_sum - R_star_sum : R_star_sum - S_star_sum)/nTpts;
        output -= (aDiff*aDiff)*(*steadyStateConstraintPrecision);

        if (!std::isfinite(output))
        {
            std::cout << "Combinatorics Error Detected\n";
            std::cout << "S_star: " << Sstar_val << "\n";
            std::cout << "E_star: " << Estar_val << "\n";
            std::cout << "R: " << R_val << "\n";
            std::cout << "S: " << S_val << "\n";
            std::cout << "p_se: " << p_se_val << "\n";
            std::cout << "p_rs: " << p_rs_val << "\n"; 
            std::cout << "S_star_sum: " << S_star_sum << "\n";
            std::cout << "R_star_sum: " << R_star_sum << "\n";
            return;
        }
        return;

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
        this -> sampleCompartment_CPU(*context,
                                  *S_star,*sliceWidth);
        return 0;
    }
    int FC_S_Star::sampleOCL()
    {
        //NOT IMPLEMENTED
        return -1;
    }
    long double FC_S_Star::getValue()
    {
        return(*(this -> value));
    }
    void FC_S_Star::setValue(long double val)
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
        value = new long double;
        
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

    int FC_E_Star::evalCPU(int startLoc, int startTime)
    {
        int j, compIdx;
        int nTpts = *((*S) -> nrow); 
        double ln_1m_p_ei = std::log(1-(**p_ei));    
        double p_se_val;
        int S_val, E_val, Estar_val, Istar_val;
        long double output = 0.0;
        long unsigned int E_star_sum;
        long unsigned int I_star_sum;
        int64_t aDiff; 

        compIdx = startLoc*nTpts + startTime;
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
                *value = -INFINITY;
                return(-1);
            }
            else
            {
                output += (((*context) -> random -> dbinom(Estar_val, S_val, p_se_val)) +    
                             ln_1m_p_ei*E_val + 
                             ((*context) -> random -> choosePartial(E_val, Istar_val)));
            }
            compIdx++;
        }

        E_star_sum = (*E_star)->marginSum(2,startLoc);
        I_star_sum = (*I_star)->marginSum(2,startLoc);
        aDiff = (E_star_sum > I_star_sum ? E_star_sum - I_star_sum : I_star_sum - E_star_sum)/nTpts;
        output -= (aDiff*aDiff)*(*steadyStateConstraintPrecision);

        if (!std::isfinite(output))
        {
            *value = -INFINITY;
            return(-1);
        }
        else
        {
            *value = output;
        }
    
       return 0;
    }

    void FC_E_Star::printDebugInfo(int loc, int tpt)
    {
        std::cout << "E_star debug info, location " << loc << ", time " << tpt << "\n";     
        std::cout << "...looking for problems\n";

        int j, compIdx;
        int nTpts = *((*S) -> nrow); 
        double ln_1m_p_ei = std::log(1-(**p_ei));    
        double p_se_val;
        int S_val, E_val, Estar_val, Istar_val;
        long double output = 0.0;
        long unsigned int E_star_sum;
        long unsigned int I_star_sum;
        int64_t aDiff; 


        compIdx = loc*nTpts + tpt;
        for (j = tpt; j < nTpts; j++)
        {
            std::cout << "time " << j << "\n";
            Estar_val = ((*E_star) -> data)[compIdx];
            S_val = ((*S) -> data)[compIdx];
            E_val = ((*E) -> data)[compIdx];
            Istar_val = ((*I_star) -> data)[compIdx];
            p_se_val = (*p_se)[compIdx];

            if (Estar_val < 0 || Estar_val > S_val || 
                    Istar_val > E_val)

            {
                std::cout << "Bounds Error Detected:\n";
                std::cout << "E_star: " << Estar_val << "\n";
                std::cout << "S: " << S_val << "\n";
                std::cout << "I_star: " << Istar_val << "\n";
                return;
            }
            else
            {
                output = (std::log(p_se_val)*Estar_val +
                             std::log(1-p_se_val)*(S_val - Estar_val) +
                             ln_1m_p_ei*E_val + 
                             ((*context) -> random -> choosePartial(E_val, Istar_val)) + 
                             ((*context)->random->choose(S_val, Estar_val)));
                if (!std::isfinite(output))
                {
                    std::cout << "Computation Error Detected: " << tpt << "\n";
                    std::cout << "E_star: " << Estar_val << "\n";
                    std::cout << "S: " << S_val << "\n";
                    std::cout << "E: " << E_val << "\n";
                    std::cout << "I_star: " << Istar_val << "\n";
                    std::cout << "p_se: " << p_se_val << "\n";
                    std::cout << "ln_1m_p_ei: " << ln_1m_p_ei << "\n";
                    return; 
                }
            }
            compIdx++;
        }

        E_star_sum = (*E_star)->marginSum(2,loc);
        I_star_sum = (*I_star)->marginSum(2,loc);
        aDiff = (E_star_sum > I_star_sum ? E_star_sum - I_star_sum : I_star_sum - E_star_sum)/nTpts;
        output -= (aDiff*aDiff)*(*steadyStateConstraintPrecision);

        if (!std::isfinite(output))
        {
            std::cout << "Combinatorics Error Detectd:\n";
            std::cout << "Computation Error Detected:\n";
            std::cout << "E_star: " << Estar_val << "\n";
            std::cout << "S: " << S_val << "\n";
            std::cout << "E: " << E_val << "\n";
            std::cout << "I_star: " << Istar_val << "\n";
            std::cout << "p_se: " << p_se_val << "\n";
            std::cout << "ln_1m_p_ei: " << ln_1m_p_ei << "\n";
            std::cout << "E_star_sum: " << E_star_sum << "\n";
            std::cout << "I_star_sum: " << I_star_sum << "\n";
            return;
        }    
       return;
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
        this -> sampleCompartment_CPU(*context,
                                  *E_star,*sliceWidth);
        return 0;
    }

    int FC_E_Star::sampleOCL()
    {
        //NOT IMPLEMENTED
        return -1;
    }

    long double FC_E_Star::getValue()
    {
        return(*(this -> value));
    }
    void FC_E_Star::setValue(long double val)
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
        value = new long double;

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

    int FC_R_Star::evalCPU(int startLoc, int startTime)
    {
        int i,j, compIdx;
        int nTpts = *((*R) -> nrow);
        int nLoc = *((*R) -> ncol);

        long double output = 0.0;
        
        double p_se_val;
        double p_rs_val;
        double ln_p_ir = std::log(**p_ir);
        double ln_1m_p_ir = std::log(1-(**p_ir));
        int Rstar_val, Sstar_val, Estar_val, R_val, I_val, S_val;   
        long unsigned int I_star_sum;
        long unsigned int R_star_sum;
        int64_t aDiff; 

        compIdx = startLoc*nTpts + startTime;
        // Is p_rs meaningful?
        if ((*context) -> config -> reinfectionMode <= 2)
        {
            for (j = startTime; j < nTpts; j++)
            {
                Rstar_val = ((*R_star) -> data)[compIdx];
                Sstar_val = ((*S_star)->data)[compIdx];
                R_val = ((*R) ->data)[compIdx];
                I_val = ((*I) ->data)[compIdx];
                p_rs_val = (*p_rs)[j];

                if (Rstar_val < 0 || Rstar_val > I_val || 
                        Sstar_val > R_val)
                {
                    *value = -INFINITY;
                    return(-1);
                }
                else
                {
                    output += (ln_p_ir*Rstar_val +
                                ln_1m_p_ir*(I_val - Rstar_val) +
                                std::log(1-p_rs_val)*(R_val) +
                                ((*context) -> random -> choosePartial(R_val, Sstar_val)) +
                                ((*context) -> random -> choose(I_val, Rstar_val)));
                }
                compIdx++;
            } 
        }
        else
        {
            for (j = startTime; j < nTpts; j++)
            {
                Rstar_val = ((*R_star) -> data)[compIdx];
                R_val = ((*R) ->data)[compIdx];
                I_val = ((*I) ->data)[compIdx];

                if (Rstar_val < 0 || Rstar_val > I_val)
                {
                    *value = -INFINITY;
                    return(-1);
                }
                else
                {
                    output += (ln_p_ir*Rstar_val +
                                ln_1m_p_ir*(I_val - Rstar_val) +
                                ((*context) -> random -> choose(I_val, Rstar_val)));
                }
                compIdx++;
            } 

        }

        // p_se changes, so need to look at p_se component for all locations and 
        // time points after startTime
        for (i = 0; i < nLoc; i++)
        {
            compIdx = i*nTpts + startTime;
            for (j = startTime; j< nTpts; j++)
            {

                p_se_val = (*p_se)[compIdx];
                Estar_val = ((*E_star) -> data)[compIdx];
                S_val = ((*S)->data)[compIdx];
                if (p_se_val > 1 || p_se_val < 0)
                {
                    *value = -INFINITY;
                    return(-1);
                }

                output += (*context) -> random -> dbinom(Estar_val,S_val, p_se_val);
                compIdx ++; 
            }
        }

        I_star_sum = (*I_star)->marginSum(2,startLoc);
        R_star_sum = (*R_star)->marginSum(2,startLoc);
        aDiff = (I_star_sum > R_star_sum ? I_star_sum - R_star_sum : R_star_sum - I_star_sum)/nTpts;
        output -= (aDiff*aDiff)*(*steadyStateConstraintPrecision);

        if (!std::isfinite(output))
        {
            *value = -INFINITY;   
            return(-1);
        }
        else
        {
            *value = output;
        }

        return 0;
    }

    void FC_R_Star::printDebugInfo(int startLoc, int startTime)
    {
        std::cout << "R_star debug info, location " << startLoc << ", time " << startTime << "\n";     
        std::cout << "...looking for problems\n";


        int i,j, compIdx;
        int nTpts = *((*R) -> nrow);
        int nLoc = *((*R) -> ncol);

        long double output = 0.0;
        
        double p_se_val;
        double p_rs_val;
        double ln_p_ir = std::log(**p_ir);
        double ln_1m_p_ir = std::log(1-(**p_ir));
        int Rstar_val, Sstar_val, Estar_val, R_val, I_val, S_val;   
        long unsigned int I_star_sum;
        long unsigned int R_star_sum;
        int64_t aDiff; 

        compIdx = startLoc*nTpts + startTime;
        if ((*context) -> config -> reinfectionMode <= 2)
        {
            for (j = startTime; j < nTpts; j++)
            {
                Rstar_val = ((*R_star) -> data)[compIdx];
                Sstar_val = ((*S_star)->data)[compIdx];
                Estar_val = ((*E_star) -> data)[compIdx];
                R_val = ((*R) ->data)[compIdx];
                I_val = ((*I) ->data)[compIdx];
                S_val = ((*S)->data)[compIdx];
                p_rs_val = (*p_rs)[j];

                if (Rstar_val < 0 || Rstar_val > I_val || 
                        Sstar_val > R_val)
                {

                    std::cout << "Bounds error detected, time " << j << "\n";
                    std::cout << "S_star: " << Sstar_val << "\n";
                    std::cout << "R_star: " << Rstar_val << "\n";
                    std::cout << "I: " << I_val << "\n";
                    std::cout << "R: " << R_val << "\n";
                    return;
                }
                else
                {
                    output += (ln_p_ir*Rstar_val +
                                ln_1m_p_ir*(I_val - Rstar_val) +
                                std::log(1-p_rs_val)*(R_val) +
                                ((*context) -> random -> choosePartial(R_val, Sstar_val)) +
                                ((*context) -> random -> choose(I_val, Rstar_val)));
                    if (!std::isfinite(output))
                    {
                        std::cout << "Calculation Error Detected, time " << j << "\n";
                        std::cout << "S_star: " << Sstar_val << "\n";
                        std::cout << "R_star: " << Rstar_val << "\n";
                        std::cout << "I: " << I_val << "\n";
                        std::cout << "R: " << R_val << "\n";
                        std::cout << "log(p_ir): " << ln_p_ir << "\n";
                        return;
                    }
                }
                compIdx++;
            } 
        }
        else 
        {
            for (j = startTime; j < nTpts; j++)
            {
                Rstar_val = ((*R_star) -> data)[compIdx];
                Estar_val = ((*E_star) -> data)[compIdx];
                R_val = ((*R) ->data)[compIdx];
                I_val = ((*I) ->data)[compIdx];
                S_val = ((*S)->data)[compIdx];

                if (Rstar_val < 0 || Rstar_val > I_val)
                {

                    std::cout << "Bounds error detected, time " << j << "\n";
                    std::cout << "R_star: " << Rstar_val << "\n";
                    std::cout << "I: " << I_val << "\n";
                    std::cout << "R: " << R_val << "\n";
                    return;
                }
                else
                {
                    output += (ln_p_ir*Rstar_val +
                                ln_1m_p_ir*(I_val - Rstar_val) +
                                ((*context) -> random -> choose(I_val, Rstar_val)));
                    if (!std::isfinite(output))
                    {
                        std::cout << "Calculation Error Detected, time " << j << "\n";
                        std::cout << "R_star: " << Rstar_val << "\n";
                        std::cout << "I: " << I_val << "\n";
                        std::cout << "R: " << R_val << "\n";
                        std::cout << "log(p_ir): " << ln_p_ir << "\n";
                        return;
                    }
                }
                compIdx++;
            } 

        }

        // p_se changes, so need to look at p_se component for all locations and 
        // time points after startTime
        for (i = 0; i < nLoc; i++)
        {
            compIdx = i*nTpts + startTime;
            for (j = startTime; j< nTpts; j++)
            {

                p_se_val = (*p_se)[compIdx];
                Estar_val = ((*E_star) -> data)[compIdx];
                S_val = ((*S)->data)[compIdx];
                if (p_se_val > 1 || p_se_val < 0)
                {
                    std::cout << "Non-Local Bounds Error Detected, loc, time: " << i << ", " << j << "\n";
                    std::cout << "E_star: " << Estar_val << "\n";
                    std::cout << "S: " << S_val << "\n";
                    std::cout << "p_se: " << p_se_val << "\n";
                    return;
                }

                output += (*context) -> random -> dbinom(Estar_val, S_val, p_se_val);

                if (!std::isfinite(output))
                {
                    std::cout<< "Non-Local Calculation Error Detected, loc, time: " << i << ", " << j << "\n";
                    std::cout << "E_star: " << Estar_val << "\n";
                    std::cout << "S: " << S_val << "\n";
                    std::cout << "p_se: " << p_se_val << "\n";
                    std::cout << "log(p_se): " << std::log(p_se_val) << "\n";
                    std::cout << "p_se == 0: " << (p_se_val == 0.0) << "\n";
                    std::cout << "Binomial term: " << ((*context) -> random -> 
                            dbinom(Estar_val, S_val, p_se_val)) << "\n";


                    return;
                }

                compIdx ++; 
            }
        }

        I_star_sum = (*I_star)->marginSum(2,startLoc);
        R_star_sum = (*R_star)->marginSum(2,startLoc);
        aDiff = (I_star_sum > R_star_sum ? I_star_sum - R_star_sum : R_star_sum - I_star_sum)/nTpts;
        output -= (aDiff*aDiff)*(*steadyStateConstraintPrecision);

        if (!std::isfinite(output))
        {
            std::cout << "Combinatorics Error Detected.\n";
            return;
        }

        std::cout << "No problems detected?\n";
        return;
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
        this -> sampleCompartment_CPU(*context,
                                  *R_star,*sliceWidth);
        return(0);
    }
    int FC_R_Star::sampleOCL()
    {
        //NOT IMPLEMENTED
        return -1;
    }

    long double FC_R_Star::getValue()
    {
        return(*(this -> value));
    }
    void FC_R_Star::setValue(long double val)
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
        value = new long double;

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
        int nLoc = *((*S) -> ncol);
        int nTpts = *((*S) -> nrow);
        double term1, term2, term3;
        term1 = 0.0; term2 = 0.0; term3 = 0.0;

        for (i = 0; i < nLoc; i++)    
        {
            compIdx = i*nTpts;
            for (j = 0; j < nTpts; j++)     
            {
                tmp = ((*E_star) -> data)[compIdx];
                term1 += std::log((*p_se)[compIdx])*tmp; 
                term2 += std::log(1-(*p_se)[compIdx])*(((*S) -> data)[compIdx] - tmp);
                compIdx++;
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

    long double FC_Beta::getValue()
    {
        return(*(this -> value));
    }
    void FC_Beta::setValue(long double val)
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
        value = new long double;

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
        int nTpts = *((*R) -> nrow);
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
            a = ((*S_star)-> marginSum(1,j));
            b = ((*R) -> marginSum(1,j)); 
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

    long double FC_Beta_P_RS::getValue()
    {
        return(*(this -> value));
    }
    void FC_Beta_P_RS::setValue(long double val)
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
        value = new long double;

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
        int nLoc = *((*S) -> ncol);
        int nTpts = *((*S) -> nrow);
        double term1, term2, term3;
        term1 = 0.0; term2 = 0.0; term3 = 0.0;

        for (i = 0; i < nLoc; i++)    
        {
            compIdx = i*nTpts;
            for (j = 0; j < nTpts; j++)     
            {
                Es = ((*E_star) -> data)[compIdx];
                pse = (*p_se)[compIdx];
                term1 += std::log(pse)*Es; 
                term2 += std::log(1-pse)*(((*S) -> data)[compIdx] - Es);
                compIdx++;
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

    long double FC_Rho::getValue()
    {
        return(*(this -> value));
    }
    void FC_Rho::setValue(long double val)
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
        value = new long double;

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
        int nLoc = *((*S) -> ncol);
        int nTpts = *((*S) -> nrow);
        double term1, term2, term3;
        term1 = 0.0; term2 = 0.0; term3 = 0.0;

        for (i = 0; i < nLoc; i++)    
        {
            compIdx = i*nTpts;
            for (j = 0; j < nTpts; j++)     
            {
                Es = ((*E_star) -> data)[compIdx];
                pse = (*p_se)[compIdx];
                term1 += std::log(pse)*Es; 
                term2 += std::log(1-pse)*(((*S) -> data)[compIdx] - Es);
                compIdx++;
            }
        } 
        for (j = 0; j < nTpts; j++)
        {
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

    long double FC_Gamma::getValue()
    {
        return(*(this -> value));
    }
    void FC_Gamma::setValue(long double val)
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
        value = new long double;

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
        int i_star_sum = (*I_star) -> marginSum(3,-1);
        int e_sum = (*E) -> marginSum(3,-1);
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

    long double FC_P_EI::getValue()
    {
        return(*(this -> value));
    }
    void FC_P_EI::setValue(long double val)
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
        value = new long double;

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
        int r_star_sum = (*R_star) -> marginSum(3,-1);
        int i_sum = (*I) -> marginSum(3,1);
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

    long double FC_P_IR::getValue()
    {
        return(*(this -> value));
    }
    void FC_P_IR::setValue(long double val)
    {
        *(this -> value) = val;
    }
}



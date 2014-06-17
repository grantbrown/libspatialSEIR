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

    int CompartmentFullConditional::sampleEntireCompartment_CPU(ModelContext* context,
                                                       CompartmentalModelMatrix* starCompartment,
                                                       double width)
    {
        ((*samples)) += 1;
        double initProposal = 0.0;
        double newProposal = 0.0;
        int i;
        int x0, x1;
        int totalPoints = (*(starCompartment -> nrow))*(*(starCompartment -> ncol));
        // Backup Compartment
        memcpy(context -> tmpContainer -> data, starCompartment -> data, totalPoints*sizeof(int)); 
        this -> calculateRelevantCompartments(); 
        this -> evalCPU();
        double initVal = (this -> getValue());
        if (! std::isfinite(initVal))
        {
            std::cerr << "Compartment sampler starting from value of zero probability!\n";
            throw(-1);
        }
        for (i = 0; i < totalPoints; i++)
        {
            x0 = (starCompartment -> data)[i];
            x1 = std::floor((context -> random -> normal(x0 + 0.5, width))); 
            (starCompartment -> data)[i] = x1;
            newProposal += (context -> random -> dnorm(x1, x0,width));
            initProposal += (context -> random -> dnorm(x0, x1,width));
        }
        this -> calculateRelevantCompartments(); 
        this -> evalCPU();
        double newVal = (this->getValue());
        double criterion = (newVal - initVal) + (initProposal - newProposal);

        if (std::log((context -> random -> uniform())) < criterion)
        {
            // Accept new values
            (*accepted) += 1;
        }
        else
        {
            // Keep Original Value
            memcpy(starCompartment -> data, context -> tmpContainer -> data, totalPoints*sizeof(int)); 
            this -> calculateRelevantCompartments();
            this -> setValue(initVal); 
        }                

        if (!std::isfinite(this -> getValue()))
        {
            std::cout << "Impossible value selected.\n";
            throw(-1);
        }

        return(0);
    }



    int CompartmentFullConditional::sampleCompartment_OCL(ModelContext* context,
                                                       CompartmentalModelMatrix* starCompartment,
                                                       double width)
    {
        (*samples) += 1;
        double initProposal = 0.0;
        double newProposal = 0.0;
        int i;
        int x0, x1;
        
        int totalPoints = (*(starCompartment -> nrow))*(*(starCompartment -> ncol));
        // Backup Compartment
        memcpy(context -> tmpContainer -> data, starCompartment -> data, totalPoints*sizeof(int)); 
        this -> calculateRelevantCompartments_OCL(); 
        this -> evalOCL();
        double initVal = (this -> getValue());
        if (! std::isfinite(initVal))
        {
            std::cerr << "Compartment sampler starting from value of zero probability!\n";
            throw(-1);
        }
        for (i = 0; i < totalPoints; i++)
        {
            x0 = (starCompartment -> data)[i];
            x1 = std::floor((context -> random -> normal(x0 + 0.5, width))); 
            (starCompartment -> data)[i] = x1;
            newProposal += (context -> random -> dnorm(x1, x0,width));
            initProposal += (context -> random -> dnorm(x0, x1,width));
        }
        this -> calculateRelevantCompartments_OCL(); 
        this -> evalOCL();
        double newVal = (this->getValue());
        double criterion = (newVal - initVal) + (initProposal - newProposal);

        if (std::log((context -> random -> uniform())) < criterion)
        {
            // Accept new values
            (*accepted) += 1;

        }
        else
        {
            // Keep Original Value
            memcpy(starCompartment -> data, context -> tmpContainer -> data, totalPoints*sizeof(int)); 
            this -> calculateRelevantCompartments_OCL();
            this -> setValue(initVal); 
        }                

        if (!std::isfinite(this -> getValue()))
        {
            std::cout << "Impossible value selected.\n";
            throw(-1);
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
   
            (*samples) += 1;
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
                (*accepted) += 1; 
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
            (*samples) += 1;
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
                (*accepted) += 1; 
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
                (*samples) += 1;
                x0 = ((context -> random -> uniform())*(r-l) + l);
                (starCompartment -> data)[compIdx] = std::floor(x0);
                this -> calculateRelevantCompartments(i,j); 
                this -> evalCPU(i,j);
                if (x0 >= x){r=x0;}
                else{l=x0;}
            } while (y >= (this -> getValue()));
            compIdx ++;
            (*accepted) += 1; 
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

    int InitCompartmentFullConditional::sampleEntireCompartment_CPU(ModelContext* context,
                                                                    int* initCompartment,
                                                                       double width)
    {
        (*samples) += 1;
        double initProposal = 0.0;
        double newProposal = 0.0;
        int i;
        int x0, x1;
        int nLoc = *(context -> S -> ncol);        
        // Backup Compartment
        memcpy(context -> tmpContainer -> data, initCompartment, nLoc*sizeof(int)); 
        this -> calculateRelevantCompartments(); 
        this -> evalCPU();
        double initVal = (this -> getValue());
        if (! std::isfinite(initVal))
        {
            std::cerr << "Compartment sampler starting from value of zero probability!\n";
            throw(-1);
        }
        for (i = 0; i < nLoc; i++)
        {
            x0 = (initCompartment)[i];
            x1 = std::floor((context -> random -> normal(x0 + 0.5, width))); 
            (initCompartment)[i] = x1;
            newProposal += (context -> random -> dnorm(x1, x0,width));
            initProposal += (context -> random -> dnorm(x0, x1,width));
        }
        this -> calculateRelevantCompartments(); 
        this -> evalCPU();
        double newVal = (this->getValue());
        double criterion = (newVal - initVal) + (initProposal - newProposal);

        if (std::log((context -> random -> uniform())) < criterion)
        {
            // Accept new values
            (*accepted) += 1;
        }
        else
        {
            // Keep Original Value
            memcpy(initCompartment, context -> tmpContainer -> data, nLoc*sizeof(int)); 
            this -> calculateRelevantCompartments();
            this -> setValue(initVal); 
        }                

        if (!std::isfinite(this -> getValue()))
        {
            std::cout << "Impossible value selected.\n";
            throw(-1);
        }

        return(0);
    }

    int InitCompartmentFullConditional::sampleEntireCompartment_OCL(ModelContext* context,
                                                                    int* initCompartment,
                                                                       double width)
    {
        (*samples) += 1;
        double initProposal = 0.0;
        double newProposal = 0.0;
        int i;
        int x0, x1;
        int nLoc = *(context -> S -> ncol);        
        // Backup Compartment
        memcpy(context -> tmpContainer -> data, initCompartment, nLoc*sizeof(int)); 
        this -> calculateRelevantCompartments_OCL(); 
        this -> evalOCL();
        double initVal = (this -> getValue());
        if (! std::isfinite(initVal))
        {
            std::cerr << "Compartment sampler starting from value of zero probability!\n";
            throw(-1);
        }
        for (i = 0; i < nLoc; i++)
        {
            x0 = (initCompartment)[i];
            x1 = std::floor((context -> random -> normal(x0 + 0.5, width))); 
            (initCompartment)[i] = x1;
            newProposal += (context -> random -> dnorm(x1, x0,width));
            initProposal += (context -> random -> dnorm(x0, x1,width));
        }
        this -> calculateRelevantCompartments_OCL(); 
        this -> evalOCL();
        double newVal = (this->getValue());
        double criterion = (newVal - initVal) + (initProposal - newProposal);

        if (std::log((context -> random -> uniform())) < criterion)
        {
            // Accept new values
            (*accepted) += 1;
        }
        else
        {
            // Keep Original Value
            memcpy(initCompartment, context -> tmpContainer -> data, nLoc*sizeof(int)); 
            this -> calculateRelevantCompartments();
            this -> setValue(initVal); 
        }                

        if (!std::isfinite(this -> getValue()))
        {
            std::cout << "Impossible value selected.\n";
            throw(-1);
        }

        return(0);
    }



    int InitCompartmentFullConditional::sampleCompartmentLocation(int i, ModelContext* context, 
                                                                  int* initCompartment,
                                                                  double width)
    {
        (*samples) += 1;
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
            (*accepted)+=1;
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
                (*samples) += 1;
                x0 = ((context -> random -> uniform())*(r-l) + l);
                variable[i] = x0;
                this -> calculateRelevantCompartments();
                this -> evalCPU();
                l = (x0 >= x ? l : x0);
                r = (x0 < x ? r : x0);  
            } while (y >= (this -> getValue()));
            (*accepted)++;
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
            (*samples) += 1;
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
                (*accepted)+=1;

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

    int ParameterFullConditional::sampleEntireDouble_CPU(ModelContext* context,
                                                         double* variable, 
                                                         int varLen, 
                                                         double width)
    {

        ((*samples)) += 1;
        double initProposal = 0.0;
        double newProposal = 0.0;
        int i;
        double x0, x1;

        // Backup Parameters (use the compartmentCache)
        memcpy((context -> compartmentCache), variable, varLen*sizeof(double)); 

        this -> calculateRelevantCompartments(); 
        this -> evalCPU();
        double initVal = (this -> getValue());
        if (! std::isfinite(initVal))
        {
            std::cerr << "Compartment sampler starting from value of zero probability!\n";
            throw(-1);
        }
        for (i = 0; i < varLen; i++)
        {

            x0 = (variable)[i];
            x1 = ((context -> random -> normal(x0, width))); 
            (variable)[i] = x1;
            newProposal += (context -> random -> dnorm(x1, x0,width));
            initProposal += (context -> random -> dnorm(x0, x1,width));
        }
        this -> calculateRelevantCompartments(); 
        this -> evalCPU();
        double newVal = (this->getValue());
        double criterion = (newVal - initVal) + (initProposal - newProposal);


        if (std::log((context -> random -> uniform())) < criterion)
        {
            // Accept new values
            (*accepted) += 1;
        }
        else
        {
            // Keep Original Value
            memcpy(variable, (context -> compartmentCache), varLen*sizeof(double)); 
            this -> calculateRelevantCompartments();
            this -> setValue(initVal); 
        }                

        if (!std::isfinite(this -> getValue()))
        {
            std::cout << "Impossible value selected.\n";
            throw(-1);
        }
        return(0);
    }

    int ParameterFullConditional::sampleEntireDouble_OCL(ModelContext* context,
                                                         double* variable, 
                                                         int varLen, 
                                                         double width)
    {

        ((*samples)) += 1;
        double initProposal = 0.0;
        double newProposal = 0.0;
        int i;
        double x0, x1;
        // Backup Parameters (use the compartmentCache cache)
        memcpy(context -> compartmentCache, variable, varLen*sizeof(double)); 
        this -> calculateRelevantCompartments_OCL(); 
        this -> evalOCL();
        double initVal = (this -> getValue());
        if (! std::isfinite(initVal))
        {
            std::cerr << "Compartment sampler starting from value of zero probability!\n";
            throw(-1);
        }
        for (i = 0; i < varLen; i++)
        {
            x0 = (variable)[i];
            x1 = ((context -> random -> normal(x0, width))); 
            (variable)[i] = x1;
            newProposal += (context -> random -> dnorm(x1, x0,width));
            initProposal += (context -> random -> dnorm(x0, x1,width));
        }
        this -> calculateRelevantCompartments_OCL(); 
        this -> evalOCL();
        double newVal = (this->getValue());
        double criterion = (newVal - initVal) + (initProposal - newProposal);

        if (std::log((context -> random -> uniform())) < criterion)
        {
            // Accept new values
            (*accepted) += 1;
        }
        else
        {
            // Keep Original Value
            memcpy(variable, context -> compartmentCache, varLen*sizeof(double)); 
            this -> calculateRelevantCompartments_OCL();
            this -> setValue(initVal); 
        }                

        if (!std::isfinite(this -> getValue()))
        {
            std::cout << "Impossible value selected.\n";
            throw(-1);
        }

        return(0);
    }

    int ParameterFullConditional::sampleDouble_OCL(ModelContext* context,
                                                         double* variable, 
                                                         int varLen, 
                                                         double width)
    {
        // Declare required variables
        int i;
        double x0,x1;
        double initVal, newVal, initProposal, newProposal;

        // Update the relevant CompartmentalModelMatrix instances
        this -> calculateRelevantCompartments_OCL();

        // Set the "value" attribute appropriately
        this -> evalOCL();
   
        // Main loop: 
        for (i = 0; i < varLen; i++)
        { 
            (*samples) += 1;
            x0 = variable[i];
            this -> calculateRelevantCompartments_OCL(); 
            this -> evalOCL();
            initVal = (this->getValue());

            x1 = (context -> random -> normal(x0, width));
            variable[i] = x1;
            this -> calculateRelevantCompartments_OCL();
            this -> evalOCL();
            newVal = (this->getValue());
            initProposal = (context->random->dnorm(x0, x1, width));
            newProposal = (context->random->dnorm(x1, x0, width)); 

            if (std::log((context -> random -> uniform())) < ((newVal - initVal) + (initProposal - newProposal)))
            {
                // Accept the new value. 
                (*accepted)+=1;

            }
            else
            {
                // Keep original value
                variable[i] = x0;
                this -> calculateRelevantCompartments_OCL();
                this -> setValue(initVal);
            }
        }
        return 0;
    }
}

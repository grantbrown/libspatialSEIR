#include<math.h>
#include<cstring>
#include<vector>
#include<cmath>
#include<algorithm>
#include<LSS_Samplers.hpp>
#include<LSS_FullConditional.hpp>
#include<ModelContext.hpp>
#include<OCLProvider.hpp>
#include<CompartmentalModelMatrix.hpp>
#include<CovariateMatrix.hpp>
#include<RandomNumberProvider.hpp>
#include<IOProvider.hpp>

namespace SpatialSEIR
{

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
        delete[] S0;
        delete[] E0;
        delete[] I0;
        delete[] R0;
        delete numLocations;
    }

    double FullConditional::acceptanceRatio()
    {
        return((*accepted*1.0)/(*samples));
    }

    void FullConditional::setSamplerType(int samplerType)
    {
        unsigned int i;
        for (i = 0; i < (this -> samplers) -> size(); i++)
        {
            if (((*samplers)[i] -> getSamplerType()) == samplerType)
            {
                *currentSampler = (*samplers)[i];
                return;
            }
        }
        lssCout << "Sampler type not found.\n";
        throw(-1);
    }


    double ParameterFullConditional::acceptanceRatio(int i)
    {
        return((accepted[i]*1.0)/(*samples));
    }

    int CompartmentFullConditional::getFullConditionalType()
    {
        return(LSS_COMPARTMENT_FULL_CONDITIONAL_TYPE);
    }

    int InitCompartmentFullConditional::getFullConditionalType()
    {
        return(LSS_INIT_COMPARTMENT_FULL_CONDITIONAL_TYPE);
    }

    int ParameterFullConditional::getFullConditionalType()
    {
        return(LSS_PARAMETER_FULL_CONDITIONAL_TYPE);
    }


    /*
     *
     * Auto-tuning logic
     *
     */

    void ParameterFullConditional::updateSamplingParameters(double desiredRatio, double targetWidth, double proportionChange)
    {
        if (*samples == 0)
        {
            // Do nothing
            return;
        }
        if (proportionChange <= 0 || proportionChange >= 1)
        {
            lssCout << "Invalid Proportion: " << proportionChange << "\n";
            throw(-1);
        }
        int i;
        for (i = 0; i < *varLen; i++)
        {
            double currentRatio = (this -> acceptanceRatio(i));
            if ((currentRatio > desiredRatio) && (std::abs(currentRatio - desiredRatio) > targetWidth))
            {
                (sliceWidth[i])*=(1+proportionChange); 
            }
            else if ((currentRatio < desiredRatio) && (std::abs(currentRatio - desiredRatio) > targetWidth))
            {
                (sliceWidth[i])*=(1-proportionChange);           
            }
            accepted[i] = 0;
        }
        
        *samples = 0;
    }

    void CompartmentFullConditional::updateSamplingParameters(double desiredRatio, double targetWidth, double proportionChange)
    {
        if (*samples == 0)
        {
            // Do nothing
            return;
        }
        if (proportionChange <= 0 || proportionChange >= 1)
        {
            lssCout << "Invalid Proportion: " << proportionChange << "\n";
            throw(-1);
        }
        double currentRatio = (this -> acceptanceRatio());
        if ((currentRatio > desiredRatio) && (std::abs(currentRatio - desiredRatio) > targetWidth))
        {
            (*sliceWidth)*=(1+proportionChange); 
        }
        else if ((currentRatio < desiredRatio) && (std::abs(currentRatio - desiredRatio) > targetWidth))
        {
            (*sliceWidth)*=(1-proportionChange);           
        }

        *samples = 0;
        *accepted = 0;
    }

    void InitCompartmentFullConditional::updateSamplingParameters(double desiredRatio, double targetWidth, double proportionChange)
    {
        if (*samples == 0)
        {
            // Do nothing
            return;
        }
        if (proportionChange <= 0 || proportionChange >= 1)
        {
            lssCout << "Invalid Proportion: " << proportionChange << "\n";
            throw(-1);
        }
        double currentRatio = (this -> acceptanceRatio());
        if ((currentRatio > desiredRatio) && (std::abs(currentRatio - desiredRatio) > targetWidth))
        {
            (*sliceWidth)*=(1+proportionChange); 
        }
        else if ((currentRatio < desiredRatio) && (std::abs(currentRatio - desiredRatio) > targetWidth))
        {
            (*sliceWidth)*=(1-proportionChange);           
        }

        *samples = 0;
        *accepted = 0;
    }
}




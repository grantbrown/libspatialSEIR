#include <ModelContext.hpp>
#include <FullConditional.hpp>
#include <OCLProvider.hpp>
#include <CompartmentalModelMatrix.hpp>
#include <CovariateMatrix.hpp>
#include <DistanceMatrix.hpp>
#include <RandomNumberProvider.hpp>
#include <IOProvider.hpp>


#ifndef BLAS_INC
#define BLAS_INC
#include<cblas.h> 
#endif

#include<cmath>
#include<ctime>

namespace SpatialSEIR
{
    using std::cout;
    using std::endl;
    int ModelContext::totalS()
    {
        return(this->S->marginSum(3,-1));
    }
    int ModelContext::totalE()
    {
        return(this->E->marginSum(3,-1));
    }
    int ModelContext::totalI()
    {
        return(this->I->marginSum(3,-1));
    }
    int ModelContext::totalR()
    {
        return(this->R->marginSum(3,-1));
    }

    int ModelContext::totalS(int tpt)
    {
        return(this->S->marginSum(2,tpt));
    }
    int ModelContext::totalE(int tpt)
    {
        return(this->E->marginSum(2,tpt));
    }
    int ModelContext::totalI(int tpt)
    {
        return(this->I->marginSum(2,tpt));
    }
    int ModelContext::totalR(int tpt)
    {
        return(this->R->marginSum(2,tpt));
    }

    double ModelContext::avgP_SE()
    {
        double out = 0.0;
        int numVals = (*(S->nrow))*(*(S->ncol));
        int i;
        for (i = 0; i < numVals; i ++)
        {
            out += p_se[i];
        }
        return(out/numVals);
    }

    double ModelContext::avgP_SE(int tpt)
    {
        double out = 0.0;
        int numLoc = *(S->nrow);
        int startVal = numLoc*tpt; 
        int i;
        for (i = 0; i < numLoc; i ++)
        {
            out += p_se[startVal + i];
        }
        return(out/numLoc);
    }

    double ModelContext::avgP_RS()
    {
        int numTpt = *(S->nrow);
        int i;
        double out = 0.0;
        for (i = 0; i< numTpt; i++)
        {
            out += p_rs[i];
        }
        return(out/numTpt);
    }






}



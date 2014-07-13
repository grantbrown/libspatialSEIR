#include <ModelContext.hpp>
#include <LSS_FullConditionalList.hpp>
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

    int ModelContext::totalS_star()
    {
        return(this->S_star->marginSum(3,-1));
    }
    int ModelContext::totalE_star()
    {
        return(this->E_star->marginSum(3,-1));
    }
    int ModelContext::totalI_star()
    {
        return(this->I_star->marginSum(3,-1));
    }
    int ModelContext::totalR_star()
    {
        return(this->R_star->marginSum(3,-1));
    }

    int ModelContext::totalS_star(int tpt)
    {
        return(this->S_star->marginSum(2,tpt));
    }
    int ModelContext::totalE_star(int tpt)
    {
        return(this->E_star->marginSum(2,tpt));
    }
    int ModelContext::totalI_star(int tpt)
    {
        return(this->I_star->marginSum(2,tpt));
    }
    int ModelContext::totalR_star(int tpt)
    {
        return(this->R_star->marginSum(2,tpt));
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

    double ModelContext::estimateR0()
    {
        return(-1.0);
    }

    double ModelContext::estimateR0(int j)
    {
        return(-1.0);
    }

    double* ModelContext::calculateG(int j)
    {
        int i, l;

        //Update Eta
        this -> X -> calculate_eta_CPU(eta, beta);

        int iIndex, lIndex, GIndex;
        int nLoc = *(S -> ncol);
        int nTpt = *(S -> nrow);
        double* G = new double[nLoc*nLoc];
        //Exponentiate
        int nrowz = *(X->nrow_z);
        for (i = 0; i < nrowz; i++)
        {
            eta[i] = std::exp(eta[i]);
        }
        for (i = 0; i < nLoc; i++) 
        {
            for (l = 0; l < nLoc; l++)
            { 
                iIndex = i*nTpt + j;
                lIndex = l*nTpt + j;

                GIndex = l*nLoc + i;
                if (i != l)
                {
                    G[GIndex] = ((I->data)[lIndex] != 0 ?
                                    (-(((N[iIndex])/((I->data)[lIndex]))
                                    * (1-std::exp((*rho)*((scaledDistMat->data)[GIndex]) 
                                    * (((I -> data)[lIndex] * (eta[lIndex]))/N[lIndex]))))) :
                                        0.0 );
                }
                else
                { 
                    G[GIndex] = ((I->data)[lIndex] != 0 ?
                                    -(((N[iIndex])/((I->data)[lIndex]))
                                    * (1-std::exp((((I -> data)[lIndex] * (eta[lIndex]))/N[lIndex])))) : 
                                     0.0 );
                }
            }
        }
        return(G);
    }
}

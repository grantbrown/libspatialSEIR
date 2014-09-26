#include <ModelContext.hpp>
#include <LSS_FullConditionalList.hpp>
#include <OCLProvider.hpp>
#include <CompartmentalModelMatrix.hpp>
#include <CovariateMatrix.hpp>
#include <DistanceMatrix.hpp>
#include <RandomNumberProvider.hpp>
#include <IOProvider.hpp>
#include <Eigen/Core>
#include <Eigen/Eigenvalues>


#ifndef BLAS_INC
#define BLAS_INC
#include<cblas.h> 
#endif

#include<cmath>
#include<ctime>

namespace SpatialSEIR
{
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
        int numTpt = *(S->nrow);
        double output = 0.0;
        int i;
        for (i = 0; i < numTpt; i++)
        {
            output += estimateR0(i);
        }
        return(output/numTpt);
    }

    double ModelContext::estimateR0(int j)
    {
        double* G; 
        int nLoc = *(S -> ncol);
        G = calculateIntegratedG(j);
        MatrixMapType Gmap(G, nLoc, nLoc);
        Eigen::EigenSolver<Eigen::MatrixXd> Es(Gmap);
        double output =  Es.eigenvalues()[0].real();
        delete[] G;
        return(output);
    }

    double ModelContext::estimateR0(int i, int j)
    {
        // Column sum
        double* G = calculateIntegratedG(j); 
        double out = 0.0;
        int nLoc = *(S -> ncol);
        int k;

        for (k = 0; k < nLoc; k++)
        {
            out += G[i*nLoc + k];
        }
        delete[] G;
        return(out); 
    }

    double* ModelContext::calculateIntegratedG(int j)
    {

        int i,k;
        int p_ir_idx = j;
        int nLoc = *(S -> ncol);
        int nTpt = *(S -> nrow);
        if (j >= nTpt)
        {
            lssCout << "Invalid time point: " << j << "\n";
            throw(-1);
        }
        calculateP_IR_CPU(); 
        double* outG = calculateG(j);
        double* newG; 
        double pIR;
        double _1mpIR_cum = 1.0;
        // Need to do prediction here for later time points?
        for (i = j+1; ((i < nTpt) && (_1mpIR_cum >= 1e-8)); i++)
        {
            newG = calculateG(i);
            pIR = (p_ir)[p_ir_idx];
            _1mpIR_cum *= (1-pIR);
            for (k = 0; k < nLoc*nLoc; k++)
            {
                outG[k] += _1mpIR_cum*newG[k];
            }
            p_ir_idx++;
            delete[] newG;
        }
        // If we're truncated at end of time period
        // do a simple extrapolation
        i = nTpt - 1;
        newG = calculateG(i); 
        pIR = (p_ir)[p_ir_idx];
        while (_1mpIR_cum >= 1e-8)
        {
            _1mpIR_cum *= (1-pIR);
            for (k = 0; k < nLoc*nLoc; k++)
            {
                outG[k] += _1mpIR_cum*newG[k];
            }
        }
        delete[] newG;          

        return(outG);
    }

    double* ModelContext::calculateG(int j)
    {
        int i, l;

        //Update Eta
        this -> X -> calculate_eta_CPU(eta, beta);

        int iIndex, lIndex, GIndex;
        int nLoc = *(S -> ncol);
        int nTpt = *(S -> nrow);
        if (j >= nTpt)
        {
            lssCout << "Invalid time point: " << j << "\n";
            throw(-1);
        }

        double* G = new double[nLoc*nLoc];
        double component1, component2;
        unsigned int k;
        //Exponentiate
        int nrowz = *(X->nrow_z);
        for (i = 0; i < nrowz; i++)
        {
            eta[i] = std::exp(eta[i]);
        }
        // Out: rows
        for (i = 0; i < nLoc; i++) 
        {

            // Out: columns
            for (l = 0; l < nLoc; l++)
            { 
                iIndex = i*nTpt + j;
                lIndex = l*nTpt + j;
                GIndex = l*nLoc + i;
                if (i != l)
                {
                    component1 = (((I -> data)[lIndex] * (eta[lIndex]))/N[lIndex]);
                    component2 = 0.0;
                    for (k = 0; k < scaledDistMatrices -> size(); k++)
                    {
                        component2 += ((rho)[k])*(((*scaledDistMatrices)[k] -> data)[GIndex])*component1;
                    }
                    G[GIndex] = ((I->data)[lIndex] != 0 ?
                                    ((((N[lIndex])/((I -> data)[lIndex]))
                                    * (1-std::exp(-(offset[j])*component2)))) :
                                        0.0 );
                }
                else
                { 
                    G[GIndex] = ((I->data)[lIndex] != 0 ?
                                    (((N[iIndex])/((I -> data)[iIndex]))
                                    * (1-std::exp(-(offset[j])*(((I -> data)[lIndex] * (eta[lIndex]))/N[lIndex])))) : 
                                     0.0 );
                }
            }
        }
        return(G);
    }
}

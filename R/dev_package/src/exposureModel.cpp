#include <Rcpp.h>
#include <exposureModel.hpp>
#include <DistanceMatrix.hpp>


using namespace Rcpp;
using namespace SpatialSEIR;

exposureModel::exposureModel(SEXP _X, SEXP _ntpt, SEXP _nloc, SEXP _paramInit, SEXP _priorMean, SEXP _prec)
{
    int i;
    Rcpp::NumericMatrix inX(_X);
    Rcpp::NumericVector inTpt(_ntpt);
    Rcpp::NumericVector inLoc(_nloc);
    Rcpp::NumericVector inPrecision(_prec);
    Rcpp::NumericVector priorMeans(_priorMean);
    Rcpp::NumericVector initParams(_paramInit);

    xDim = new int[2];
    nTpt = new int;
    nLoc = new int;
    *nTpt = inTpt[0];
    *nLoc = inLoc[0];
    xDim[0] = inX.nrow();
    xDim[1] = inX.ncol();

    X = new double[xDim[0]*xDim[1]];
    int nBeta = xDim[1];
    beta = new double[nBeta];
    betaPriorPrecision = new double[nBeta];
    betaPriorMean = new double[nBeta];

    if (xDim[0] != ((*nTpt) *(*nLoc)))
    {
        Rcpp::Rcout << "Error: Covariate matrix has invalid dimensions." 
                    << " X should be (n*t) by p\n";
        ::Rf_error("error building exposure model");
    }
    if (initParams.length() != nBeta || inPrecision.length() != nBeta || priorMeans.length() != nBeta)
    {
        Rcpp::Rcout << "Initial parameters have different length then number of columns of X\n";
        Rcpp::Rcout << "ncol(X): " << xDim[1] << "\n";
        Rcpp::Rcout << "length(beta): " << initParams.length() << "\n";
        Rcpp::Rcout << "length(betaPriorMean): " << priorMeans.length() << "\n";
        Rcpp::Rcout << "length(betaPriorPrecision): " << inPrecision.length() << "\n";
        ::Rf_error("error building exposure model");
    }
    offset = new double[*nTpt];
    for (i = 0; i < *nTpt; i++)
    {
        offset[i] = 1.0;
    }
    for (i = 0; i < xDim[0]*xDim[1]; i++)
    {
        X[i] = inX[i];
    }
    for (i = 0; i < xDim[1]; i++)
    {
        beta[i] = initParams[i]; 
        betaPriorMean[i] = priorMeans[i];
        betaPriorPrecision[i] = inPrecision[i];
    }
}

int exposureModel::getModelComponentType()
{
    return(LSS_EXPOSURE_MODEL_TYPE);
}


void exposureModel::setOffset(NumericVector offsets)
{
    if (offsets.length() != (*nTpt))
    {
        Rcpp::Rcout << "Error: offsets must have length equal to the number of time points.\n";
        ::Rf_error("Invalid offsets.");
    }
    memcpy(offset, offsets.begin(), (*nTpt)*sizeof(double));
}

Rcpp::NumericVector exposureModel::getOffset()
{
    Rcpp::NumericVector output((*nTpt));
    memcpy(output.begin(), offset, (*nTpt)*sizeof(double));
    return(output);
}

void exposureModel::summary()
{
    Rcpp::Rcout << "Covariate matrix dimensions: " << xDim[0] << ", " << xDim[1] << "\n";
    Rcpp::Rcout << "Number of time points: " << *nTpt << "\n";
    Rcpp::Rcout << "Number of locations: " << *nLoc << "\n";
    Rcpp::Rcout << "Linear model parameter prior precision: " << "\n";
    int i;
    for (i = 0; i < xDim[1]; i++)
    {
        Rcpp::Rcout << betaPriorPrecision[i] << ", ";
    }
    Rcpp::Rcout << "\n";

    Rcpp::Rcout << "Linear model parameter prior mean: " << "\n";
    for (i = 0; i < xDim[1]; i++)
    {
        Rcpp::Rcout << betaPriorMean[i] << ", ";
    }
    Rcpp::Rcout << "\n";
}

exposureModel::~exposureModel()
{
    delete[] xDim;
    delete[] X;
    delete nLoc;
    delete nTpt;
    delete[] beta;
    delete[] betaPriorPrecision;
    delete[] betaPriorMean;
}

RCPP_MODULE(mod_exposureModel)
{
    using namespace Rcpp;
    class_<exposureModel>( "exposureModel" )
    .constructor<SEXP,SEXP,SEXP,SEXP,SEXP,SEXP>()
    .method("summary", &exposureModel::summary)
    .property("offsets", &exposureModel::getOffset, &exposureModel::setOffset);
}



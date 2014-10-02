#include <Rcpp.h>
#include <exposureModel.hpp>
#include <DistanceMatrix.hpp>


using namespace Rcpp;
using namespace SpatialSEIR;

exposureModel::exposureModel(SEXP _X, SEXP _Z, SEXP _paramInit, SEXP _prec, SEXP _hasZ)
{
    int i;
    Rcpp::NumericMatrix inX(_X);
    Rcpp::NumericMatrix inZ(_Z); 
    Rcpp::NumericVector inPrecision(_prec);
    Rcpp::NumericVector initParams(_paramInit);
    Rcpp::IntegerVector hasZVec(_hasZ);
    int hasZ = hasZVec[0];

    betaPriorPrecision = new double;
    *betaPriorPrecision = *(inPrecision.begin()); 
    xDim = new int[2];
    zDim = new int[2];
    xDim[0] = inX.nrow();
    xDim[1] = inX.ncol();
    if (hasZ)
    {
        zDim[0] = inZ.nrow();
        zDim[1] = inZ.ncol();
        Z = new double[zDim[0]*zDim[1]];
    }
    else
    {
        zDim[0] = inZ[0]*xDim[0];
        zDim[1] = 0;
        Z = new double[zDim[0]];
        for (i = 0; i < zDim[0]; i++){Z[i] = 0.0;} 
    }
    X = new double[xDim[0]*xDim[1]];

    beta = new double[xDim[1] + zDim[1]];
    if (zDim[0] % xDim[0] != 0)
    {
        Rcpp::Rcout << "Error: Covariate matrices have invalid dimensions." 
                    << " X should have one row per spatial unit, and"
                    << " Z should have nrow(X) blocks of length equal"
                    <<" to the number of time points.\n";
        throw(-1);
    }
    if (initParams.length() != xDim[1] + zDim[1])
    {
        Rcpp::Rcout << "Initial parameters have different length then number of columns of X and Z\n";
        Rcpp::Rcout << "ncol(X): " << xDim[1] << "\n";
        Rcpp::Rcout << "ncol(Z): " << zDim[1] << "\n";
        Rcpp::Rcout << "length(beta): " << initParams.length() << "\n";
        Rcpp::Rcout << initParams << "\n";
        throw(-1);
    }
    offset = new double[(zDim[0])/(xDim[0])];
    for (i = 0; i < (zDim[0])/(xDim[0]); i++)
    {
        offset[i] = 1.0;
    }
    for (i = 0; i < xDim[0]*xDim[1]; i++)
    {
        X[i] = inX[i];
    }
    for (i = 0; i < zDim[0]*zDim[1]; i++)
    {
        Z[i] = inZ[i];
    }
    for (i = 0; i < (xDim[1]+zDim[1]); i++)
    {
        beta[i] = initParams[i]; 
    }
}


void exposureModel::setOffset(NumericVector offsets)
{
    if (offsets.length() != ((zDim[0])/(xDim[0])))
    {
        Rcpp::Rcout << "Error: offsets must have length equal to the number of time points.\n";
        throw(-1);
    }
    memcpy(offset, offsets.begin(), ((zDim[0])/(xDim[0]))*sizeof(double));
}

Rcpp::NumericVector exposureModel::getOffset()
{
    Rcpp::NumericVector output(((zDim[0])/(xDim[0])));
    memcpy(output.begin(), offset, ((zDim[0])/(xDim[0]))*sizeof(double));
    return(output);
}

void exposureModel::summary()
{
    Rcpp::Rcout << "Time invariant covariate dimensions: " << xDim[0] << ", " << xDim[1] << "\n";
    Rcpp::Rcout << "Time varying covariate dimensions: " << zDim[0] << ", " << zDim[1] << "\n";
    Rcpp::Rcout << "Linear model parameter prior precision: " << *betaPriorPrecision << "\n";

}

exposureModel::~exposureModel()
{
    delete[] xDim;
    delete[] zDim;
    delete[] X;
    delete[] Z;
    delete[] beta;
    delete betaPriorPrecision;
}

RCPP_MODULE(mod_exposureModel)
{
    using namespace Rcpp;
    class_<exposureModel>( "exposureModel" )
    .constructor<SEXP,SEXP,SEXP,SEXP,SEXP>()
    .method("summary", &exposureModel::summary)
    .property("offsets", &exposureModel::getOffset, &exposureModel::setOffset);
}



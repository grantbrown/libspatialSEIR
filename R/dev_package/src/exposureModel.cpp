#include <Rcpp.h>
#include <exposureModel.hpp>
#include <DistanceMatrix.hpp>


using namespace Rcpp;
using namespace SpatialSEIR;

exposureModel::exposureModel(SEXP _X, SEXP _Z, SEXP _prec)
{
    Rcpp::NumericMatrix inX(_X);
    Rcpp::NumericMatrix inZ(_Z); 
    Rcpp::NumericVector inPrecision(_prec);

    betaPriorPrecision = new double;
    *betaPriorPrecision = *(inPrecision.begin()); 
    xDim = new int[2];
    zDim = new int[2];
    xDim[0] = inX.nrow();
    xDim[1] = inX.ncol();
    zDim[0] = inZ.nrow();
    zDim[1] = inZ.ncol();
    X = new double[xDim[0]*xDim[1]];
    Z = new double[zDim[0]*zDim[1]];
    if (zDim[0] % xDim[0] != 0)
    {
        Rcpp::Rcout << "Error: Covariate matrices have invalid dimensions." 
                    << " X should have one row per spatial unit, and"
                    << " Z should have nrow(X) blocks of length equal"
                    <<" to the number of time points.\n";
        throw(-1);
    }
    offset = new int[(zDim[0])/(xDim[0])];
    memset(offset, 1, (zDim[0])/(xDim[0])*sizeof(int));
    memcpy(X, inX.begin(), xDim[0]*xDim[1]*sizeof(double));
    memcpy(Z, inZ.begin(), zDim[0]*zDim[1]*sizeof(double));
}


void exposureModel::setOffset(SEXP _offs)
{
    Rcpp::IntegerVector offsets(_offs);
    if (offsets.length() != ((zDim[0])/(xDim[0])))
    {
        Rcpp::Rcout << "Error: offsets must have length equal to the number of time points.\n";
        throw(-1);
    }
    memcpy(offset, offsets.begin(), ((zDim[0])/(xDim[0]))*sizeof(int));
}

Rcpp::IntegerVector exposureModel::getOffset()
{
    Rcpp::IntegerVector output(((zDim[0])/(xDim[0])));
    memcpy(output.begin(), offset, ((zDim[0])/(xDim[0]))*sizeof(int));
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
    delete betaPriorPrecision;
}

RCPP_MODULE(mod_exposureModel)
{
    using namespace Rcpp;
    class_<exposureModel>( "exposureModel" )
    .constructor<SEXP,SEXP,SEXP>()
    .method("summary", &exposureModel::summary)
    .property("offsets", &exposureModel::getOffset, &exposureModel::setOffset, "Vector of temporal offsets.");
}



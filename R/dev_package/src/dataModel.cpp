#include <Rcpp.h>
#include <dataModel.hpp>
#include <DistanceMatrix.hpp>


using namespace Rcpp;
using namespace SpatialSEIR;

dataModel::dataModel(SEXP _Y, SEXP type)
{
    Rcpp::NumericMatrix input(_Y);
    nLoc = new int; *nLoc = input.ncol();
    nTpt = new int; *nTpt = input.nrow();
    dataModelType = new int;
    Rcpp::StringVector inputType(type);
    if (inputType[0] == "identity")
    {
        *dataModelType = 0;
    }
    else if (inputType[0] == "overdispersion")
    {
        *dataModelType = 1;
    }
    else
    {
        Rcpp::Rcout << "Unrecognized data model type: " << type << "\n";
        Rcpp::Rcout << "Falling back to default: identity\n";
        *dataModelType = 0;
    }
    Y = new int[(*nLoc)*(*nTpt)];
    memcpy(Y, input.begin(), (*nLoc)*(*nTpt)*sizeof(int));
}

void dataModel::summary()
{
    Rcpp::Rcout << "Number of locations: " << *nLoc << "\n";
    Rcpp::Rcout << "Number of time points: " << *nLoc << "\n";
    Rcpp::Rcout << "Data Model Type: ";
    if (*dataModelType == 0)
    {
        Rcpp::Rcout << "identity\n";
    }
    else if (*dataModelType == 1)
    {
        Rcpp::Rcout << "overdispersion\n";
    }
    else
    {
        Rcpp::Rcout << "Unknown\n";
    }
}
dataModel::~dataModel()
{
    delete[] Y;
    delete nLoc;
    delete nTpt;
}

RCPP_MODULE(mod_dataModel)
{
    using namespace Rcpp;
    class_<dataModel>( "dataModel" )
    .constructor<SEXP,SEXP>()
    .method("summary", &dataModel::summary);
}



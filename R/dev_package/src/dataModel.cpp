#include <Rcpp.h>
#include <dataModel.hpp>
#include <DistanceMatrix.hpp>


using namespace Rcpp;
using namespace SpatialSEIR;

dataModel::dataModel(SEXP _Y, SEXP type, SEXP compartment)
{
    setMode = new int; *setMode = -1;
    // We don't actually need 10 parameters, but leave room for future data models. 
    priorParameters = new double[10]; memset(priorParameters, -1.0, 10*sizeof(double));
    initialParameterValues = new double[10]; memset(initialParameterValues, -1.0, 10*sizeof(double));
    Rcpp::NumericMatrix input(_Y);
    compartmentDimensions = new Rcpp::IntegerVector(2);
    nLoc = new int; *nLoc = input.ncol();
    nTpt = new int; *nTpt = input.nrow();
    (*compartmentDimensions)[0] = *nTpt;
    (*compartmentDimensions)[1] = *nLoc;
    dataModelType = new int;
    dataModelCompartment = new int;
    Rcpp::StringVector inputType(type);
    Rcpp::StringVector inputCompartment(compartment);
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

    if (inputCompartment[0] == "I_star")
    {
        *dataModelCompartment = 0;
    }
    else if (inputCompartment[0] == "R_star")
    {
        *dataModelCompartment = 1;
    }
    else
    {
        Rcpp::Rcout << "Unrecognized data model compartment: " << compartment << "\n";
        Rcpp::Rcout << "Falling back to default: I_star\n";
        *dataModelCompartment = 0;
    }
    Y = new int[(*nLoc)*(*nTpt)];
    int i;
    for (i = 0; i < (*nLoc)*(*nTpt); i++)
    {
        Y[i] = input[i]; 
    }
}

void dataModel::setOverdispersionParameters(SEXP priorAlpha, SEXP priorBeta, SEXP initialValue)
{
    *setMode = 1;
    Rcpp::NumericVector alpha(priorAlpha);
    Rcpp::NumericVector beta(priorBeta);
    Rcpp::NumericVector init(initialValue);
    
    priorParameters[0] = alpha[0];
    priorParameters[1] = beta[0];
    initialParameterValues[0] = init[0];
}

int dataModel::getModelComponentType()
{
    return(LSS_DATA_MODEL_TYPE);
}

void dataModel::summary()
{
    Rcpp::Rcout << "Number of locations: " << *nLoc << "\n";
    Rcpp::Rcout << "Number of time points: " << *nLoc << "\n";
    Rcpp::Rcout << "Data Model Compartment: ";
    if (*dataModelCompartment == 0)
    {
        Rcpp::Rcout << "I_star\n";
    }
    else if (*dataModelCompartment == 1)
    {
        Rcpp::Rcout << "R_star\n";
    }
    else
    {
        Rcpp::Rcout << "Unknown\n";
    }
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
    delete[] priorParameters;
    delete[] initialParameterValues;
    delete setMode;
    delete compartmentDimensions;
    delete nLoc;
    delete nTpt;
}

RCPP_MODULE(mod_dataModel)
{
    using namespace Rcpp;
    class_<dataModel>( "dataModel" )
    .constructor<SEXP,SEXP,SEXP>()
    .method("summary", &dataModel::summary)
    .method("setOverdispersionParameters",&dataModel::setOverdispersionParameters);
}



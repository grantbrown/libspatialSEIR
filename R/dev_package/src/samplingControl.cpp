#include <Rcpp.h>
#include <samplingControl.hpp>
#include <DistanceMatrix.hpp>


using namespace Rcpp;
using namespace SpatialSEIR;

samplingControl::samplingControl()
{
    iterationStride = new int; *iterationStride = 100;
    steadyStateConstraintPrecision = new double; *steadyStateConstraintPrecision = -1.0;
    verbose = new int; *verbose = 0;
    debug = new int; *debug = 0;
    sliceWidths = new double[11];
    memset(sliceWidths, 0.1, 11*sizeof(double));
}

void samplingControl::setIterationStride(SEXP stride)
{
    Rcpp::IntegerVector input(stride);
    *iterationStride = input[0];
}

Rcpp::IntegerVector samplingControl::getIterationStride()
{
    Rcpp::IntegerVector output(1);
    output[0] = *iterationStride;
    return(output);
}

void samplingControl::setVerbose(SEXP vb)
{
    Rcpp::IntegerVector input(vb);
    *verbose = input[0];
}

Rcpp::IntegerVector::getVerbose()
{
    Rcpp::IntegerVector output(1);
    output[0] = *verbose;
    return(output);
}

void samplingControl::setDebug(SEXP dbg)
{
    Rcpp::IntegerVector input(dbg);
    *debug = input[0];
}

Rcpp::IntegerVector samplingControl::getDebug()
{
    Rcpp::IntegerVector output(1);
    output[0] = *debug;
    return(output);
}

void samplingControl::setSliceWidths(SEXP widths)
{
    Rcpp::NumericVector inWidths(widths);
    if (inWidths.length() != 11)
    {
        Rcpp::Rcout << "Slice widths must have length 11.\n";
        throw(-1);
    }
    memcpy(sliceWidths, inWidths.begin(), 11*sizeof(double));
}

Rcpp::NumericVector samplingControl::getSliceWidths()
{
    Rcpp::NumericVector output(11);
    memcpy(output.begin(), sliceWidths, 11*sizeof(double));
    return(output);
}

void samplingControl::setSteadyStateConstraintPrecision(SEXP widths)
{
    Rcpp::NumericVector input(widths);
    *steadyStateConstraintPrecision = widths[0];
}

Rcpp::NumericVector samplingControl::getSteadyStateConstraintPrecision()
{
    Rcpp::NumericVector output(1);
    output[0] = *steadyStateConstraintPrecision;
    return(output);
}


void samplingControl::summary()
{
    Rcpp::Rcout << "I should really write a summary function to go here, huh?\n";
}
samplingControl::~samplingControl()
{
    delete[] xDim;
    delete[] zDim;
    delete[] X;
    delete[] Z;
    delete betaPriorPrecision;
}

RCPP_MODULE(mod_samplingControl)
{
    using namespace Rcpp;
    class_<samplingControl>( "samplingControl" )
    .constructor()
    .method("summary", &samplingControl::summary)
    .property("verbose", &samplingControl::getVerbose, &samplingControl::setVerbose, "Verbose level output?")
    .property("debug", &samplingControl::getDebug, &samplingControl::setDebug, "Debug level output?")
    .property("steadyStateConstraintPrecision", &samplingControl::getSteadyStateConstraintPrecision, &samplingControl::setSteadyStateConstraintPrecision, "Precision of normal constraint on compartment flows. Set to negative value to suppress.")
    .property("sliceWidths", &samplingControl::getSliceWidths, &samplingControl::setSliceWidths, "MCMC sampling parameters: S*, E*, R*, S0, I0, beta, betaPrs, rho, gamma_ei, gamma_ir, phi")
    .property("iterationStride", &samplingControl::getIterationStride, &samplingControl::setIterationStride, "Iterations between write operations to output file."); 
}



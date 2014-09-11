#ifndef SPATIALSEIR_SAMPLING_CONTROL
#define SPATIALSEIR_SAMPLING_CONTROL
#include <Rcpp.h>
#include<ModelContext.hpp>


using namespace Rcpp;
using namespace SpatialSEIR;

RCPP_EXPOSED_CLASS(samplingControl)
class samplingControl
{
    public:
        samplingControl();
        virtual void summary();
        ~samplingControl();

        virtual Rcpp::IntegerVector getIterationStride();
        virtual void setIterationStride(SEXP stride);

        virtual Rcpp::NumericVector getSteadyStateConstraintPrecision();
        virtual void setSteadyStateConstraintPrecision(SEXP prec);

        virtual Rcpp::IntegerVector getVerbose();
        virtual void setVerbose(SEXP vb);

        virtual Rcpp::IntegerVector getDebug();
        virtual void setDebug(SEXP dbg);

        virtual Rcpp::NumericVector getSliceWidths();
        virtual void setSliceWidths();

        int* iterationStride;
        int* verbose;
        int* debug;
        double* sliceWidths;
        double* steadyStateConstraintPrecision;       
};


#endif

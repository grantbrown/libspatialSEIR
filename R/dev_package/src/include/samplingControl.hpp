#ifndef SPATIALSEIR_SAMPLING_CONTROL
#define SPATIALSEIR_SAMPLING_CONTROL
#include <Rcpp.h>
#include<ModelContext.hpp>
#include<modelComponent.hpp>

using namespace Rcpp;
using namespace SpatialSEIR;

RCPP_EXPOSED_CLASS(samplingControl)
class samplingControl : public modelComponent
{
    public:
        samplingControl();
        virtual void summary();
        int getModelComponentType();
        ~samplingControl();

        virtual Rcpp::IntegerVector getIterationStride();
        virtual void setIterationStride(Rcpp::IntegerVector stride);

        virtual Rcpp::NumericVector getSteadyStateConstraintPrecision();
        virtual void setSteadyStateConstraintPrecision(Rcpp::NumericVector prec);

        virtual Rcpp::IntegerVector getVerbose();
        virtual void setVerbose(Rcpp::IntegerVector vb);

        virtual Rcpp::IntegerVector getDebug();
        virtual void setDebug(Rcpp::IntegerVector dbg);

        virtual Rcpp::NumericVector getSliceWidths();
        virtual void setSliceWidths(Rcpp::NumericVector widths);

        int* iterationStride;
        int* verbose;
        int* debug;
        double* sliceWidths;
        double* steadyStateConstraintPrecision;       
};


#endif

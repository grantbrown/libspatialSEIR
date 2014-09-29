#ifndef SPATIALSEIR_MODEL_INC
#define SPATIALSEIR_MODEL_INC
#include <Rcpp.h>

class dataModel;
class exposureModel;
class distanceModel;
class initialValueContainer;
class reinfectionModel;
class samplingControl;
class transitionPriors;

namespace SpatialSEIR
{
    class ModelContext;
}
class distanceModel;

using namespace Rcpp;
using namespace SpatialSEIR;

class spatialSEIRModel
{

    private:
        //Attributes: 
        ModelContext* context;
        std::string* chainOutputFile;
        int* verbose;
        int* debug;

    public: 

        spatialSEIRModel(SEXP outFileName);
        int buildSpatialSEIRModel(dataModel& dataModel_,
                                  exposureModel& exposureModel_,
                                  reinfectionModel& reinfectionModel_,
                                  distanceModel& distanceModel_,
                                  transitionPriors& transitionPriors_,
                                  initialValueContainer& initialValueContainer_,
                                  samplingControl& samplingControl_);
        // Simulation Functions
        virtual int setRandomSeed(int seedVal);
        virtual int simulate(int iters);
        virtual void setPredictionTraces();
        virtual int setTrace(int locationIndex);
        virtual int setTrace2(int locationIndex, int timeIndex);
        virtual void setDevice(int platformId, int deviceId);
        virtual void setCompartmentSamplingMode(int mode);
        virtual int getCompartmentSamplingMode();
        virtual void setParameterSamplingMode(int mode);
        virtual int getParameterSamplingMode();

        // Calculation Functions
        virtual int printDebugInfo();
        virtual int calculateS();
        virtual int calculateE();
        virtual int calculateI();
        virtual int calculateR();
        virtual int calculateP_SE(); 
        virtual int calculateP_SE2(int i, int j); 
        virtual int calculateP_SE_OCL();
        virtual double estimateEffectiveR0();
        virtual Rcpp::NumericVector estimateEffectiveR02(int t);
        virtual double estimateR0();
        virtual Rcpp::NumericVector estimateR02(int t);
        virtual int calculateP_RS();

        
        // Property Getter Functions
        // (All of this stuff is read only, should 
        // be changed only by calls to libspatialSEIR)
        virtual void updateSamplingParameters(double desiredRatio, double targetWidth, double proportionChange);
        virtual void printSamplingParameters();
        virtual void printAcceptanceRates();       
        virtual void printOCLSummary();
        virtual Rcpp::IntegerMatrix getS();
        virtual Rcpp::IntegerMatrix getE();
        virtual Rcpp::IntegerMatrix getI();
        virtual Rcpp::IntegerMatrix getR();

        virtual Rcpp::IntegerMatrix getS_star();
        virtual Rcpp::IntegerMatrix getE_star();
        virtual Rcpp::IntegerMatrix getI_star();
        virtual Rcpp::IntegerMatrix getR_star();

        virtual Rcpp::IntegerMatrix getY();

        virtual Rcpp::IntegerVector getS0();
        virtual Rcpp::IntegerVector getE0();
        virtual Rcpp::IntegerVector getI0();
        virtual Rcpp::IntegerVector getR0();

        virtual Rcpp::NumericMatrix getP_SE();
        virtual Rcpp::NumericVector getP_EI();
        virtual Rcpp::NumericVector getP_IR();
        virtual Rcpp::NumericVector getP_RS();
        virtual Rcpp::NumericVector getGenerationMatrix(int t);
        virtual Rcpp::NumericVector getIntegratedGenerationMatrix(int t);
        virtual Rcpp::NumericVector getBeta();
        virtual Rcpp::NumericVector getBetaP_RS();
        virtual Rcpp::NumericVector getRho();
        virtual Rcpp::NumericVector getPhi();


        virtual int getDebug();
        virtual void setDebug(int debug_);

        virtual int getVerbose();
        virtual void setVerbose(int verbose_);

        virtual int getUseDecorrelation();
        virtual void setUseDecorrelation(int val);

        virtual void standardizeDistanceMatrices();

        //Destructor
        ~spatialSEIRModel();
};

#endif

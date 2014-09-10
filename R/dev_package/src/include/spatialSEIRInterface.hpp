#ifndef SPATIALSEIR_INTERFACE_INC
#define SPATIALSEIR_INTERFACE_INC
#include <Rcpp.h>


namespace SpatialSEIR
{
    class ModelContext;
}
class distanceModel;

using namespace Rcpp;
using namespace SpatialSEIR;
struct LSSCout {};
extern LSSCout lssCout;

template <typename T>
    LSSCout& operator<< (LSSCout &s, const T &x){
        Rcpp::Rcout << x;
        return s;
    }

class spatialSEIRInterface
{

    private:
        //Attributes: 
        ModelContext* context;
        std::string* chainOutputFile;
        int* verbose;
        int* debug;

    public: 

        spatialSEIRInterface();
        int buildSpatialSEIRInterface(SEXP compMatDim,
                     SEXP xDim,
                     SEXP zDim,
                     SEXP xPrsDim,
                     SEXP S0_,
                     SEXP E0_,
                     SEXP I0_,
                     SEXP R0_,
                     SEXP Y,
                     SEXP Sstar, 
                     SEXP Estar, 
                     SEXP Istar, 
                     SEXP Rstar, 
                     SEXP offset_,
                     SEXP X_,
                     SEXP Z_,
                     SEXP X_pRS_,
                     const distanceModel& DM, 
                     SEXP rho_,
                     SEXP phi_,
                     SEXP priorAlpha_pEI_,
                     SEXP priorBeta_pEI_,
                     SEXP priorAlpha_pIR_,
                     SEXP priorBeta_pIR_,
                     SEXP priorAlpha_phi_,
                     SEXP priorBeta_phi_,
                     SEXP beta_,
                     SEXP betaPriorPrecision_,
                     SEXP betaPrs_,
                     SEXP betaPrsPriorPrecision_,
                     SEXP gamma_ei_,
                     SEXP gamma_ir_,
                     SEXP N_,
                     SEXP outFile,
                     SEXP iterationStride,
                     SEXP steadyStateConstraintPrecision_,
                     SEXP verboseFlag,
                     SEXP debugFlag,
                     SEXP sliceWidths,
                     SEXP reinfectionMode,
                     SEXP dataModel_);
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
        virtual double estimateR0();
        virtual double estimateR02(int t);
        virtual double estimateR03(int i, int t);
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

        const distanceModel* distModel;
 
        //Destructor
        ~spatialSEIRInterface();
};

#endif

\name{spatialSEIRModel}
\alias{spatialSEIRModel}
\title{R object pointing to a libSpatialSEIR model object.}
\description{
    The spatialSEIRModel object is created either implicitly, using qSEIR, qSpatialSEIR, 
    fit.qSEIR, or fit.qSpatialSEIR, or explicitly using the buildSEIRModel function. 
    Once created, it provides a set of methods used to draw additional MCMC samples, 
    investigate the current values of model parameters, and estimate several 
    reproductive numbers based on the current iteration. The object also
    allows the user to re-configure the MCMC samplers for the various parameter types. 
}




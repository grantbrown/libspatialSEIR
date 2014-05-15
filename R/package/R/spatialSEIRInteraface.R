spatialSEIRModel = function(compMatDim, xDim, zDim, xPrsDim, S0_, E0_, I0_, R0_, Sstar0, Estar0, Istar0, Rstar0, Sstar, Estar, Istar, Rstar, X_, Z_, X_pRS_, DistMat_, rho_, gamma_, priorAlpha_gamma_, priorBeta_gamma_, priorAlpha_pEI_, priorBeta_pEI_, priorAlpha_pIR_, priorBeta_pIR_, beta_, betaPriorPrecision_, betaPrs_, betaPrsPriorPrecision_, p_ei_, p_ir_, N_, outFile, logVarList, iterationStride,steadyStateConstraintPrecision, verboseFlag, debugFlag, sliceWidths)
{
    interface = new( spatialSEIRInterface )
    interface$buildSpatialSEIRInterface(compMatDim, xDim, zDim, xPrsDim, S0_, E0_, I0_, R0_, Sstar0, Estar0, Istar0, Rstar0, Sstar, Estar, Istar, Rstar, X_, Z_, X_pRS_, DistMat_, rho_, gamma_, priorAlpha_gamma_, priorBeta_gamma_, priorAlpha_pEI_, priorBeta_pEI_, priorAlpha_pIR_, priorBeta_pIR_, beta_, betaPriorPrecision_, betaPrs_, betaPrsPriorPrecision_, p_ei_, p_ir_, N_, outFile, logVarList,iterationStride,steadyStateConstraintPrecision, verboseFlag, debugFlag, sliceWidths)
    return(interface)
}

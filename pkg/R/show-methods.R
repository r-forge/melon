#############################################################
### show :: show a MELON 
#############################################################


setMethod("show","MELON",function(object) {
	cat("\n")
	cat("              MELON Normalized Data \n")
	cat("\n")
	cat("Reference Sample : ")
	cat(names(getNormFactors(object))[getRefId(object)])
	cat("\n")
	cat("\n")
	cat("Normalized Data (use getNormData() to extract all data)\n")
	print(head(getNormData(object)))
	cat("\n")
	cat("Normalization Factors (use getNormFactor() to extract Normalization Factors)\n")
	print(getNormFactors(object))
	cat("\n")
	cat("Indices stable loci (use getIdStableLoci() and getDataIdStableLoci())\n ")
	cat(getIdStableLoci(object))
	cat("\n")
})



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
	cat("Normalization Factors (use getNormFactors() to extract Normalization Factors)\n")
	print(getNormFactors(object))
	cat("\n")
	cat("Amount of optimal considered stable loci in the MELON-algorithm (use getIdStabLoci() for indices and getDataIdStabLoci() for data) : ")
	cat(length(getIdStabLoci(object)))
	cat("\n")
})



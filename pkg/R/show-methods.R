#############################################################
### show :: show a MELON 
#############################################################


setMethod("show","MELON",function(object) {
	cat("\n")
	cat("              object from class MELON \n")
	cat("\n")
	cat("Reference Sample : ")
	cat(names(getNormFactor(object))[getRefId(object)])
	cat("\n")
	cat("\n")
	cat("Normalized Data (use getNormData() to extract all data)\n")
	print(head(getNormData(object)))
	cat("\n")
	cat("Normalization Factors (use getNormFactor() to extract Normalization Factors)\n")
	print(getNormFactor(object))
	cat("\n")
	cat("Trim Intensity : ")
	cat(getTrimIntensity(object))
	cat("\n")
	cat("Trim Value : ")
	cat(getTrimValue(object))
	cat("\n")
})


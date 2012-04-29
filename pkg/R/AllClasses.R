######################
### CLASS : MELON ####
######################

setClass("MELON", 
	representation=representation(normData="matrix",normFactors="numeric",idStabLoci="numeric",dataStabLoci="matrix",stabFeatures="matrix",refId="numeric"),
	validity=function(object) {
		if (ncol(object@normData)!=length(object@normFactor)) {
			stop("Dimension of normalized data and factors do not match!")
		}
		return (TRUE)
	}
)


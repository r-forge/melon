######################
### CLASS : MELON ####
######################

setClass("MELON", 
	representation=representation(normData="matrix",normFactor="numeric",trimIntensity="numeric",trimValue="numeric",refId="numeric",medianResults="list"),
	validity=function(object) {
		if (ncol(object@normData)!=length(object@normFactor)) {
			stop("Dimension of normalized data and factors do not match!")
		}
		return (TRUE)
	}
)


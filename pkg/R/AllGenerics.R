setGeneric("getNormData",function(object){standardGeneric ("getNormData")})

setGeneric("getNormFactors",function(object){standardGeneric ("getNormFactors")})

setGeneric("getIdStabLoci",function(object){standardGeneric ("getIdStabLoci")})

setGeneric("getDataStabLoci",function(object){standardGeneric ("getDataStabLoci")})

setGeneric("getRefId",function(object){standardGeneric ("getRefId")})

setGeneric("getStabFeatures",function(object){standardGeneric ("getStabFeatures")})

setGeneric("normScatterPlot", function(object, nres = 5, datatype = c("Count","Continuous"), savedir=NULL, plotabline=T, tracefile=NULL) {
			standardGeneric("normScatterPlot")
})

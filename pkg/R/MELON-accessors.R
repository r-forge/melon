#########################################
### Accessors (MELON) 
#########################################

setMethod("getNormData","MELON", function(object){ 
	return(object@normData)
})

setMethod("getNormFactors","MELON", function(object){ 
	return(object@normFactors)
})

setMethod("getIdStabLoci","MELON", function(object){ 
	return(object@idStabLoci)
})

setMethod("getDataStabLoci","MELON", function(object){ 
	return(object@dataStabLoci)
})

setMethod("getRefId","MELON", function(object){ 
	return(object@refId)
})

setMethod("getStabFeatures","MELON", function(object){ 
	return(object@stabFeatures)
})



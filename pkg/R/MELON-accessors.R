#########################################
### Accessors (MELON) 
#########################################

setMethod("getNormData","MELON", function(object){ 
	return(object@normData)
})

setMethod("getNormFactor","MELON", function(object){ 
	return(object@normFactor)
})

setMethod("getTrimIntensity","MELON", function(object){ 
	return(object@trimIntensity)
})

setMethod("getTrimValue","MELON", function(object){ 
	return(object@trimValue)
})

setMethod("getRefId","MELON", function(object){ 
	return(object@refId)
})

setMethod("getMedianResults","MELON", function(object){ 
	return(object@medianResults)
})

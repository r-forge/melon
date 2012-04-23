setGeneric("getNormData",function(object){standardGeneric ("getNormData")})

setGeneric("getNormFactor",function(object){standardGeneric ("getNormFactor")})

setGeneric("getTrimIntensity",function(object){standardGeneric ("getTrimIntensity")})

setGeneric("getTrimValue",function(object){standardGeneric ("getTrimValue")})

setGeneric("getRefId",function(object){standardGeneric ("getRefId")})

setGeneric("getMedianResults",function(object){standardGeneric ("getMedianResults")})

setGeneric("setMedianResults<-",function(this,value){standardGeneric ("setMedianResults<-")})

setGeneric("makeMELON",
		function(object, ...)
		#function(object, mindatapointfraction=1/2, trimvalue=0, trimoption="sum", refsample=0, verbose=F)
			standardGeneric("makeMELON")
	  )

setGeneric("plotTrimValue", function(object, ...) {
			standardGeneric("plotTrimValue")
})

setGeneric("plotMedianRatio", function(object, ...) {
			standardGeneric("plotMedianRatio")
})

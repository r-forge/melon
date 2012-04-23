#########################################
### Constructor (MELON) 
#########################################

setMethod(f="initialize",signature="MELON",definition=function(.Object,normData,normFactor,trimIntensity, trimValue, refId){
	.Object@normFactor<-normFactor
	.Object@normData<-normData
	.Object@trimIntensity<-trimIntensity
	.Object@trimValue<-trimValue		
	.Object@refId<-refId		
	.Object@medianResults<-list()
	return (.Object)
})


#########################################
### Constructor (MELON) 
#########################################

setMethod(f="initialize",signature="MELON",definition=function(.Object,normData, normFactors, idStabLoci, dataStabLoci, stabFeatures, refId){
	.Object@normData<-normData
	.Object@normFactors<-normFactors
  .Object@idStabLoci<-idStabLoci
	.Object@dataStabLoci<-dataStabLoci		
	.Object@stabFeatures<-stabFeatures  
	.Object@refId<-refId		
	return (.Object)
})



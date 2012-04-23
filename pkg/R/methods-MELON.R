#########################################
### Plot Trim Value Selection (MELON) 
#########################################

setMethod("plotTrimValue","MELON",function(object,mindatapointfraction = 1/2) {
	    allresults<-getMedianResults(object)
	    refid<-getRefId(object)
	    if (refid != 1) mdresults<-matrix(ncol=length(getNormFactor(object)),nrow=nrow(allresults[[1]]))
	    else {
	    	mdresults<-matrix(ncol=length(getNormFactor(object)),nrow=nrow(allresults[[2]]))
	    }
	    for (i in (1:length(getNormFactor(object)))[-refid]){
		mdresults[,i]<-allresults[[i]][,3]
	    }
	    mdresults<-mdresults[,-refid]
	    trimval<-getTrimValue(object)
	    if(is.null(ncol(mdresults))==F){
		       mdsummed<-rowSums(mdresults)
	    }
	    else {    
		       mdsummed<-mdresults
	    }
	    if (is.null(nrow(mdresults))==F){
		mindatapoints<-ceiling(mindatapointfraction*nrow(mdresults))
	    }
	    else {
		mindatapoints<-ceiling(mindatapointfraction*length(mdresults))    
	    }
            xlims<-c(0,length(mdsummed))
            ylims<-c(0,1.1*max(mdsummed))
            plot(mdsummed,xlim=xlims,ylim=ylims,ylab="Number of trailing points",type="l",col="blue")
            points((length(mdsummed)-mindatapoints+2):length(mdsummed),mdsummed[(length(mdsummed)-mindatapoints+2):length(mdsummed)],pch=4,col="blue",cex=0.5)
            lines(c(trimval,trimval),c(min(mdsummed),max(mdsummed)),lty="dotted",col="red")
}
)


setMethod("plotMedianRatio","MELON",function(object,mindatapointfraction = 1/2, sampleId=1, legend=T, minY=NA, maxY=NA) {
	    allresults<-getMedianResults(object)
	    refid<-getRefId(object)
	    if (refid != sampleId) mdresults<-matrix(ncol=length(getNormFactor(object)),nrow=nrow(allresults[[sampleId]]))
	    else return ("Sample is reference sample");	
	    for (i in (1:length(getNormFactor(object)))[-refid]){
		mdresults[,i]<-allresults[[i]][,3]
	    }
	    mdresults<-mdresults[,-refid]
	    trimval<-getTrimValue(object)
	    if(is.null(ncol(mdresults))==F){
		       mdsummed<-rowSums(mdresults)
	    }
	    else {    
		       mdsummed<-mdresults
	    }
	    if (is.null(nrow(mdresults))==F){
		mindatapoints<-ceiling(mindatapointfraction*nrow(mdresults))
	    }
	    else {
		mindatapoints<-ceiling(mindatapointfraction*length(mdresults))    
	    }


            #plot local ratios
            results<-allresults[[sampleId]]
	    if (is.na(minY)) minY<- min(sqrt(results[,6]),na.rm=T)-0.1
	    if (is.na(maxY)) maxY<- max(sqrt(results[,6]),na.rm=T)+0.1

            plot(sqrt(results[,6]),ylab="sqrt(Ratio sample/reference)",main=paste("Sample ",colnames(getNormData(object))[sampleId]," (ref Sample : ", colnames(getNormData(object))[refid]," )",sep=""),ylim=c(minY,maxY))
 
            #plot final ratio (over different reference sample intensities)
            #title(main=colnames(countdata)[sampleId])
            lines(c(trimval,nrow(results)),sqrt(c(results[trimval,2],results[trimval,2])),lty="dashed")
            #plot sliding ratio
            lines(sqrt(results[,2]))
            
            #plot consecutive deviations from median
            maxrat<-results[2:nrow(results),6]
            maxrat<-max(maxrat[maxrat<Inf])
            
                lines(results[,3]*(sqrt(maxrat)/max(results[,3])),col="blue")
           
            #plot selected threshold for maximally enriched loci
            lines(c(trimval,trimval),c(0,maxrat),col="red",lty="dashed")
	    if (legend) legend("topright",c("Local Ratios","Final Ratio","Smoothed Ratio","Median Deviation","Threshold"),col=c("black","black","black","blue","red"),pch=c(1,-1,-1,-1,-1), lty=c(-1,2,1,1,2), merge=T)	
        
})


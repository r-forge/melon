#########################################
### Scatter Plots (MELON) 
#########################################

setMethod("normScatterPlot","MELON",function(object, nres = 5, datatype = c("Count","Continuous"), savedir=NULL, plotabline=TRUE, tracefile=NULL) {
  samples <- getNormData(object)
  refSample <- getRefId(object)
  stabLoci <- getIdStabLoci(object)
  if(is.null(tracefile)==F){sink(file = tracefile, append = TRUE, type = c("output", "message"),split = FALSE);writeLines(date())}
  if(is.null(savedir)==F){if(substr(savedir,nchar(savedir),nchar(savedir))!="/"){savedir<-paste(savedir,"/",sep="")}}
  if(nres<0|nres>100){warnings("nres should be between 0 and 100 (percentages), it has been put at 5%");nres<-5}
  firstplot<-1
  if((datatype%in%c("Count","Continuous"))==F){print("Datatype should be either 'Count' or 'Continuous', datatype has been set to 'Count'");datatype<-"Count"}
  if(datatype=="Continuous"){samples<-round(1000*samples)/1000}
  
  for (i in (1:ncol(samples))[-refSample]){
    
    writeLines(paste("Creating plot for sample ",colnames(samples)[i],sep=""))
    sampleredtable<-table(apply(samples[,c(refSample,i)],1,paste,sep="@",collapse="@"))
    allreddata<-cbind(t(apply(apply(as.matrix(strsplit(names(sampleredtable),split="@")),1,unlist),2,as.numeric)),as.vector(sampleredtable))
    colnames(allreddata)<-c(colnames(samples[,c(refSample,i)]),"counts")
    
    if(datatype=="Count"){allreddata[,1:2]<-log(allreddata[,1:2]+1)}
    
    if(is.null(savedir)==F){pdf(file=paste(savedir,"scatterplot",colnames(samples)[refSample],"_",colnames(samples)[i],".pdf",sep="",collapse=""))}
    
    normsamplered<-allreddata[,1:3]
    avs<-(normsamplered[,1]+normsamplered[,2])/2
    dists<-(avs-normsamplered[,1])*sqrt(2)
    normsamplered<-cbind(avs,normsamplered,dists)
    colnames(normsamplered)<-c("Average intensity",colnames(normsamplered)[2:4],"Dists")
    normsamplered<-normsamplered[order(normsamplered[,1]),]
    
    if(firstplot==0&is.null(savedir)==T){par(ask=T)}
    
    for (j in 1:nrow(normsamplered)){
      plot(normsamplered[j,2],normsamplered[j,3],xlim=c(0,max(normsamplered[,2:3])*1.1),ylim=c(0,max(normsamplered[,2:3])*1.1),cex=0.8+4*log(normsamplered[j,4])/log(max(normsamplered[,4])),pch=20,xlab="",ylab="")
      par(new=T)
    }
    if(is.null(stabLoci)==F){
      if(datatype=="Count"){
        for (j in 1:nrow(samples[stabLoci,])){
          plot(log(samples[stabLoci[j],refSample]+1),log(samples[stabLoci[j],i]+1),pch=20,col="blue",xlab="",ylab="",xlim=c(0,max(normsamplered[,2:3])*1.1),ylim=c(0,max(normsamplered[,2:3])*1.1),cex=0.8)
          par(new=T)    
        }
      } else {
        for (j in 1:nrow(samples[stabLoci,])){
          plot(samples[stabLoci[j],refSample],samples[stabLoci[j],i],pch=20,col="blue",xlab="",ylab="",xlim=c(0,max(normsamplered[,2:3])*1.1),ylim=c(0,max(normsamplered[,2:3])*1.1),cex=0.8)
          par(new=T)    
        }                
      }
    }
    
    if(datatype=="Count"){
      mtext(paste("ln ",colnames(normsamplered)[2]),side=1,font=2,line=3)
      mtext(paste("ln ",colnames(normsamplered)[3]),side=2,font=2,line=3)
    } else {
      mtext(colnames(normsamplered)[2],side=1,font=2,line=3)
      mtext(colnames(normsamplered)[3],side=2,font=2,line=3)            
    }
    
    if(plotabline==T){abline(0,1)}
    
    subsets<-round(seq(0,100,nres)*(nrow(normsamplered)-1)/100)+1
    subsetsums<-rep(0,length(subsets)-1)
    subsetsums<-cbind(subsetsums,subsetsums,subsetsums,subsetsums)
    colnames(subsetsums)<-c("x=y","dists","xwp","ywp")
    
    for (j in 1:(length(subsets)-1)){
      subsetsums[j,2]<-sum(normsamplered[(subsets[j]+1):subsets[j+1],5]*normsamplered[(subsets[j]+1):subsets[j+1],4])/sum(normsamplered[(subsets[j]+1):subsets[j+1],4])
      subsetsums[j,1]<-sum(normsamplered[(subsets[j]+1):subsets[j+1],1]*normsamplered[(subsets[j]+1):subsets[j+1],4])/sum(normsamplered[(subsets[j]+1):subsets[j+1],4])
    }
    
    subsetsums[,3]<-subsetsums[,1]-subsetsums[,2]/sqrt(2)
    subsetsums[,4]<-subsetsums[,1]+subsetsums[,2]/sqrt(2)
    #print(subsetsums[,2])
    points(subsetsums[,3],subsetsums[,4],pch=4,col="red",lwd=4)
    #y=-x+2a i.e. orthogonal to y=x, and through the point (a,a)
    #distance towards point (a,a) should be d, therefore: (x-a)?+(y-a)?=d? or (x-a)?+(-x+2a-a)?=d? or (x-a)?=d?/2 or x=a+abs(d)/sqrt(2) and y=a-abs(d)/sqrt(2) (sign )
    
    if(is.null(savedir)==F){dev.off()}
    firstplot<-0
  }
  
  if(is.null(savedir)==F){
    writeLines(paste("\nNormalization figures have been saved in ",savedir,"\n",sep=""))
  } else {par(ask=F)}
  if(sink.number()>0){sink()}
}
)



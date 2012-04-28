#setMethod("makeMELON",
 #   signature(object = "ANY"), function(object, method = c("Max","Clean"), regrwindow = 201, minC = 50, manC = NULL, maxL = 20000, refsample = 0, doplot = T, savedir = NULL, tracefile = NULL) {
makeMELON<- function(object, method = c("Max","Clean"), regrwindow = 201, minC = 50, manC = NULL, maxL = 20000, refsample = 0, doplot = T, savedir = NULL, tracefile = NULL) {
## object == countdata
	    stopifnot(is.data.frame(object) | is.matrix(object))
      
      if(is.null(tracefile)==F){sink(file = tracefile, append = TRUE, type = c("output", "message"),split = FALSE);writeLines(date())}  
      if(doplot==F & is.null(savedir)==F){writeLines("Warning: Save directory is ignored when doplot = FALSE\n")}
      if((method%in%c("Max","Clean"))==F){writeLines("Warning: 'method' should be either 'Max' or 'Clean', 'method' has been set to the latter\n");method<-'Clean'}
      if(is.null(savedir)==F){if(substr(savedir,nchar(savedir),nchar(savedir))!="/"){savedir<-paste(savedir,"/",sep="")}}
      if(is.null(savedir)==F){if(is.na(file.info(savedir)$size)){dir.create(savedir)}}
      if(method=="Max"&is.null(regrwindow)==F){writeLines("Local regression (regrwindow) is ignored for method = 'Max'")}
      if(method=="Clean"&is.null(regrwindow)==T){regrwindow<-201;writeLines("Local regression (regrwindow) should be a positive odd integer for method = 'Clean', and has been set to 201 (default)\n")}
      if(method=="Clean"&regrwindow%%2!=1){regrwindow<-201;writeLines("Local regression (regrwindow) should be a positive odd integer ans has been set to 201 (default)\n")}
      
      object<-as.matrix(object)
      
      #Identify sample with largest resolution (i.e. number of unique intensities)
      sampleresolution<-rep(NA,ncol(object))
      names(sampleresolution)<-colnames(object)
      for (i in 1:length(sampleresolution)){
        sampleresolution[i]<-length(round(unique(object[,i])*100)/100)
      }
      
      #Report which sample is used as reference sample, provides warning when other sample has higher resolution
      if (refsample==0){
        refid<-which(sampleresolution==max(sampleresolution))[1]
        writeLines(paste("Sample",colnames(object)[refid],"has the maximum of all sample resolutions and is therefore used as reference sample \n"))
      } else {
        refid<-refsample
        writeLines(paste("Sample",colnames(object)[refid],"(user provided) is considered as reference sample\n"))
      }
      
      #Determine maximally enriched loci
      writeLines("Identifying stable set of massively enriched loci\n")
      countorders<-apply(object,2,order,decreasing=T)
      
      allinds<-vector("list",ncol(object))
      tempind<-1:nrow(object)
      maxzero<-0
      
      for (i in 1:length(allinds)){
        maxzero<-max(c(maxzero,sum(object[,i]==0)))
        allinds[[i]]<-tempind[countorders[,i]]    
      }
      allinds<-do.call(cbind,allinds)
      
      evalthreshold<-nrow(allinds)-maxzero
      if(is.null(maxL)==F){if(maxL<=nrow(allinds)){evalthreshold<-maxL}else{evalthreshold<-nrow(allinds);writeLines("maxL should be lower than or equal to the number of rows in countdata. maxL has been altered automatically to this number of rows, but consider a further decrease for shorter calculation time.\n")}}
      
      if(evalthreshold!=nrow(allinds)){
        nonzeroinds<-setdiff((1:nrow(object)),unique(as.vector(allinds[(evalthreshold+1):nrow(allinds),])))
      } else {nonzeroinds<-(1:nrow(object))}
      nonzerocounts<-matrix(rep(0,2*length(nonzeroinds)),nrow=length(nonzeroinds),ncol=2)
      nonzerocounts[,1]<-nonzeroinds
      
      stabnumloci<-rep(0,evalthreshold)
      for (i in 1:evalthreshold){
        pickedup<-table(allinds[i,])
        pickeduploci<-which(nonzerocounts[,1]%in%names(pickedup))
        nonzerocounts[pickeduploci,2]<-nonzerocounts[pickeduploci,2]+pickedup[names(pickedup)%in%nonzerocounts[pickeduploci,1]]
        stabnumloci[i]<-sum(nonzerocounts[,2]==ncol(object))
      }
      mincutoff<-which(stabnumloci>=minC)[1]

      if(is.null(manC)==F){
        if(manC>max(stabnumloci)){manC<-max(stabnumloci);writeLines("manC was too large and has been reduced to the number of loci in the full dataset\n")}
        selcutoff<-which(stabnumloci>=manC)[i]
        writeLines("As manC was not NULL, this number of stable loci has been selected (minC is ignored)\n")    
      } else {
        stabvector<-stabnumloci/(1:length(stabnumloci))-((1:length(stabnumloci))/nrow(object))^(ncol(object)-1)
        if(method=="Max"){    
          selcutoff<-which.max(stabvector[mincutoff:evalthreshold])+mincutoff-1
        } else {
          firstdiff<-stabvector[2:length(stabvector)]-stabvector[1:(length(stabvector)-1)]
          firstdiffw<-rep(0,length(firstdiff))
          halfwindow<-(regrwindow-1)/2
          lt<-max(c(halfwindow+1,mincutoff))
          ut<-min(c(evalthreshold-halfwindow))
          for (i in lt:ut){
            firstdiffw[i]<-mean(firstdiff[(i-halfwindow):(i+halfwindow)])
          }
          selcutoff<-which.max(firstdiffw)
        }
      }      
      
      stabloci<-table(as.vector(allinds[1:selcutoff,]))
      stabloci<-as.numeric(names(stabloci[stabloci==ncol(object)]))
      writeLines(paste("The resulting amount of stable loci was ",length(stabloci),", this number can be altered by setting 'manC' manually\n",sep=""))
      
      #create plot
      if(doplot==1){
        if(is.null(savedir)==F){pdf(file=paste(savedir,"StabilityPlot.pdf",sep=""))}
        
        parset<-par(mar=c(5.1,4.1,2.1,4.1))
        plot(1:length(stabvector),stabvector,pch=20,axes=F,ylab="",xlab="",cex=0.8)
        axis(2,pretty(range(stabvector),10))
        mtext("Number of loci evaluated (L)",side=1,line=3,font=2)
        mtext("Consistent methylation enrichment (CME)",side=2,line=3,font=2)
        
        par(new=T)
        plot(1:length(stabvector),stabnumloci,type="l",axes=F,ylab="",xlab="",col="blue",lwd=4)
        axis(4,pretty(range(stabnumloci),20))
        mtext("Number of consistently methylated loci (C)",side=4,line=3,font=2)
        
        axis(1,pretty(range(1:length(stabvector)),20))
        
        lines(c(selcutoff,evalthreshold),c(length(stabloci),length(stabloci)),lty="dashed",col="green",lwd=2)
        lines(c(selcutoff,selcutoff),c(0,max(stabnumloci)),lty="solid",col="green",lwd=2)
        lines(c(mincutoff,evalthreshold),c(minC,minC),lty="dashed",col="red",lwd=2)
        
        box()
        par(parset)
        
        if(is.null(savedir)==F){dev.off()}
      }      
      
      #Determine normalization factors
      writeLines("Normalizing data\n")
      normfactors<-rep(1,ncol(object))
      names(normfactors)<-colnames(object)
      stabdata<-object[stabloci,]
      for(i in (1:ncol(stabdata))[-refid]){
        normfactors[i]<-median(stabdata[,i]/stabdata[,refid])
      }
      
      normdata<-object
      for (i in (1:length(normfactors))[-refid]){
        normdata[,i]<-normdata[,i]/normfactors[i]
      }
      rownames(normdata)<-rownames(object)
      stabmatrix<-cbind(1:length(stabvector),stabvector,stabnumloci)
      colnames(stabmatrix)<-c("Included amount of ranks","Stability indices","Number of stable loci")
      names(refid)<-colnames(object)[refid]
      
      writeLines("Normalization procedure completed. Normalized data and associated parameters are stored as a list.\n")
      if(is.null(savedir)==F){
        writeLines(paste("Stabilization figure has been saved in ",savedir,"\n",sep=""))
      }            
      if(sink.number()>0){sink()}
     return (new("MELON", normdata, normfactors,stabloci,object[stabloci,],stabmatrix,refid));	
    }




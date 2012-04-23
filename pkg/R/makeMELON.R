setMethod("makeMELON",
    signature(object = "ANY"),
    function (object, mindatapointfraction = 1/2, trimvalue = 0, trimoption = "sum", refsample = 0, verbose = F) {
	    stopifnot(is.data.frame(object) | is.matrix(object))
	    #Prepare variables
	    trimval<-trimvalue
	    countdata<-as.matrix(object)
    
	    #Identify sample with largest resolution (i.e. number of unique intensities)
	    sampleresolution<-rep(NA,ncol(countdata))
	    names(sampleresolution)<-colnames(countdata)
	    for (i in 1:length(sampleresolution)){
	        sampleresolution[i]<-length(unique(countdata[,i]))
	    }
	    if (verbose) { 
			writeLines("The different sample resolutions (i.e. max. number of unique intensity values) are:")
    			print(sampleresolution)
    			writeLines("\n")
    	    }
	    #Report which sample is used as reference sample, provides warning when other sample has higher resolution
	    if (refsample==0){
		refid<-which(sampleresolution==max(sampleresolution))[1]
		writeLines(paste("Sample",colnames(countdata)[refid],"has the maximum of all sample resolutions and is therefore used as reference sample \n"))
	    } else {
		refid<-refsample
		writeLines(paste("Sample",colnames(countdata)[refid],"(user provided) is considered as reference sample\n"))
		if ((refsample%in%which(sampleresolution==max(sampleresolution)))==F){
		    writeLines(paste("Sample",colnames(countdata)[which(sampleresolution==max(sampleresolution))],"has highest resolution and might be a better alternative"))
		    writeLines("\n")
		}
	    }

	    #Determine normalization factors for different thresholds (i.e. reference sample intensities) to define the maximally enriched loci (>= reference sample intensity)
	    #Data preparation
	    normfactors<-rep(1,ncol(countdata))
	    allresults<-as.list(rep("NULL",ncol(countdata)))
	    names(allresults)<-colnames(countdata)
	    allresults[[refid]]<-"Control sample"
	    mdresults<-matrix(ncol=ncol(countdata),nrow=sampleresolution[refid])
    
	    #Invoke endmedian function for each sample (+ reference sample)
	    for (i in (1:ncol(countdata))[-refid]){
		if (verbose) { writeLines(paste("Calculating normalization factor for sample",colnames(countdata)[i])) }
		results<-endmedian(countdata[,c(refid,i)])
		allresults[[i]]<-results
		mdresults[,i]<-results[,3]
	    }

	    mdresults<-mdresults[,-refid]
	    if (is.null(nrow(mdresults))==F){
		mindatapoints<-ceiling(mindatapointfraction*nrow(mdresults))
	    }
	    else {
		mindatapoints<-ceiling(mindatapointfraction*length(mdresults))    
	    }
    
    if (trimvalue==0){
        
        if (trimoption%in%c("sum","robust")==0){
            writeLines("Trimoption parameter should be either 'sum' or 'robust', analysis continues with option 'sum'\n")
            trimoption<-"sum"
        }
        
        if (trimoption=="sum"|is.null(ncol(mdresults))==T){
            allresults[[i]][,3]
            if(is.null(ncol(mdresults))==F){
               mdsummed<-rowSums(mdresults)
            }
            else {    
               mdsummed<-mdresults
            }
               
            mdmaxuser<-max(mdsummed[(length(mdsummed)-mindatapoints+1):length(mdsummed)])
        
            #From trimval on, consecutive deviations should be below or equal to maximum number of consecutive deviations in user defined region (mdmaxuser), and should fall outside of this region 
            devstarts<-which((mdsummed[1:(length(mdsummed)-mindatapoints+1)]-mdmaxuser)>0) #which points outside the user defined region have a higher amount of consecutive deviations
            if(length(devstarts)==0){devstarts<-0} #in case no such points are found
            trimvaltemp<-max(devstarts)+1 #first point for which number of consecutive deviations is below mdmaxuser
            trimvalinc<-which(mdsummed[1:(length(mdsummed)-1)]-mdsummed[2:length(mdsummed)]<0) #increasing mdsummed values
            trimval<-trimvalinc[trimvalinc>trimvaltemp][1]
            if(trimoption=="robust"){writeLines("Only two samples provided, trimoption 'robust' yields identical results as trimoption 'sum'\n")}
            } else {
            #Procedure is identical as above, but for each sample individually, finally the maximum value is selected             
            trimval<-0
            for (i in 1:ncol(mdresults)){
                mdsummed<-mdresults[,i]
                mdmaxuser<-max(mdsummed[(length(mdsummed)-mindatapoints+1):length(mdsummed)])
                devstarts<-which((mdsummed[1:(length(mdsummed)-mindatapoints+1)]-mdmaxuser)>0)
                if(length(devstarts)==0){devstarts<-0}
                trimvaltemp<-max(devstarts)+1 
                trimvalinc<-which(mdsummed[1:(length(mdsummed)-1)]-mdsummed[2:length(mdsummed)]<0) 
                trimvaltemp<-trimvalinc[trimvalinc>trimvaltemp][1]
                trimval<-max(c(trimval,trimvaltemp))
            }
        }
        
        #the summed variant infra is used for graphical reasons only: the robust variant is visible on sample specific plots
        if(is.null(ncol(mdresults))==F){
            mdsummed<-rowSums(mdresults)
            }
        else {    
            mdsummed<-mdresults
            }
        
        if (verbose) { writeLines(paste("The optimal trimvalue was",trimval,"corresponding with an intensity of",results[trimval,4],"(intensity value is usually close to trim value) \n")) }
        if (trimval>=(length(mdsummed)-mindatapoints+1)){
            writeLines("Warning! The optimal trim value falls within the defined region of maximally enriched loci.")
            writeLines(paste("Please decrease resolution fraction to at least",floor(100*(length(mdsummed)-trimval)/length(mdsummed))/100,"\n"))
        }
    }
    
    normfactors[refid]<-1
    for (i in (1:length(normfactors))[-refid]){
        normfactors[i]<-allresults[[i]][trimval,2]
    }
    names(normfactors)<-colnames(countdata)
    
    
    writeLines("\n")
    
    normfactors<-1/normfactors

    MELONResult<-new("MELON", t(normfactors*t(countdata)), normfactors,as.numeric(results[trimval,4]),trimval,refid)
    setMedianResults(MELONResult)<-allresults;
    return (MELONResult);	
    }
)


endmedian <-
function(rawsamples){
    #n x 2 matrix, with in first column ref sample
    #prepare variables

    #determine number of unique reference sample intensities
    refvals<-sort(unique(rawsamples[,1]))
    resl<-length(refvals) #resolution
    
    #prepare variables
    medians<-rep(NA,resl)
    meddevs<-rep(NA,resl)
    numvars<-rep(0,resl)
    diffs<-rep(0,resl)
    samplevals<-rep(NA,resl)
    
    #we start with the highest reference intensity value
    samplevals[resl]<-median(rawsamples[rawsamples[,1]>=refvals[resl],2])
    medians[resl]<-samplevals[resl]/refvals[resl]
    meddevs[resl]<-0
    numvars[resl]<-sum(rawsamples[,1]==refvals[resl])
    
    #for decreasing reference intensities, larger sets of sample subsets (putative maximally enriched loci) are created
    #by the addition of the ratio's for the specific ref. int. value (intsubset)
    #median ratios for the sample subset are calculated, as well as for the intsubset (local median ratio)
    #numvars stores the amount of the variables for that specific reference intensity(=nrow(intsubset))
    #diffs is only a temporary variable that tracks whether the local median ratio is higher or lower than the samplesubset median ratio
    #meddevs subsequently stores the number of consecutive (i.e. over the different reference intensity values) times the local median ratio is higher (or lower) than the samplesubset median ratio 

    for (i in ((resl-1):1)){
        samplesubset<-rawsamples[rawsamples[,1]>=refvals[i],]
        intsubset<-samplesubset[samplesubset[,1]==refvals[i],]
        medians[i]<-median(samplesubset[,2]/samplesubset[,1])
        numvars[i]<-length(intsubset)/2
        if (numvars[i]>1){
            diffs[i]<-sign(median(intsubset[,2]/intsubset[,1])-medians[i+1])
            samplevals[i]<-median(intsubset[,2])
        } else{
            diffs[i]<-sign(intsubset[2]/intsubset[1]-medians[i+1])
            samplevals[i]<-intsubset[2]
        }
        meddevs[i]<-which(abs(sign(diffs[(i+1):length(diffs)]-diffs[i:(length(diffs)-1)]))==1)[1]
    }
    
    #local median ratio calculation by local sample intensity/reference sample intensity (same denominator, so median(a/b) = median(a)/b)
    allres<-cbind(1:resl,medians,meddevs,refvals,samplevals,samplevals/refvals,numvars)
    colnames(allres)<-c("index","sample median ratio","consecutive deviations","reference median intensity","sample median intensity","local median ratio","number of variables (for reference intensity)")
    allres
}


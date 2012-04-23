#########################################
### Mutators (MELON) 
#########################################

setReplaceMethod("setMedianResults","MELON", function(this,value){ 
	this@medianResults<-value
	return (this)	
})

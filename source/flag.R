#Find out the length of the sequence (W/O deletion)
usualLength<-function(seq){
	max=length(seq[[1]])
	for(i in 2:length(seq)){
		if(max<length(seq[[i]])){
			max=length(seq[[i]])
		}
	}
return(max)
}

#Find out the allele which is shorter (with deletion)
flagAllele<-function(seq){
	ignoreList=vector('character')
	normal=usualLength(seq)
	for(i in names(seq)){
		if(length(seq[[i]])<normal){
			ignoreList=c(ignoreList, i)
		}
	}
return(ignoreList)
}


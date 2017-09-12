#' \code{usualLength} is used to find out 
#' the length of the sequence (W/O deletion).
#' @param seq a list of SeqFastadna. To keep it simple, 
#' use read.fasta from seqinr to import the fasta file.
#' @examples 
#' library(seqinr)
#' chlamydia <- read.fasta('Chlamydia.fas')
#' usualLength(chlamydia)
#' @return Will return the maximum length of all the allelic profiles.
usualLength<-function(seq){
	max=length(seq[[1]])
	for(i in 2:length(seq)){
		if(max<length(seq[[i]])){
			max=length(seq[[i]])
		}
	}
return(max)
}

#' \code{flagAllele} is used to find out a list of allelic profiles
#' that has been flagged and will not be included in computation
#' of minimum SNPs. 
#' @inheritParams usualLength
#' @examples
#' flagAllele(Chlamydia)
#' @return Will return a list of ignored allelic profiles.
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

#' \code{processAllele} is used to returned the processed allelic
#' profiles. 
#' @inheritParams usualLength
#' @examples 
#' processAllele(Chlamydia)
#' @return Will return the processed allelic profiles.
processAllele<-function(seq){
	ignore<-flagAllele(seq)
	processed=seq
	for (a in ignore){
		processed[[a]]<-NULL
	}
	return(processed)
}



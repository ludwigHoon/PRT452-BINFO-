#' \code{usualLength} is used to find out 
#' the length of the sequence (W/O deletion).
#' @param seq a list of SeqFastadna. To keep it simple, 
#' use read.fasta from seqinr to import the fasta file.
#' @return Will return the maximum length of all the allelic profiles.
#' @example tests/integrated.R
#' @export
usualLength<-function(seq){
	seqLen=list()
	for (i in 1:length(seq)){
		if( !(as.character(length(seq[[i]])) %in% names(seqLen)) ){
			seqLen[[ as.character(length(seq[[i]])) ]]=1
		}
		else{
			seqLen[[ as.character(length(seq[[i]])) ]]= seqLen[[ as.character(length(seq[[i]])) ]]+1
		}
	}
	maxInd=0
	max=0
	for (i in names(seqLen))	{
		if(seqLen[[i]]>max){
			maxInd=as.numeric(i)
			max=seqLen[[i]]
		}
	}
return(maxInd)
}

#' \code{flagAllele} is used to find out a list of allelic profiles
#' that has been flagged and will not be included in computation
#' of minimum SNPs. 
#' @inheritParams usualLength
#' @return Will return a list of ignored allelic profiles.
#' @export
flagAllele<-function(seq){
	ignoreList=vector('character')
	normal=usualLength(seq)
	for(i in names(seq)){
		if(length(seq[[i]])!=normal){
			ignoreList=c(ignoreList, i)
		}
	}
return(ignoreList)
}

#' \code{processAllele} is used to returned the processed allelic
#' profiles. 
#' @inheritParams usualLength
#' @return Will return the processed allelic profiles.
#' @export
processAllele<-function(seq){
	ignore<-flagAllele(seq)
	processed=seq
	for (a in ignore){
		print(a)
		processed[[a]]<-NULL
	}
	return(processed)
}

#' \code{flagPosition} is used to find out positions that 
#' will be ignored in calculation (either not A,C,G,T or '-')
#' @param proSeq Sequences after processed, i.e. all with the same length
#' @param dashIgnore whether to treat '-' as another type
#' @return Will return a list of positions that need to be ignored.
#' @export
flagPosition<-function(proSeq, dashIgnore=TRUE){
	acceptedChar=c("A","a","C","c","T","t","G","g")
	ignoredPos=vector('numeric')
	if (dashIgnore==FALSE) acceptedChar=c(acceptedChar, "-")
	
	for (i in 1:length(proSeq[[1]])){
		for (i2 in 1:length(proSeq)){
			if( !is.element(proSeq[[i2]][i],acceptedChar) ){
				ignoredPos=c(ignoredPos,i)
				break
			}

		}
	}
	print(paste('Ignored Positions:', toString(ignoredPos)))
	return(ignoredPos)
}

#' \code{readFasta} is used to read fasta file
#' @import seqinr
#' @param file file path
#' @return Will return fastaDNA.
#' @export
readFasta<-function(file){
    return(read.fasta(file))
}



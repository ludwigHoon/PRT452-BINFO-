#' \code{similar.percent} is used to find the calculate the percentage of 
#' similarity at alleles.
#' @param seq a list of SeqFastadna. To keep it simple, 
#' use read.fasta from seqinr to import the fasta file.
#' @param ref the specific allele to be identified.
#' @return Will return the a list of SNPs that can be used 

similar.percent <-function(seq, ref){
	percentList=list()
	target=match(ref, getName(seq))
	
	if(is.na(target)){
	print("Can't find reference in allelic profiles")
	return(NULL)
	}
	
	numSeq=length(seq)
	numGene=length(seq[[target]])
	
	for(position in 1:numGene){
		same=0
		for (profile in 1:numSeq){
			if(profile!=target){
				if(seq[[target]][position]==seq[[profile]][position]){
				same=same+1
				}
			}
		}
		percent= (1-(same/(length(seq)-1)))*100
		percentList[[position]]=list(position=position, percent=percent)
	}
	return (percentList)
}

#' \code{present.percent} is used to find present and filter the 
#' similarity calculated using \code{similar.percent}.
#' @import rlist
#' @param result the result from \code{similar.percent}.
#' @param percent minimum percentage to be included
#' @param number number of results to be displayed
#' @return Will return the a list of SNPs (as specified) that can be used 
#' and the associated percentage at the particular location.

present.percent <-function(result, percent=100, number=100){
	if(percent>100){
	print('Percent error')
	return(NULL)
	}
	result=list.sort(result, (percent))
	if(number<length(result)){
		result=result[1:number]
	}
	result2=list()
	found=list.which(result, percent>=percent)
	for(a in found){
		result2[[a]]=result[[a]]
	}
	return(result2)
}
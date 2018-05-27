#' \code{percent.calculate} is used to find the calculate the percentage of 
#' similarity at alleles.
#' @param seq the pattern generated from \code{percent.pattern}
#' @param step the pattern generated either from previous function invocation or empty list
#' @param tName the position to add to the pattern
#' @return Will return the pattern for analysis later
#' @export
percent.pattern<-function(seq, step, positions){
    for (a in seq){
        b=getName(a)
        newSeq=""
        for (p in positions){
            newSeq=paste(newSeq, a[p], sep='')
        }
        step[[b]]=paste(step[[b]], newSeq, sep='')
    }
    return(step)
}

#' \code{percent.calculate} is used to calculate the percentage of dissimilarity.
#' @param pattern is the structure generated from \code{percent.pattern}.
#' @param tNames is the name of the reference allele.
#' @return Will returns percentage of dissimilarity.
#' @export
percent.calculate<-function(pattern, tNames){
    for (n in tNames){
		if (is.null(pattern[[n]])){
			print("Target not found")
		} 
	}
        targetSeqs=vector('character')
		for (n in tNames){
			targetSeqs=c(targetSeqs, pattern[[n]])
		}
		sum=0
        for (a in pattern){
            if (a %in% targetSeqs){
                sum=sum+1
            }
        }
        sum = sum - length(targetSeqs)
        freq = (1 - (sum/(length(pattern)-length(tNames))))*100
        return(round(freq, 5))
    
}

#' \code{similar.percent} is used to calculate the percentage of dissimilarity and list the position with highest percentage.
#' @import seqinr
#' @import rlist
#' @param seq is fastaDNA object to analysed.
#' @param targets is the names of the reference alleles
#' @param level is the number of positions.
#' @param included is the included position for analysed, these will force the computation
#' to compute the percentage of dissimilarity at the position no matter what
#' @param excluded is the positions that are excluded from computation
#' @return Will returns percentage of dissimilarity and the position.
#' @export
similar.percent<-function(seq, targets, level=1, included=NULL, excluded=NULL){
    for (n in targets){
		if (is.null(seq[[n]])){
			return(NULL)
		} 
	}
    N=length(seq)
	result=list()
	lvl=1
	type=list()
    remnant=list()
	explored=vector('numeric')
	if(!is.null(included)){
		for (inc in included){
			type=percent.pattern(seq, type, inc)
			percent=percent.calculate(type, targets)
			explored=c(explored, inc)
			result[[lvl]]=list(position=setdiff(explored, excluded), percent=percent)
			lvl=lvl+1
		}
	}
	for(a in 1:level){
	curRes=list()
	explored=c(explored, excluded)
	for(position in 1:length(seq[[1]])){
		if (position %in% explored){next}
		ctype=percent.pattern(seq, type, position)
		percent=percent.calculate(ctype, targets)
		curRes[[position]]=list(position=position, value=percent)	
	}	
	for (a in explored){ curRes[[a]]=list(position=a, value=0)}
		pos=list.order(curRes, (value))[1]
		index=as.numeric(curRes[[pos]]['value'])
		actPos=curRes[[pos]]$position
        explored=c(explored, actPos)
		result[[lvl]]=list(position=setdiff(explored, excluded),percent=index)
		type=percent.pattern(seq, type, actPos)
		lvl=lvl+1
	}
	return(result)
}

#' \code{percent.residual} is used to find the residual.
#' @param seq is fastaDNA object to analysed.
#' @param positions is the positions from the result.
#' @param targets is the name of the reference allele.
#' @param excluded is the excluded positions
#' @return Will returns the residual positions.
#' @export
percent.residual<-function(seq, targets, positions=NULL, excluded=NULL){
 	for (n in targets){
		if (is.null(seq[[n]])){
			return(NULL)
		} 
	}
    
	result=list()
	for (a in 1:length(seq[[1]])){
		type=list()
		if(a %in% c(positions, excluded)){next}
		type=percent.pattern(seq, type, c(positions, a))
		percent=percent.calculate(type, targets)
		if(percent==100){
			result=c(result, a)
		}
	}
	return (list(prepended=toString(positions), residual=toString(result)))
}

#' \code{branch.percent} is used to calculate 1 or more percentage of dissimilarity and list the position with highest percentage.
#' @import seqinr
#' @import rlist
#' @param seq is fastaDNA object to analysed.
#' @param targets is the name of the reference allele.
#' @param level is the number of positions.
#' @param included is the included position for analysed, these will force the computation
#' to compute the percentage at the position no matter what
#' @param excluded is the positions that are excluded from computation
#' @param numRes is the number of result to be returned.
#' @return Will returns percentage of dissimilarity and the position
#' as well as the residual positions (i.e. those positions with 100% on their own).
#' @export
branch.percent<-function(seq, targets, level=1, included=NULL, excluded=NULL, numRes=1){
    if(numRes<=0){numRes=1}
    result=list()
	num=1
	res=similar.percent(seq, targets, level, included, excluded)
	result[[paste('result', num, sep= ' ')]]=res
	excluded=c(excluded, res[[1]]$'position')
	num=num+1

	while(num<=numRes){
		res=similar.percent(seq, targets, level, included, excluded)
		result[[paste('result', num, sep= ' ')]]=res
		excluded=c(excluded, res[[1]]$'position')
		num=num+1	
	}
	return(result)
}

#' \code{present.percent} is used to present the result of percent calculation.
#' @param result the result from \code{branch.percent}
#' @return Will print out the result
#' @export
present.percent<-function(result){
    print("Result:")
    print(result$'result')
    print("Residual")
    print(result$'remnant')
}
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
#' @param seq is fastaDNA object to analysed.
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
#' @param last force the function to return whole current result when it's the last round = to find the residual
#' @return Will returns percentage of dissimilarity and the position.
#' @export
similar.percent<-function(seq, targets, level=1, included=NULL, excluded=NULL, last=FALSE){
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
	appended=list()
	explored=vector('numeric')
	if(!is.null(included)){
		for (inc in included){
			type=percent.pattern(seq, appended, inc)
			percent=percent.calculate(type, targets)
			for(al in 1:length(seq)){
			appended[[getName(seq[[al]])]]=paste(appended[[getName(seq[[al]])]], seq[[al]][inc], sep="")
			}
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
		type=percent.pattern(seq, appended, position)
		percent=percent.calculate(type, targets)
		curRes[[position]]=list(position=position, value=percent)	
		}	
		explored=sort(explored, decreasing=TRUE)
		for (a in explored){ curRes[[a]]<-NULL}
		position=list.order(curRes, (value))[1]
		index=as.numeric(curRes[[position]]['value'])
		position=curRes[[position]]$position
        explored=c(explored, position)
		result[[lvl]]=list(position=setdiff(explored, excluded),percent=index)
		lvl=lvl+1
		for(al in 1:length(seq)){
			appended[[getName(seq[[al]])]]=paste(appended[[getName(seq[[al]])]], seq[[al]][position], sep="")
			}
	}
	if(last){
		sorted=list.order(curRes,(value))
		residual=vector('numeric')
		for (p in sorted){
			if(as.numeric(curRes[[p]]['value'])==100){
				residual=c(residual, curRes[[p]]$position)
			}else{
				break
			}
		}
		return(list(result=result, prepended=setdiff(explored, c(excluded, position)), residual=setdiff(residual, position)))
	}
	return(result)
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
	res=similar.percent(seq, targets, level, included, excluded, (num==numRes))
	if(num==numRes){
			result[[paste('result', num, sep= ' ')]]=res$'result'
			excluded=c(excluded, res$'result'[[1]]$'position')
			remnant=res$'residual'
		}else{
			result[[paste('result', num, sep= ' ')]]=res
			excluded=c(excluded, res[[1]]$'position')
		}
	num=num+1
	while(num<=numRes){
		res=similar.percent(seq, targets, level, included, excluded, (num==numRes))
		if(num==numRes){
			result[[paste('result', num, sep= ' ')]]=res$'result'
			excluded=c(excluded, res$'result'[[1]]$'position')
			remnant=list(prepended=res$'prepended', position=res$'residual')
			num=num+1
		}else{
			result[[paste('result', num, sep= ' ')]]=res
			excluded=c(excluded, res[[1]]$'position')
			num=num+1
		}
	}
    RES=list(result=result, remnant=remnant)
	return(RES)
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
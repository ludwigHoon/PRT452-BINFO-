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
#' @param tName is the name of the reference allele.
#' @return Will returns percentage of dissimilarity.
#' @export
percent.calculate<-function(pattern, tName){
    if (is.null(pattern[[tName]])){
        print("Target not found")
    } else{
        targetSeq=pattern[[tName]]
        sum=0
        for (a in pattern){
            if (a==targetSeq){
                sum=sum+1
            }
        }
        sum = sum - 1
        freq = (1 - (sum/(length(pattern)-1)))*100
        return(round(freq, 5))
    }
}

#' \code{similar.percent} is used to calculate the percentage of dissimilarity and list the position with highest percentage.
#' @import seqinr
#' @import rlist
#' @param seq is fastaDNA object to analysed.
#' @param target is the name of the reference allele
#' @param level is the number of positions.
#' @param included is the included position for analysed, these will force the computation
#' to compute the percentage of dissimilarity at the position no matter what
#' @param excluded is the positions that are excluded from computation
#' @return Will returns simpson's index and the position.
#' @export
similar.percent<-function(seq, target, level=1, included=NULL, excluded=NULL){
    if(is.null(seq[[target]])){
        return(NULL)
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
			percent=percent.calculate(type, target)
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
		percent=percent.calculate(type, target)
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
	return(result)
}

#' \code{branch.percent} is used to calculate 1 or more percentage of dissimilarity and list the position with highest percentage.
#' @import seqinr
#' @import rlist
#' @param seq is fastaDNA object to analysed.
#' @param target is the name of the reference allele.
#' @param level is the number of positions.
#' @param included is the included position for analysed, these will force the computation
#' to compute the percentage at the position no matter what
#' @param excluded is the positions that are excluded from computation
#' @param numRes is the number of result to be returned.
#' @return Will returns percentage of dissimilarity and the position
#' as well as the residual positions (i.e. those positions with 100% on their own).
#' @export
branch.percent<-function(seq, target, level=1, included=NULL, excluded=NULL, numRes=1){
    if(numRes<=0){numRes=1}
    result=list()
	num=1
	res=similar.percent(seq, target, level, included, excluded)
	result[[paste('result', num, sep= ' ')]]=res
	num=num+1
	excluded=c(excluded, res[[1]]$'position')
	
	while(num<=numRes){
		res=similar.percent(seq, target, level, included, excluded)
		result[[paste('result', num, sep= ' ')]]=res
		num=num+1
		excluded=c(excluded, res[[1]]$'position')
	}
    remnant=vector('numeric')
    for (p in 1:length(seq[[1]])){
        if(p %in% excluded){next}
        type=percent.pattern(seq, list(), p)
        percent=percent.calculate(type, target)
        if(percent==100){
            remnant=c(remnant, p)
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
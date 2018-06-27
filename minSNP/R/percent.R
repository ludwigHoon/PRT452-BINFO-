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
			return(-1)
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

distinct<-function(seq, tNames, exc=NULL){
	for (name in tNames){
		if(! name %in% names(seq)){
			print('Sequence not found')
			return(NULL)
		}
	}
	newExc = list()
	i<-1
	for (tName in tNames){
		for (p in 1:length(seq[[tName]])){
			if(p %in% exc){next}
			if (percent.calculate(percent.pattern(seq, list(), p) , tName) ==100){
				newExc[[i]]<-p
				i<-i+1
			}
		}
	}
	newExc = c(newExc, exc)
	return(newExc)
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
	if(level>0){
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
	}
	return(result)
}

#' \code{percent.residual} is used to find the residual.
#' @param seq is fastaDNA object to analysed.
#' @param positions is the positions from the result.
#' @param targets is the name of the reference allele.
#' @return Will returns the residual alleles.
#' @export
percent.residual<-function(seq, targets, positions=NULL){
	if (is.null(seq[[targets]])){
		return(NULL)
	} 

	if (is.null(positions)){
		return(NULL)
	}
	
	result=list()
	type=percent.pattern(seq, list(), positions)

	for (t in names(type)){
		if (t==targets){
			next
		}
		if (type[t] == type[[targets]]){
			result=c(result, t)
		}
	}
	if(length(result)>0){
		return (toString(result))
	}else{
		return(NULL)
	}
	
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
#' @param target is the name of the allele targeted
#' @param seq is fastaDNA object to analysed.
#' @return Will print out the result
#' @export
present.percent<-function(seq, target, result){
	numRes=length(result)
	if(numRes>=1){
		level=length(result[[1]])
	}
	for (a in 1: numRes){
		print(paste("Result ", a, ":", sep=''))
   		print(result[[a]])
		print("Residual:")
		res=percent.residual(seq, target, result[[a]][[level]]$position)
		print(res)
		cat('----------\n\n')
	}
}
#' \code{result.csv} is used to output the result to csv file. 
#' @param result is from \code{branch.percent} 
#' @param target is the name of the allele targeted
#' @param filename is the name of the output csv file, default to datetime of when the file is generated
#' @return will generate a csv file with the filename in the current directory
#' @export 
result.csv<-function(result, target, filename=NULL){
	#file name is generated automatically from the current system date 
	if (filename ==NULL){filename=format(Sys.time(), "%Y-%m-%d_%H:%M");}
	filename=paste(filename,".csv", sep = "")

	#Convert the result into suitable format for export to csv
	for (a in names(result)){

	}

	#Write to CSV file

	#Adding the details & metadata of the analysis
	cat("\n", "Details", target, sep = "\n", append = TRUE, file = filename)

}
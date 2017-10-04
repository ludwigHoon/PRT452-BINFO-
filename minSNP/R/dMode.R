#' \code{simpson.calculate} is used to calculate the simpson's index given a pattern.
#' @param pattern is a pattern, can be a vector or a list.
#' @param N is the total number of entities that are in the pattern.
#' @return Will returns the simpson's index for the pattern.
simpson.calculate<-function(pattern, N){
	tsum=0
	for(t in pattern){
		tsum=tsum+(t*(t-1))
	}
	if((tsum==0)||(N==1)){
		simpIndex=1
	}
	else{
		simpIndex=1-(tsum/(N*(N-1)))
	}
	return(simpIndex)
}

#' \code{simpson.pattern} is used to generate pattern for calculation at a later stage.
#' @import seqinr
#' @param seq is fastaDNA object to analysed.
#' @param position is the position of the sequences used to generate pattern.
#' @param appended is the pattern that the current operation will appended onto
#' @return Will returns the generated pattern.
simpson.pattern<-function(seq, position, appended=NULL)
{
	type=list()
			for(al in 1:length(seq)){
				name= paste(appended[[getName(seq[[al]])]], seq[[al]][position], sep="")
				if (is.null(type[[name]])){
					type[[name]]=1}
				else{
					type[[name]]=type[[name]]+1
				}
		}
	return(type)
}

#' \code{similar.simpson} is used to calculate the simpson index and list the position with highest index.
#' @import seqinr
#' @import rlist
#' @param seq is fastaDNA object to analysed.
#' @param level is the number of positions.
#' @param included is the included position for analysed, these will force the computation
#' to compute the simpson's index at the position no matter what
#' @param excluded is the positions that are excluded from computation
#' @return Will returns simpson's index and the position.
similar.simpson <-function(seq, level=1, included=NULL, excluded=NULL){
	N=length(seq)
	result=list()
	lvl=1
	type=list()
	appended=list()
	explored=vector('numeric')
	if(!is.null(included)){
		for (inc in included){
			type=simpson.pattern(seq, inc, appended)
			simpIndex=simpson.calculate(type, N)
			for(al in 1:length(seq)){
			appended[[getName(seq[[al]])]]=paste(appended[[getName(seq[[al]])]], seq[[al]][inc], sep="")
			}
			explored=c(explored, inc)
			result[[lvl]]=list(position=inc,index=simpIndex)
			lvl=lvl+1
		}
	}
	for(a in 1:level){
	curRes=list()
	explored=c(explored, excluded)
	for(position in 1:length(seq[[1]])){
		if (position %in% explored){next}
		type=simpson.pattern(seq, position, appended)
		simpIndex=simpson.calculate(type, N)
		curRes[[position]]=list(position=position, value=simpIndex)	
		}	
		explored=sort(explored, decreasing=TRUE)
		for (a in explored){ curRes[[a]]<-NULL}
		position=list.order(curRes, (value))[1]
		index=as.numeric(curRes[[position]]['value'])
		position=curRes[[position]]$position
		result[[lvl]]=list(position=position,index=index)
		lvl=lvl+1
		for(al in 1:length(seq)){
			appended[[getName(seq[[al]])]]=paste(appended[[getName(seq[[al]])]], seq[[al]][position], sep="")
			}
		explored=c(explored, position)
	}
	return(result)
	
}

#' \code{branch.simpson} is used to calculate 1 or more simpson index and list the position with highest index.
#' @import seqinr
#' @import rlist
#' @param seq is fastaDNA object to analysed.
#' @param level is the number of positions.
#' @param included is the included position for analysed, these will force the computation
#' to compute the simpson's index at the position no matter what
#' @param excluded is the positions that are excluded from computation
#' @param numRes is the number of result to be returned.
#' @return Will returns simpson's index and the position.
branch.simpson <-function(seq, level=1, included=NULL, excluded=NULL, numRes=1){
	if(numRes<=0){numRes=1}
	result=list()
	num=1
	res=similar.simpson(seq, level, included, excluded)
	result[[paste('result', num, sep= ' ')]]=res
	num=num+1
	
	while(num<=numRes){
		excluded=c(excluded, res[[1]]$'position')
		res=similar.simpson(seq, level, included, excluded)
		result[[paste('result', num, sep= ' ')]]=res
		num=num+1
	}
	return(result)
}

#' \code{present.simpson} is used to present the result from the calculation of simpson's index.
#' @import seqinr
#' @param seq is fastaDNA object to analysed.
#' @param result is the result from \code{branch.simpson}.
#' @return Will returns the presentation of the result.
present.simpson <-function(seq, result){
	numRes=length(result)
	level=length(result[[1]])
	output=list()
	for (num in 1:numRes){
	appended=list()
	descrip='At position:'
	for (lvl in 1:level){
			po=result[[num]][[lvl]]$'position'
			descrip=paste(descrip, po, sep='-')
			for (al in 1:length(seq)){
				appended[[getName(seq[[al]])]]=paste(appended[[getName(seq[[al]])]], seq[[al]][po], sep="")
			}
	}
	type=list()
	for (a in 1:length(appended)){
		type[[appended[[a]]]]=c(type[[appended[[a]]]], names(appended)[a])
	}
	type[['Description']]=descrip
	type[['Index']]=result[[num]][[level]]$'index'
	output[[num]]=type
	}
	return(output)
}		
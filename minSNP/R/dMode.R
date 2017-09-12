#' \code{similar.simpson} is used to find the calculate the percentage of 
#' similarity at alleles.
#' @inheritParams usualLength
#' @params ref the specific allele to be identified.
#' @examples 
#' similar.simpson(chlamydia)
#' @return Will return the a list of SNPs that can be used 
#' and the associated percentage at the particular location.
library(rlist)
similar.simpson <-function(seq, level=1, resFromUp=NULL){
	N=length(seq)
	result=list()
	if ((level==1)&&is.null(resFromUp)){
		curRes=list()
		for(position in 1:length(seq[[1]])){
		type=list()
			for(al in 1:length(seq)){
				if (is.null(type[[seq[[al]][position]]])){
					type[[seq[[al]][position]]]=1}
				else{
					type[[seq[[al]][position]]]=type[[seq[[al]][position]]]+1
				}		
		}
		tsum=0
		for(t in type){
			tsum=tsum+(t*(t-1))
		}
		simpIndex=1-(tsum/(N*(N-1)))
		curRes[[position]]=list(value=simpIndex)
	}
	position=list.order(curRes, (value))[1]
	index=as.numeric(curRes[[list.order(curRes, (value))[1]]]['value'])
	
	result[[1]]=list(position=position, index=index)
	}
	
	else{
		if (length(resFromUp) < (level-1)){
			result=similar.simpson(seq, (level-1), resFromUp)
			resFromUp=result
		}
		
		explored=vector('numeric')
		appended=list()
		for (a in 1:length(resFromUp)){
			explored=c(explored, resFromUp[[a]]$'position')
			po=resFromUp[[a]]$'position'
			for (al in 1:length(seq)){
				appended[[getName(seq[[al]])]]=paste(appended[[getName(seq[[al]])]], seq[[al]][po], sep="")
			}
		}
		#print(appended)
		curRes=list()
		for(position in 1:length(seq[[1]])){
		type=list()
		if (position %in% explored){next}
			for(al in 1:length(seq)){
				name= paste(appended[[getName(seq[[al]])]], seq[[al]][position], sep="")
				if (is.null(type[[name]])){
					type[[name]]=1}
				else{
					type[[name]]=type[[name]]+1
				}
		}
		tsum=0
		for(t in type){
			tsum=tsum+(t*(t-1))
		}
		simpIndex=1-(tsum/(N*(N-1)))
		curRes[[position]]=list(position=position, value=simpIndex)	
		}
		explored=sort(explored, decreasing=TRUE)
		for (a in explored){ curRes[[a]]<-NULL}
		position=list.order(curRes, (value))[1]
		index=as.numeric(curRes[[position]]['value'])
		position=curRes[[position]]$position
		result[[level]]=list(position=position,index=index)
	}
	return(result)
}

present.simpson <-function(seq, result){
	level=length(result)
	appended=list()
	descrip='At position:'
	for (lvl in 1:level){
			po=result[[lvl]]$'position'
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
	return(type)
}
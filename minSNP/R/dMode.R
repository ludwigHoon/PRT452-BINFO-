#' \code{similar.simpson} is used to find the calculate the percentage of 
#' similarity at alleles.
#' @inheritParams usualLength
#' @params ref the specific allele to be identified.
#' @examples 
#' similar.simpson(chlamydia)
#' @return Will return the a list of SNPs that can be used 
#' and the associated percentage at the particular location.
library(rlist)

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
	# excl=''
	# incl=''
	# for(exc in excluded){excl=paste(excl, exc)}
	# for(inc in included){incl=paste(incl, inc)}
	# result['excluded']=excl
	# result['included']=incl
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
	#inc=paste('included:', result['included'])
	#exc=paste('excluded:', result['excluded'])
	#descrip=paste(descrip, inc, exc, sep='\n')
	type[['Description']]=descrip
	return(type)
}

simpson.calculate<-function(pattern, N){
	tsum=0
	for(t in pattern){
		tsum=tsum+(t*(t-1))
	}
	simpIndex=1-(tsum/(N*(N-1)))
	return(simpIndex)
}

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

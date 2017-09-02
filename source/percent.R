library(seqinr)

similar <-function(seq){
	a=seq[[1]][1598]
	diff=0
	same=0
	for( i in 2:length(seq)){
		if (seq[[i]][1598]!=a){
		diff<-diff+1}
	}
	for( i in 2:length(seq)){
		if (seq[[i]][1598]==a){
		same<-same+1}
	}
	print(c(diff, same, diff/(length(seq)-1)*100, same/(length(seq)-1)*100))
}


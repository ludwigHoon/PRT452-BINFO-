sample.case1 <-function(){
	#Read the file
	Chlamydia <- read.fasta(file='../data/Chlamydia_mapped.txt')
	
	#STEP 1. Process the file
	Chlamydia <- processAllele(Chlamydia)
	untempered <-read.fasta(file='../data/Chlamydia_mapped.txt')
	
	#Since Chlamydia is normal
	checkIdentical(Chlamydia, untempered)
	
	#PERCENT MODE
	result=similar.percent(Chlamydia, 'A_D213')
	present=present.percent(result, 98, 100)
	
	#All result should have percent higher or equal to 98
	for (a in 1:100){
	checkTrue(present[[a]]$percent>=98)
	}
	
	#SIMPSON MODE
	result=branch.simpson(Chlamydia, level=1, numRes=3)
	output=present.simpson(Chlamydia, result)
	
	#Should have 3 results
	checkTrue(length(output)==3)
	
	checkEquals(output[[1]]$'Index',  0.7344, tolerance=0.00016)
	Description=paste('At position:', '1988', sep='-')
	checkEquals(output[[1]]$'Description', Description)
	
	checkEquals(output[[2]]$'Index', 0.7318, tolerance=0.00016)
	Description=paste('At position:', '2044', sep='-')
	checkEquals(output[[2]]$'Description', Description)
	
	checkEquals(output[[3]]$'Index', 0.7266, tolerance=0.00016)
	Description=paste('At position:', '2034', sep='-')
	checkEquals(output[[3]]$'Description', Description)
}

sample.case2 <-function(){
	#Read the file
	Chlamydia <- read.fasta(file='../data/Chlamydia_1.txt')
	
	#STEP 1. Process the file
	Chlamydia <- processAllele(Chlamydia)
	
	#Since 1 sequence has deletion and is ignored
	checkEquals(length(Chlamydia), 55)
	
	#PERCENT MODE
	result=similar.percent(Chlamydia, 'A_D213')
	present=present.percent(result, 98, 100)
	
	#Should have no result because A_D213 is ignored
	checkEquals(result, NULL)
	checkEquals(present, list())
	
	#PERCENT MODE
	result=similar.percent(Chlamydia, 'H_S1432')
	present=present.percent(result, 98, 100)
	
	#Check result
	checkEquals(length(present), 100)
	checkEquals(present[[1]]$'position', 171)
	checkEquals(present[[1]]$'percent', 100)	
	
	#SIMPSON MODE
	result=branch.simpson(Chlamydia, level=1, numRes=3)
	output=present.simpson(Chlamydia, result)
	
	#Should have 3 results
	checkTrue(length(output)==3)
	
	checkEquals(output[[1]]$'Index',  0.7347, tolerance=0.00016)
	Description=paste('At position:', '1988', sep='-')
	checkEquals(output[[1]]$'Description', Description)
	
	checkEquals(output[[2]]$'Index', 0.7306, tolerance=0.00016)
	Description=paste('At position:', '2044', sep='-')
	checkEquals(output[[2]]$'Description', Description)
	
	checkEquals(output[[3]]$'Index', 0.7199, tolerance=0.00016)
	Description=paste('At position:', '2034', sep='-')
	checkEquals(output[[3]]$'Description', Description)
}

sample.case3 <-function(){
	#Read the file
	Chlamydia <- read.fasta(file='../data/Chlamydia_2.txt')
	
	#STEP 1. Process the file
	Chlamydia <- processAllele(Chlamydia)
	
	#Since 1 sequence has deletion and is ignored
	checkEquals(length(Chlamydia), 53)
	
	#PERCENT MODE
	result=similar.percent(Chlamydia, 'A_D213')
	present=present.percent(result, 98, 100)
	
	#Should have no result because A_D213 is ignored
	checkEquals(result, NULL)
	checkEquals(present, list())
	
	#PERCENT MODE
	result=similar.percent(Chlamydia, 'H_S1432')
	present=present.percent(result, 98, 100)
	
	#Check result
	checkEquals(length(present), 100)
	checkEquals(present[[1]]$'position', 171)
	checkEquals(present[[1]]$'percent', 100)	
	
	#SIMPSON MODE
	result=branch.simpson(Chlamydia, level=1, numRes=3)
	output=present.simpson(Chlamydia, result)
	
	#Should have 3 results
	checkTrue(length(output)==3)
	
	checkEquals(output[[1]]$'Index',  0.7343, tolerance=0.00016)
	Description=paste('At position:', '2044', sep='-')
	checkEquals(output[[1]]$'Description', Description)
	
	checkEquals(output[[2]]$'Index', 0.7329, tolerance=0.00016)
	Description=paste('At position:', '1988', sep='-')
	checkEquals(output[[2]]$'Description', Description)
	
	checkEquals(output[[3]]$'Index', 0.7263, tolerance=0.00016)
	Description=paste('At position:', '2034', sep='-')
	checkEquals(output[[3]]$'Description', Description)
}

test.deactivation <- function()
{
  DEACTIVATED('Deactivating integration test function')
}
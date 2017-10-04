sample.case1 <-function(){
	#Read the file
	Chlamydia <- read.fasta(file='../data/Chlamydia_mapped.txt')
	
	#STEP 1. Process the file
	Chlamydia <- processAllele(Chlamydia)
	untempered <-read.fasta(file='../data/Chlamydia_mapped.txt')
	
	#Since Chlamydia is normal
	checkIdentical(Chlamydia2, untempered)
	
	#PERCENT MODE
	result=similar.percent(Chlamydia, 'A_D213')
	present=present.percent(result, 98, 100)
	
	#All result should have percent higher or equal to 98
	for (a in 1:100){
	checkTrue(present[[a]]$percent>=98)
	}
	
	#SIMPSON MODE
	
}

sample.case2 <-function(){

}

test.deactivation <- function()
{
  DEACTIVATED('Deactivating integration test function')
}

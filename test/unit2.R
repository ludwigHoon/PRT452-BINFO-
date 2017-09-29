# Location of testing resources: ../resource/
# Resources used to test:
# 1. Chlamydia_mapped_Frankenstein_SNPs_Ns_to_Gs_oneofeachaussie.fas

#source('../minSNP/R/percent.R')
library('minSNP')
test.setUp <-function(){
Chlamydia<<- read.fasta(file='../resource/Chlamydia_mapped_Frankenstein_SNPs_Ns_to_Gs_oneofeachaussie.fas')
res<<-read.csv(file='../resource/result2.txt')
resp<<-list()
for (a in 1:length(res[[1]])){
resp[[a]]<<-list(position=res[[1]][a], percent=res[[2]][a])
}
}

test.similar.percent <- function()
{
	test.setUp()
	result=similar.percent(Chlamydia, 'ASD32')
	#Return NULL when the reference allele can't be found
	checkIdentical(result, NULL, "Can't find reference in allelic profiles")
	
	#Returning percentage at every SNP
	result=similar.percent(Chlamydia, 'A_D213')
	checkIdentical(class(result), 'list')
	checkIdentical(length(result), length(Chlamydia[[1]]))
	
	#Check the range of percentage (within 0 and 100)
	for (a in result){
		checkTrue((a$percent<=100)&&(a$percent>=0))
	}
	
	#Check the results of 352 SNPs (result from original minSNP software)
	resCheck=list.sort(result, (percent))[1:length(resp)]
	checkIdentical(resCheck, resp)
	
}

test.present.percent <- function()
{
	test.setUp()
	result=similar.percent(Chlamydia, 'A_D213')
	
	#Check number of result expected are correct
	checkEquals(length(present.percent(result, 98, 100)),100)
	checkEquals(length(present.percent(result, 98, 50)), 50)
	#Expect no result
	checkIdentical(present.percent(result, 101, 100), NULL, 'Percent error')
	
	present=present.percent(result, 98, 100)
	for (a in 1:100){
	checkTrue(present[[a]]$percent>=98)
	}
}

test.deactivation <- function()
{
  DEACTIVATED('Deactivating similar test function')
}
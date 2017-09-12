# Location of testing resources: ../resource/
# Resources used to test:
# 1. Chlamydia_mapped_Frankenstein_SNPs_Ns_to_Gs_oneofeachaussie.fas
# 2. Chlamydia_1.fas   --- (expecting 1 allele flagged - at position 1)
# 3. Chlamydia_2.fas   --- (expecting 5 allele flagged - at position 1, 3, 56)

source('../minSNP/R/dMode.R')

test.setUp <-function(){
Chlamydia<<- read.fasta(file='../resource/Chlamydia_mapped_Frankenstein_SNPs_Ns_to_Gs_oneofeachaussie.fas')
}

test.similar.simpson <- function()
{
	test.setUp()
	
	#Result from minSNP
	result1=similar.simpson(Chlamydia, 1) 
	checkEquals(result1[[1]]$position, 1988)
	checkEquals(result1[[1]]$index, 0.7344, tolerance=0.00016)
	
	#Result from minSNP
	result2=similar.simpson(Chlamydia, 2)
	checkEquals(result2[[1]]$position, 1988)
	checkEquals(result2[[1]]$index, 0.7344, tolerance=0.00016)
	checkEquals(result2[[2]]$position, 8241)
	checkEquals(result2[[2]]$index, 0.9025, tolerance=0.00016)
	
	#Predicted result, ensuring that in level 3, index is higher, 
	#and position is not the same as in level 1 and 2
	result3=similar.simpson(Chlamydia, 3)
	checkEquals(result3[[1]]$position, 1988)
	checkEquals(result3[[1]]$index, 0.7344, tolerance=0.00016)
	checkEquals(result3[[2]]$position, 8241)
	checkEquals(result3[[2]]$index, 0.9025, tolerance=0.00016)
	checkTrue(result3[[3]]$index>result3[[2]]$index)
	checkTrue(result3[[3]]$position!=result3[[1]]$position)
	checkTrue(result3[[3]]$position!=result3[[2]]$position)
	
	#Predicted result, ensuring that in level 4, index is higher, 
	#and position is not the same as in level 1, 2, and 3
	result4=similar.simpson(Chlamydia, 4)
	checkEquals(result4[[1]]$position, 1988)
	checkEquals(result4[[1]]$index, 0.7344, tolerance=0.00016)
	checkEquals(result4[[2]]$position, 8241)
	checkEquals(result4[[2]]$index, 0.9025, tolerance=0.00016)
	checkTrue(result4[[3]]$index>result4[[2]]$index)
	checkTrue(result4[[3]]$position!=result4[[1]]$position)
	checkTrue(result4[[3]]$position!=result4[[2]]$position)
	checkTrue(result4[[4]]$index>result4[[3]]$index)
	checkTrue(result4[[4]]$position!=result4[[1]]$position)
	checkTrue(result4[[4]]$position!=result4[[2]]$position)
	checkTrue(result4[[4]]$position!=result4[[3]]$position)
}


test.deactivation <- function()
{
  DEACTIVATED('Deactivating similar test function')
}
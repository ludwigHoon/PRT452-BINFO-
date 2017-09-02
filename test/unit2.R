# Location of testing resources: ../resource/
# Resources used to test:
# 1. Chlamydia_mapped_Frankenstein_SNPs_Ns_to_Gs_oneofeachaussie.fas
# 2. Chlamydia_1.fas   --- (expecting 1 allele flagged - at position 1)
# 3. Chlamydia_2.fas   --- (expecting 5 allele flagged - at position 1, 3, 56)
	

test.similar <- function()
{
	#Returning a list of SNPs (100) with percentage
	a <<- read.fasta(file='../resource/Chlamydia_mapped_Frankenstein_SNPs_Ns_to_Gs_oneofeachaussie.fas')
	checkIdentical(flagAllele(a), vector('character'))
	
	#Returning a list of SNPs (100) with percentage (Also have to take into account of included loci)
	includedLoci <- c(3,4,5)
	checkIdentical(flagAllele(a), c('A_D213'))
	
}

test.deactivation <- function()
{
  DEACTIVATED('Deactivating similar test function')
}
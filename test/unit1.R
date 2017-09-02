# Location of testing resources: ../resource/
# Resources used to test:
# 1. Chlamydia_mapped_Frankenstein_SNPs_Ns_to_Gs_oneofeachaussie.fas
# 2. Chlamydia_1.fas   --- (expecting 1 allele flagged - at position 1)
# 3. Chlamydia_2.fas   --- (expecting 5 allele flagged - at position 1, 3, 56)
	

test.flag <- function()
{
	#When the sequence within the test file are all of the same length
	a <<- read.fasta(file='../resource/Chlamydia_mapped_Frankenstein_SNPs_Ns_to_Gs_oneofeachaussie.fas')
	checkIdentical(flagAllele(a), vector('character'))
	
	#When there is a sequence within the test file that is not of the same length
	a <<- read.fasta(file='../resource/Chlamydia_1.fas')
	checkIdentical(flagAllele(a), c('A_D213'))
	
	#When there is multiple sequences within test file that are not of the same length
	a <<- read.fasta(file='../resource/Chlamydia_2.fas')
	checkIdentical(flagAllele(a), c('A_D213', 'Ia_SotoGIa3', 'D_SotoGD1'))	
}

test.deactivation <- function()
{
  DEACTIVATED('Deactivating this test function')
}
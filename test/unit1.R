# Location of testing resources: ../resource/
# Resources used to test:
# 1. Chlamydia_mapped_Frankenstein_SNPs_Ns_to_Gs_oneofeachaussie.fas
# 2. Chlamydia_1.fas   --- (expecting 1 allele flagged - at position 1)
# 3. Chlamydia_2.fas   --- (expecting 5 allele flagged - at position 1, 3, 56)
	
source('../minSNP/R/flag.R')

test.setUp <-function(){
Chlamydia <<- read.fasta(file='../resource/Chlamydia_mapped_Frankenstein_SNPs_Ns_to_Gs_oneofeachaussie.fas')
ErrorFile1 <<- read.fasta(file='../resource/Chlamydia_1.fas')
ErrorFile2 <<- read.fasta(file='../resource/Chlamydia_2.fas')
}

test.usualLength <- function()
{
	test.setUp()
	
	#Check that the function returns the normal length correctly for all 3 different files
	checkEqualsNumeric(usualLength(Chlamydia), 19570)
	checkEqualsNumeric(usualLength(ErrorFile1), 19570)
	checkEqualsNumeric(usualLength(ErrorFile2), 19570)
}

test.flag <- function()
{
	test.setUp()
	#When the sequence within the test file are all of the same length
	checkIdentical(flagAllele(Chlamydia), vector('character'))
	
	#When there is a sequence within the test file that is not of the same length
	checkIdentical(flagAllele(ErrorFile1), c('A_D213'))
	
	#When there is multiple sequences within test file that are not of the same length
	checkIdentical(flagAllele(ErrorFile2), c('A_D213', 'Ia_SotoGIa3', 'D_SotoGD1'))	
}

test.processAllele <- function()
{
	test.setUp()
	#Check that processed Allele returns correctly
	checkIdentical(processAllele(Chlamydia), Chlamydia)
	
	#Processed allelic profiles have removed the profiles with deletion
	#Only removed correct profiles
	checkIdentical(processAllele(ErrorFile1)[['A_D213']], NULL)
	checkEqualsNumeric(length(processAllele(ErrorFile1)), 55)
	checkIdentical(processAllele(ErrorFile2)[['A_D213']], NULL)
	checkIdentical(processAllele(ErrorFile2)[['Ia_SotoGIa3']], NULL)
	checkIdentical(processAllele(ErrorFile2)[['D_SotoGD1']], NULL)
	checkEqualsNumeric(length(processAllele(ErrorFile2)), 53)
}

test.deactivation <- function()
{
  DEACTIVATED('Deactivating flag test function')
}
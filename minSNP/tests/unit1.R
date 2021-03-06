# Location of testing resources: ../data/
# Resources used to test:
# 1. Chlamydia_mapped.txt
# 2. Chlamydia_1.txt   --- (expecting 1 allele flagged - at position 1)
# 3. Chlamydia_2.txt   --- (expecting 5 allele flagged - at position 1, 3, 56)
	
# source('../R/flag.R')

test.setUp <-function(){
Chlamydia <<- read.fasta(file='../data/Chlamydia_mapped.txt')
ErrorFile1 <<- read.fasta(file='../data/Chlamydia_1.txt')
ErrorFile2 <<- read.fasta(file='../data/Chlamydia_2.txt')
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

test.flagPosition <- function()
{
	test.setUp()

	#When the sequence within test file are all valid characters
	checkIdentical(flagPosition(Chlamydia), vector('numeric'))

	#When there are invalid characters
	checkIdentical(flagPosition(processAllele(ErrorFile1)), c(22,24))
	#When there are invalid characters, but dash treated as a type
	checkIdentical(flagPosition(processAllele(ErrorFile1), FALSE), c(22))

	#When there are invalid characters
	checkIdentical(flagPosition(processAllele(ErrorFile2)), c(1,2,3,4))
	#When there are invalid characters, but dash treated as a type
	checkIdentical(flagPosition(processAllele(ErrorFile2), FALSE), c(1,2))
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

test.tearDown <- function()
{
  #DEACTIVATED('Deactivating flag test function')
}
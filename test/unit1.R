test.flag <- function()
{
	a <- read.fasta(file='../resource/Chlamydia_mapped_Frankenstein_SNPs_Ns_to_Gs_oneofeachaussie.fas')
	checkIdentical(flagAllele(a), vector('numeric'))

}
 
test.deactivation <- function()
{
  DEACTIVATED('Deactivating this test function')
}
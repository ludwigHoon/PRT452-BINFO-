# Location of testing resources: ../data/
# Resources used to test:
# 1. Chlamydia_mapped.txt

#source('../R/dMode.R')

test.setUp <-function(){
Chlamydia<<- read.fasta(file='../data/Chlamydia_mapped.txt')
}

test.simpson.calculate <- function()
{
	#Check simpson's index with normal pattern
	checkEquals(simpson.calculate(c(3,4,5,6), 18), 7/9)
	checkEquals(simpson.calculate(list(a=3,b=4,c=5,d=6), 18), 7/9)
	#Check simpson's index at extreme 
	checkEquals(simpson.calculate(c(1,1,1), 3), 1)
	checkEquals(simpson.calculate(c(1), 1), 1)
}

test.simpson.pattern <- function()
{
	test.setUp()
	
	#Generate pattern from sequences when it's just 1 level
	pattern1=simpson.pattern(Chlamydia, 1)
	checkTrue('c' %in% names(pattern1))
	checkTrue('t' %in% names(pattern1))
	checkEquals(pattern1[['c']], 54)
	checkEquals(pattern1$'t', 2)
	
	#Appending pattern to result from 1st level to 2nd level
	appended=list()
	for(al in 1:length(Chlamydia)){
			appended[[getName(Chlamydia[[al]])]]=paste(appended[[getName(Chlamydia[[al]])]], Chlamydia[[al]][2], sep="")
			}
	pattern2=simpson.pattern(Chlamydia, 1, appended)
	checkTrue('gc' %in% names(pattern2))
	checkTrue('at' %in% names(pattern2))
	checkTrue('ac' %in% names(pattern2))
	checkEquals(pattern2[['gc']], 1)
	checkEquals(pattern2$'at', 2)
	checkEquals(pattern2$'ac', 53)
}

test.similar.simpson <- function()
{
	test.setUp()
	
	#Comparing simpson's with result from minSNP at position 1
	result1=similar.simpson(Chlamydia, 1) 
	checkEquals(result1[[1]]$position, 1988)
	checkEquals(result1[[1]]$index, 0.7344, tolerance=0.00016)
	
	#Comparing simpson's with result from minSNP at position 2
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

test.branch.simpson <- function()
{
	test.setUp()
	#Getting 3 results, 2nd result should ignore the first result at level 1 & 
	#3rd result should ignore both the 1st and 2nd at level 1
	result1=branch.simpson(Chlamydia, level=1, numRes=3)
	checkEquals(result1[[1]][[1]]$'position', similar.simpson(Chlamydia)[[1]]$'position')
	checkEquals(result1[[2]][[1]]$'position', similar.simpson(Chlamydia, excluded=result1[[1]][[1]]$'position')[[1]]$'position')
	checkEquals(result1[[3]][[1]]$'position', similar.simpson(Chlamydia, excluded=c(result1[[1]][[1]]$'position', result1[[2]][[1]]$'position'))[[1]]$'position')
	
	#Check that the pair of position returned by the 3 results are correct
	result2=branch.simpson(Chlamydia, level=2, numRes=3)
	reference=similar.simpson(Chlamydia, level=2)
	checkEquals(result2[[1]][[1]]$'position', reference[[1]]$'position')
	checkEquals(result2[[1]][[2]]$'position', reference[[2]]$'position')
	reference2=similar.simpson(Chlamydia, level=2, excluded=reference[[1]]$'position')
	checkEquals(result2[[2]][[1]]$'position', reference2[[1]]$'position')
	checkEquals(result2[[2]][[2]]$'position', reference2[[2]]$'position')
	reference3=similar.simpson(Chlamydia, level=2, excluded=c(reference[[1]]$'position',reference2[[1]]$'position'))
	checkEquals(result2[[3]][[1]]$'position', reference3[[1]]$'position')
	checkEquals(result2[[3]][[2]]$'position', reference3[[2]]$'position')
}

test.present.simpson <- function()
{
	test.setUp()
	
	#Confirm the number of result
	result=branch.simpson(Chlamydia, level=3, numRes=3)
	output=present.simpson(Chlamydia, result)
	checkTrue(length(output)==3)
	
	#Check the result at 1st level for the 3 outputs
	indDescription=paste('', round(result[[1]][[1]]$'index',4), round(result[[1]][[2]]$'index',4), round(result[[1]][[3]]$'index',4), sep=" ")
	checkEquals(output[[1]]$'Index', indDescription)
	Description=paste('At position:', result[[1]][[1]]$'position', result[[1]][[2]]$'position', result[[1]][[3]]$'position', sep='-')
	checkEquals(output[[1]]$'Description', Description)
	
	indDescription=paste('', round(result[[2]][[1]]$'index',4), round(result[[2]][[2]]$'index',4), round(result[[2]][[3]]$'index',4), sep=" ")
	checkEquals(output[[2]]$'Index', indDescription)
	Description=paste('At position:', result[[2]][[1]]$'position', result[[2]][[2]]$'position', result[[2]][[3]]$'position', sep='-')
	checkEquals(output[[2]]$'Description', Description)
	
	indDescription=paste('', round(result[[3]][[1]]$'index',4), round(result[[3]][[2]]$'index',4), round(result[[3]][[3]]$'index',4), sep=" ")
	checkEquals(output[[3]]$'Index', indDescription)
	Description=paste('At position:', result[[3]][[1]]$'position', result[[3]][[2]]$'position', result[[3]][[3]]$'position', sep='-')
	checkEquals(output[[3]]$'Description', Description)
}

test.tearDown <- function()
{
  #DEACTIVATED('Deactivating simpson test function')
}
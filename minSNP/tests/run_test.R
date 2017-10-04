library('RUnit')
library(seqinr)
library(rlist)

test.suite <- defineTestSuite('Testing',
							dirs = file.path('.'),
                            testFileRegexp = '^unit.+\\.R',
							testFuncRegexp = '^test.+')

test.result <- runTestSuite(test.suite)
printTextProtocol(test.result)
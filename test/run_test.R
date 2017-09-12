library('RUnit')
library('seqinr')

test.suite <- defineTestSuite('Testing',
							dirs = file.path('.'),
                            testFileRegexp = '^unit.+\\.R',
							testFuncRegexp = '^test.+')

test.result <- runTestSuite(test.suite)
printTextProtocol(test.result)
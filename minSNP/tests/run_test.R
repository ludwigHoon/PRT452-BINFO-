library('RUnit')

pkgname <- "minSNP"
require(pkgname, quietly=TRUE, character.only=TRUE) || stop("package '", pkgname, "' not found")

pattern <- '^unit.+\\.R' 
testFunctionRegexp = '^test.+'
dir <- '.'
suite <- defineTestSuite(name=paste(pkgname, "RUnit Tests"),
                         dirs=dir,
                         testFileRegexp=pattern,
                         testFuncRegexp = testFunctionRegexp,
                         rngKind="default",
                         rngNormalKind="default")
result <- runTestSuite(suite)
printTextProtocol(result)
printJUnitProtocol(result, fileName="junit.xml")

library('RUnit')
library('acml')

source("FittingFns.R") # Validated source for comparisons
source("FittingFnsC.R")
source("SimGenDatFns.R")

tolerance <- 1e-10

test.suite <- defineTestSuite("acml", dirs = file.path("."), testFileRegexp = '^\\d+\\.R')

test.result <- runTestSuite(test.suite)

printTextProtocol(test.result)

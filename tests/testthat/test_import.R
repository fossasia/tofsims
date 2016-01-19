context("Test importers.")
library(tofsimsData)
importFile<-system.file("rawdata", 
                        "trift_test_001.RAW", 
                        package = "tofsimsData")

test_that("read_BIF yields valid data for BIF6",{
    importFile<-system.file("rawdata", "512-10peak.bif6", package = "tofsimsData")
    instrument <- 'iontof'
    imported <- read_BIF(importFile, instrument)
    
    expect_that(is.integer(imported$imageData), is_true())
    expect_that(is.vector(imported$massValues), is_true())
    expect_that(length(imported$ids), equals(imported$nIntervals))
    expect_that(length(imported$middles), equals(length(imported$lowers)))
    expect_that(length(imported$lowers), equals(length(imported$uppers)))
})

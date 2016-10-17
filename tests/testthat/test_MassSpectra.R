context("Test MassSpectra.")
library(tofsimsData)
data("tofsimsData")

test_that("MassSpectra show method works with tofsims test data", {
    expect_that(show(testSpectra), prints_text("none"))
    expect_that(show(testSpectra), not(throws_error()))
    prComp(testSpectra)
    #expect_that(show(testSpectra), not(prints_text("none")))
})

test_that("MassSpectra plot method works with tofsims test data", {
    expect_that(plot(testSpectra), not(throws_error()))
    expect_that(plot(testSpectra, mzRange=c(1,10)), not(throws_error()))
    expect_that(plot(testSpectra, mzRange=c(-1,10)), not(throws_error()))
    expect_that(plot(testSpectra, mzRange=c(10,1)), not(throws_error()))
    expect_that(plot(testSpectra, normalize=TRUE), not(throws_error()))
    expect_that(plot(testSpectra, normalize=FALSE), not(throws_error()))
})

test_that("MassSpectra scale method works with tofsims test data", {
    expect_that(scale(testSpectra), throws_error('can not autoscale a single spectra'))
    expect_that(scale(testImage), not(throws_error()))
    expect_that(scale(testImage, center = FALSE), not(throws_error()))
    expect_that(scale(testImage, scale = FALSE), not(throws_error()))
})

test_that("MassSpectra poissonScaling method works with tofsims test data", {
    expect_that(poissonScaling(testSpectra), not(throws_error()))
    expect_that(poissonScaling(testImage), not(throws_error()))
})

test_that("MassSpectra reduceSpectrumResolution method works with tofsims test data", {
    expect_that(reduceSpectrumResolution(testSpectra), not(throws_error()))
    expect_that(reduceSpectrumResolution(testImage), not(throws_error()))
    expect_that(reduceSpectrumResolution(testSpectra, mode="remove"), not(throws_error()))
    expect_that(reduceSpectrumResolution(testSpectra, mode="keep"), not(throws_error()))
    expect_that(reduceSpectrumResolution(testSpectra, mode="temp"), throws_error("Unsupported mode."))
    expect_that(reduceSpectrumResolution(testSpectra, everyN = -1), throws_error("Invalid.*"))
})

test_that("MassSpectra points method works with tofsims test data", {
    expect_that(points(testSpectra), not(throws_error()))
    expect_that(points(testImage), not(throws_error()))
    expect_that(points(testImage, mzRange=c(-1,10)), not(throws_error()))
})

test_that("MassSpectra overlayPlot method works with tofsims test data", {
    expect_that(overlayPlot(testSpectra), throws_error())
    expect_that(overlayPlot(list(testSpectra, testImage)), not(throws_error()))
})

test_that("MassSpectra makeTIC method works with tofsims test data", {
    expect_that(makeTIC(testSpectra), not(throws_error()))
    expect_that(makeTIC(testImage), not(throws_error()))
})

test_that("MassSpectra smootherGolay method works with tofsisms test data", {
    expect_that(smootherGolay(testSpectra, p = 3, n = 5), not(throws_error()))
    expect_that(smootherGolay(testImage, p = 3, n = 5), not(throws_error()))
    # expect_that(smootherGolay(testImage, p = -1, n = 9), not(throws_error()))
    # expect_that(smootherGolay(testImage, p = 3, n = -1), not(throws_error()))
    # expect_that(smootherGolay(testImage, p = c(1,2), n = 5), not(throws_error()))
    # expect_that(smootherGolay(testImage, p = 3, n = c(1,2)), not(throws_error()))
})

test_that("MassSpectra smootherSpline method works with tofsisms test data", {
    expect_that(smootherSpline(testSpectra), not(throws_error()))
    expect_that(smootherSpline(testImage), not(throws_error()))
    expect_that(smootherSpline(testSpectra, stepsize = -1), throws_error())
    #expect_that(smootherSpline(testSpectra, stepsize = c(1,2)), throws_error())
    #expect_that(smootherSpline(testSpectra, spar = c(1,2)), throws_error())
})

test_that("MassSpectra peakPick method works with tofsisms test data", {
    expect_that(peakPick(testSpectra), not(throws_error()))
    expect_that(peakPick(testSpectra, span=100), not(throws_error()))
    #expect_that(peakPick(testSpectra, span=-1), throws_error())
    #expect_that(peakPick(testSpectra, span=c(1,2)), throws_error())
})

test_that("MassSpectra getTOFs method works with tofsisms test data", {
    expect_that(getTOFs(testSpectra), not(throws_error()))
    expect_that(getTOFs(testImage), not(throws_error()))
})

test_that("MassSpectra calibPointNew method works with tofsisms test data", {
    expect_that(calibPointNew(testSpectra, mz = 15, value = 15.01551), not(throws_error()))
    expect_that(calibPointNew(testImage, mz = 15, value = 15.01551), not(throws_error()))
    expect_that(calibPointNew(testSpectra, mz = -1, value = -1.01551), not(throws_error()))
    expect_that(calibPointNew(testSpectra, mz = 5000000, value = 5000000.01551), not(throws_error()))
    expect_that(calibPointNew(testSpectra, mz = c(1,2), value = c(1.01551, 2.01551)), not(throws_error()))
})










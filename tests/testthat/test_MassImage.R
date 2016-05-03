context("Test MassImage.")
library(tofsimsData)
data("tofsimsData")

test_that("MassImage show method works with tofsims test data", {
    expect_that(show(testImage), prints_text("none"))
    expect_that(show(testImage), not(throws_error()))
    prComp(testImage)
    #expect_that(show(testImage), not(prints_text("none")))
})

test_that("MassImage plot method works with tofsims test data", {
    expect_that(plot(testImage), not(throws_error()))
})

test_that("MassImage image method works with tofsims test data", {
    expect_that(image(testImage), not(throws_error()))
    expect_that(image(testImage, mzSelect = 1), not(throws_error()))
    expect_that(image(testImage, mzSelect = c(1,2)), not(throws_error()))
})

test_that("MassImage xySpec method works with tofsims test data", {
    expect_that(xySpec(testImage, x = c(1), y = c(1)), not(throws_error()))
    expect_that(validObject(xySpec(testImage, x = c(1), y = c(1))), is_true())
    expect_that(xySpec(testImage, x = c(1,2), y = c(1,2)), not(throws_error()))
    expect_that(validObject(xySpec(testImage, 
                                   x = c(1,2), 
                                   y = c(1,2))), is_true())
    expect_that(xySpec(testImage, x = c(1, 2), y = c(1)), throws_error())
    expect_that(xySpec(testImage, 
                       x = xdim(testImage) + 1, 
                       y = c(1)), throws_error())
    expect_that(xySpec(testImage, 
                       x = c(1), 
                       y = ydim(testImage) + 1), throws_error())
})

test_that("MassImage subset method works with tofsims test data", {
        expect_that(subset(testImage, 
                       xyUpperLeft = c(100, 100), 
                       xyLowerRight = c(50, 50)), not(throws_error()))
    expect_that(subset(testImage, 
                       xyUpperLeft = c(1, 1), 
                       xyLowerRight = c(50, 50)), not(throws_error()))
#    expect_that(subset(testImage, 
#                       xyUpperLeft = c(1, 1, 1), 
#                       xyLowerRight = c(50, 50)), throws_error())
})

test_that("MassImage imageMatrix method works with tofsims test data", {
    expect_that(is.matrix(imageMatrix(testImage)), is_true())
})

test_that("MassImage binning method works with tofsims test data", {
  library(BiocParallel)
  ## the following param will cause to run non parallel
  register(SerialParam(), default=TRUE)
    expect_that(binning(testImage, 
                        binningFactor = 2), not(throws_error()))
    expect_that(binning(testImage, 
                        binningFactor = -1), throws_error())
})





















context("Test PCA (PrComp & PrinComp).")
library(tofsimsData)
nulPrCompObj <- new("PrComp")
nulPrinCompObj <- new("PrinComp")
data("tofsimsData")
prCompObj <- prComp(testImage)
prinCompObj <- prinComp(testImage)
invalidPlotType <- "invalidPlotType"

test_that("PrComp/PrinComp yield valid PCA object", {
    expect_that(validObject(prCompObj), is_true())
    expect_that(validObject(prinCompObj), is_true())
})

test_that("image method behave correctly on PrComp / PrinComp object", {
    expect_error(image(nulPrCompObj, 1), "Invalid .*")
    expect_error(image(nulPrinCompObj, 1), "Invalid .*")
    #expect_error(image(prCompObj, nComp(prCompObj) + 1), "Invalid .*")
    #expect_that(image(prCompObj), throws_error())
    #expect_that(image(prCompObj, 1), not(throws_error()))
    #expect_error(image(prinCompObj, nComp(prinCompObj) + 1), "Invalid .*")
    #expect_that(image(prinCompObj), throws_error())
    #expect_that(image(prinCompObj, 1), not(throws_error()))
})

test_that("smoothScatter method behave correctly on PrComp / PrinComp object", {
    expect_error(smoothScatter(nulPrCompObj), "Invalid .*")
    expect_error(smoothScatter(nulPrinCompObj), "Invalid .*")
    #expect_error(smoothScatter(prCompObj, 
    #                           pcType=invalidPlotType), 
    #             "Invalid plot type.")
    #expect_error(smoothScatter(prCompObj, 
    #                           pcType=c(invalidPlotType, invalidPlotType)), 
    #             "Invalid plot type.")
    #expect_error(smoothScatter(prCompObj, 
    #                           comps = c(1,2,3)), 
    #             "Invalid number of components. It should be > 0 <= 2.")
    #expect_error(smoothScatter(prCompObj, 
    #                           comps = c()), 
    #             "Invalid number of components. It should be > 0 <= 2.")
    #expect_error(smoothScatter(prCompObj, 
    #                           comps = c(nComp(prCompObj) + 1)), 
    #             "Invalid component index.*")
    #expect_that(smoothScatter(prCompObj), not(throws_error()))
    
    #expect_error(smoothScatter(prinCompObj, 
    #                           pcType=invalidPlotType), 
    #             "Invalid plot type.")
    #expect_error(smoothScatter(prinCompObj, 
    #                           pcType=c(invalidPlotType, invalidPlotType)), 
    #             "Invalid plot type.")
    #expect_error(smoothScatter(prinCompObj, 
    #                           comps = c(1,2,3)), 
    #             "Invalid number of components. It should be > 0 <= 2.")
    #expect_error(smoothScatter(prinCompObj, 
    #                           comps = c()), 
    #             "Invalid number of components. It should be > 0 <= 2.")
    #expect_error(smoothScatter(prinCompObj, 
    #                           comps = c(nComp(prinCompObj) + 1)), 
    #             "Invalid component index.*")
    #expect_that(smoothScatter(prinCompObj), not(throws_error()))
})

test_that("plot method behave correctly on PrComp / PrinComp object", {
    expect_error(plot(nulPrCompObj), "Invalid .*")
    expect_error(plot(nulPrinCompObj), "Invalid .*")
    #expect_error(plot(prCompObj, 
    #                           pcType=invalidPlotType), 
    #             "Invalid plot type.")
    #expect_error(plot(prCompObj, 
    #                           pcType=c(invalidPlotType, invalidPlotType)), 
    #             "Invalid plot type.")
    #expect_error(plot(prCompObj, 
    #                           comps = c(1,2,3,4)), 
    #             "Invalid number of components. It should be > 0 <= 3.")
    #expect_error(plot(prCompObj, 
    #                           comps = c()), 
    #             "Invalid number of components. It should be > 0 <= 3.")
    #expect_error(plot(prCompObj, 
    #                           comps = c(nComp(prCompObj) + 1)), 
    #             "Invalid component index.*")
    expect_that(plot(prCompObj), not(throws_error()))
    
    #expect_error(plot(prinCompObj, 
    #                           pcType=invalidPlotType), 
    #             "Invalid plot type.")
    #expect_error(plot(prinCompObj, 
    #                           pcType=c(invalidPlotType, invalidPlotType)), 
    #             "Invalid plot type.")
    #expect_error(plot(prinCompObj, 
    #                           comps = c(1,2,3,4)), 
    #             "Invalid number of components. It should be > 0 <= 3.")
    #expect_error(plot(prinCompObj, 
    #                           comps = c()), 
    #             "Invalid number of components. It should be > 0 <= 3.")
    #expect_error(plot(prinCompObj, 
    #                           comps = c(nComp(prinCompObj) + 1)), 
    #             "Invalid component index.*")
    #expect_that(plot(prinCompObj), not(throws_error()))
})
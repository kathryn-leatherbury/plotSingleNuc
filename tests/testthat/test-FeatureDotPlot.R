test_that("FeatureDotPlot is exported", {
  pkg <- utils::packageName()
  expect_true("FeatureDotPlot" %in% getNamespaceExports(pkg))
})

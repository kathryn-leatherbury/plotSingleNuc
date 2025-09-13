test_that("FeaturePlotMultiGene is exported", {
  pkg <- utils::packageName()
  expect_true("FeaturePlotMultiGene" %in% getNamespaceExports(pkg))
})

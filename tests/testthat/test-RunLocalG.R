library(testthat)
library(semla)

# Mock Seurat object
se_mbrain <- readRDS(system.file("extdata/mousebrain", "se_mbrain", package = "semla"))

test_that("RunLocalG.Seurat returns valid results", {
  result <- RunLocalG(se_mbrain, features = VariableFeatures(se_mbrain)[1:2])
  
  # Check if the result is a Seurat object with results in metadata
  expect_s4_class(result, "Seurat")
  expect_true(all(c("Gi[Hbb-bs]", "Gi[Hbb-bs]") %in% names(result[[]])))

  result <- RunLocalG(se_mbrain, features = VariableFeatures(se_mbrain)[1:10], store_in_metadata = FALSE, alternative = "greater")
  expect_true("GiScores" %in% names(result@assays))
  expect_true(all(c("Gi[Hbb-bs]", "Gi[Hbb-bs]") %in% rownames(result@assays$GiScores)))
})

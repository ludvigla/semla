library(testthat)
library(semla)

# Load test data
se_mbrain <- readRDS(system.file("extdata/mousebrain", "se_mbrain", package = "semla"))
featureMat <- FetchData(se_mbrain, vars = VariableFeatures(se_mbrain)[1:100])

coordfile <- system.file("extdata/mousebrain/spatial", "tissue_positions_list.csv", package = "semla")
xys <- setNames(read.csv(coordfile, header = FALSE), nm = c("barcode", "selection", "grid_y", "grid_x", "y", "x"))
xys$sampleID <- 1
xys <- xys |>
  dplyr::filter(selection == 1) |>
  dplyr::select(barcode, x, y, sampleID) |>
  tibble::as_tibble()

spatnet <- GetSpatialNetwork(xys)

result <- CorSpatialFeatures(featureMat, spatnet, nCores = 1, verbose = FALSE)
result_seurat <- CorSpatialFeatures(se_mbrain, features = VariableFeatures(se_mbrain)[1:100], verbose = FALSE)

test_that("CorSpatialFeatures.default handles errors correctly", {
  
  # Test with non-matrix objects
  expect_error(CorSpatialFeatures.default(list(), list()), 
               "Invalid format of feature matrix: 'list'", 
               fixed = TRUE)
  
  # Test with matrix objects but wrong dimensions
  expect_error(CorSpatialFeatures.default(matrix(nrow=0, ncol=0), list()), 
               "Invalid dimensions of feature matrix: '0x0'", 
               fixed = TRUE)
  
  # Test with non-list spatnet objects
  expect_error(CorSpatialFeatures.default(matrix(nrow=1, ncol=1), "string"), 
               "Invalid format of spatnet: 'character'", 
               fixed = TRUE)
  
  # Test with incorrect match between feature matrix and spatnet
  expect_error(CorSpatialFeatures.default(matrix(nrow=1, ncol=1, dimnames=list("abc", "abc")), list(data.frame(from="xyz", to="xyz"))),
               "1 spots in the spatial networks could not be found in the feature matrix.",
               fixed = TRUE)
})


# Test for valid output type
test_that("CorSpatialFeatures returns a list of tibbles", {
  expect_type(result, "list")
  expect_s3_class(result[[1]], "tbl_df")
  expect_type(result_seurat, "list")
  expect_s3_class(result_seurat[[1]], "tbl_df")
})

# Test for valid output size
test_that("CorSpatialFeatures returns the correct number of rows", {
  expect_equal(nrow(result[[1]]), ncol(featureMat))
  expect_equal(nrow(result_seurat[[1]]), ncol(featureMat))
})

# Test for valid column names
test_that("CorSpatialFeatures returns the correct column names", {
  expect_identical(colnames(result[[1]]), c("gene", "cor"))
  expect_identical(colnames(result_seurat[[1]]), c("gene", "cor"))
})

result <- CorSpatialFeatures(featureMat, spatnet, across_all = TRUE, nCores = 1, verbose = FALSE)
result_seurat <- CorSpatialFeatures(se_mbrain, features = VariableFeatures(se_mbrain)[1:100], across_all = TRUE, nCores = 1, verbose = FALSE)

# Test for valid output type when across_all = TRUE
test_that("CorSpatialFeatures returns a tibble when across_all = TRUE", {
  expect_s3_class(result, "tbl_df")
  expect_s3_class(result_seurat, "tbl_df")
})

# Test for valid output size
test_that("CorSpatialFeatures returns the correct number of rows when across_all = TRUE", {
  expect_equal(nrow(result), ncol(featureMat))
  expect_equal(nrow(result_seurat), ncol(featureMat))
})

# Test when calculate_pvalue = TRUE
result <- CorSpatialFeatures(featureMat, spatnet, nCores = 1, calculate_pvalue = TRUE, verbose = FALSE)
result_seurat <- CorSpatialFeatures(se_mbrain, features = VariableFeatures(se_mbrain)[1:100], calculate_pvalue = TRUE, verbose = FALSE)

# Test for valid output type when calculate_pvalue = TRUE and across_all = FALSE
test_that("CorSpatialFeatures returns a list of tibbles when calculate_pvalue = TRUE", {
  expect_type(result, "list")
  expect_s3_class(result[[1]], "tbl_df")
  expect_type(result_seurat, "list")
  expect_s3_class(result_seurat[[1]], "tbl_df")
})

# Test for valid output size when calculate_pvalue = TRUE
test_that("CorSpatialFeatures returns the correct number of rows when calculate_pvalue = TRUE", {
  expect_equal(nrow(result[[1]]), ncol(featureMat))
  expect_equal(nrow(result_seurat[[1]]), ncol(featureMat))
})

# Test for valid column names
test_that("CorSpatialFeatures returns the correct column names when calculate_pvalue = TRUE", {
  expect_identical(colnames(result[[1]]), c("gene", "cor", "pval", "adj.pval"))
  expect_identical(colnames(result_seurat[[1]]), c("gene", "cor", "pval", "adj.pval"))
})

# Test when calculate_pvalue = TRUE and across_all = TRUE
result <- CorSpatialFeatures(featureMat, spatnet, nCores = 1, across_all = TRUE, calculate_pvalue = TRUE, verbose = FALSE)
result_seurat <- CorSpatialFeatures(se_mbrain, features = VariableFeatures(se_mbrain)[1:100], across_all = TRUE, calculate_pvalue = TRUE, verbose = FALSE)

# Test for valid output type when calculate_pvalue = TRUE and across_all = TRUE
test_that("CorSpatialFeatures returns a tibble when calculate_pvalue = TRUE and across_all = TRUE", {
  expect_s3_class(result, "tbl_df")
  expect_s3_class(result_seurat, "tbl_df")
})

# Test for valid column names
test_that("CorSpatialFeatures returns the correct column names when calculate_pvalue = TRUE and across_all = TRUE", {
  expect_identical(colnames(result), c("gene", "cor", "pval", "adj.pval"))
  expect_identical(colnames(result_seurat), c("gene", "cor", "pval", "adj.pval"))
})


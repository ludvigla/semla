library(testthat)
library(semla)

se_mcolon <- readRDS(system.file("extdata/mousecolon", "se_mcolon", package = "semla"))
se_mcolon <- SubsetSTData(se_mcolon, expression = selection == "GALT")
spatial_network <- GetSpatialNetwork(se_mcolon)[[1]]
nodes <- spatial_network |> select(from, x_start, y_start) |> group_by(from) |> 
  slice_head(n = 1) |> rename(name = from, x = x_start, y = y_start)
edges <- spatial_network |> select(from, to) |> mutate(keep = TRUE)
spatial_network <- tidygraph::tbl_graph(nodes = nodes, edges = edges, directed = FALSE)

# Test the AdjustTissuecoordinates function
test_that("AdjustTissuecoordinates works as expected", {
  
  res <- AdjustTissueCoordinates(full_graph = spatial_network)
  expect_true(inherits(res, what = "tbl_df"))
})

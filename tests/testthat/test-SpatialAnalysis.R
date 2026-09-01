test_that("findSpatiallyVariableMetabolites validates sampling controls", {
  expect_error(
    findSpatiallyVariableMetabolites(NULL, max_spots = 1),
    "max_spots must be NULL or one finite number >= 2",
    fixed = TRUE
  )
  expect_error(
    findSpatiallyVariableMetabolites(NULL, seed = Inf),
    "seed must be one finite number",
    fixed = TRUE
  )
})

test_that("findSpatiallyVariableMetabolites reads Seurat v5 assay layers", {
  function_body <- paste(deparse(body(findSpatiallyVariableMetabolites)), collapse = "\n")

  expect_match(function_body, "SeuratObject::LayerData", fixed = TRUE)
  expect_false(grepl("GetAssayData", function_body, fixed = TRUE))
})

test_that("add_frontlines_spe requires SpatialExperiment input", {
  skip_if_not_installed("SpatialExperiment")
  skip_if_not_installed("SummarizedExperiment")
  
  # Create mock border data
  border_df <- data.frame(
    spot_id = c("spot1", "spot2"),
    interface = c("1", "1"),
    frontline_mode = c("inner", "inner")
  )
  
  # Should fail with non-SpatialExperiment input
  expect_error(add_frontlines_spe(list(), border = border_df))
})

test_that("add_frontlines_spe requires border data", {
  skip_if_not_installed("SpatialExperiment")
  skip_if_not_installed("SummarizedExperiment")
  
  data("visium_simulated_spe", package = utils::packageName())
  
  # Should fail if border is NULL
  expect_error(add_frontlines_spe(visium_simulated_spe, border = NULL))
})

test_that("add_frontlines_spe adds border spots correctly", {
  skip_if_not_installed("SpatialExperiment")
  skip_if_not_installed("SummarizedExperiment")
  
  data("visium_simulated_spe", package = utils::packageName())
  
  # Create border data with spots from the SPE
  spot_ids <- colnames(visium_simulated_spe)
  border_df <- data.frame(
    spot_id = spot_ids[1:5],
    interface = rep("1", 5),
    frontline_mode = rep("inner", 5)
  )
  
  spe_result <- add_frontlines_spe(visium_simulated_spe, border = border_df)
  
  # Check that colData contains the new columns
  cd <- SummarizedExperiment::colData(spe_result)
  expect_true("is_frontline" %in% colnames(cd))
  expect_true("frontline_mode" %in% colnames(cd))
  expect_true("interface" %in% colnames(cd))
  
  # Check that border spots are marked correctly
  border_col <- cd$is_frontline[1:5]
  expect_true(all(border_col == "border"))
})

test_that("add_frontlines_spe adds core spots correctly", {
  skip_if_not_installed("SpatialExperiment")
  skip_if_not_installed("SummarizedExperiment")
  
  data("visium_simulated_spe", package = utils::packageName())
  
  spot_ids <- colnames(visium_simulated_spe)
  border_df <- data.frame(
    spot_id = spot_ids[1:3],
    interface = rep("1", 3),
    frontline_mode = rep("inner", 3)
  )
  
  core_df <- data.frame(
    spot_id = spot_ids[6:10],
    interface = rep("1", 5),
    frontline_mode = rep("inner", 5)
  )
  
  spe_result <- add_frontlines_spe(visium_simulated_spe, border = border_df, core = core_df)
  
  cd <- SummarizedExperiment::colData(spe_result)
  
  # Check that core spots are marked correctly
  core_col <- cd$is_frontline[6:10]
  expect_true(all(core_col == "core"))
  
  # Check that non-selected spots are NA
  expect_true(all(is.na(cd$is_frontline[11:50])))
})

test_that("add_frontlines_spe preserves interface and frontline_mode info", {
  skip_if_not_installed("SpatialExperiment")
  skip_if_not_installed("SummarizedExperiment")
  
  data("visium_simulated_spe", package = utils::packageName())
  
  spot_ids <- colnames(visium_simulated_spe)
  border_df <- data.frame(
    spot_id = spot_ids[1:5],
    interface = c("1", "2", "1", "2", "1"),
    frontline_mode = c("inner", "outer", "inner", "both", "inner")
  )
  
  spe_result <- add_frontlines_spe(visium_simulated_spe, border = border_df)
  
  cd <- SummarizedExperiment::colData(spe_result)
  
  # Check that interface and frontline_mode are preserved
  expect_equal(as.character(cd$interface[1:5]), 
               c("1", "2", "1", "2", "1"))
  expect_equal(as.character(cd$frontline_mode[1:5]), 
               c("inner", "outer", "inner", "both", "inner"))
})

test_that("add_frontlines_spe does not overwrite existing battlefield columns by default", {
  skip_if_not_installed("SpatialExperiment")
  skip_if_not_installed("SummarizedExperiment")
  
  data("visium_simulated_spe", package = utils::packageName())
  
  spot_ids <- colnames(visium_simulated_spe)
  
  # Add initial border data
  border_df_1 <- data.frame(
    spot_id = spot_ids[1:5],
    interface = rep("1", 5),
    frontline_mode = rep("inner", 5)
  )
  
  spe_result_1 <- add_frontlines_spe(visium_simulated_spe, border = border_df_1)
  
  cd_1 <- SummarizedExperiment::colData(spe_result_1)
  
  # Try to add different border data without erase -> function returns with warning/message
  # and doesn't add new data (sets border/core to NULL)
  border_df_2 <- data.frame(
    spot_id = spot_ids[11:15],
    interface = rep("2", 5),
    frontline_mode = rep("outer", 5)
  )
  
  spe_result_2 <- add_frontlines_spe(spe_result_1, border = border_df_2, erase = FALSE)
  
  cd <- SummarizedExperiment::colData(spe_result_2)
  
  # When erase=FALSE and columns exist, new data is not added 
  # and all battlefield columns are reset to NA
  expect_true(all(is.na(cd$is_frontline)))
  expect_true(all(is.na(cd$frontline_mode)))
  expect_true(all(is.na(cd$interface)))
})

test_that("add_frontlines_spe erases pre-existing columns when erase=TRUE", {
  skip_if_not_installed("SpatialExperiment")
  skip_if_not_installed("SummarizedExperiment")
  
  data("visium_simulated_spe", package = utils::packageName())
  
  spot_ids <- colnames(visium_simulated_spe)
  
  # Add initial border data
  border_df_1 <- data.frame(
    spot_id = spot_ids[1:5],
    interface = rep("1", 5),
    frontline_mode = rep("inner", 5)
  )
  
  spe_result_1 <- add_frontlines_spe(visium_simulated_spe, border = border_df_1)
  
  # Add different border data with erase=TRUE
  border_df_2 <- data.frame(
    spot_id = spot_ids[11:15],
    interface = rep("2", 5),
    frontline_mode = rep("outer", 5)
  )
  
  spe_result_2 <- add_frontlines_spe(spe_result_1, border = border_df_2, erase = TRUE)
  
  cd <- SummarizedExperiment::colData(spe_result_2)
  
  # First spots should now be NA (erased)
  expect_true(all(is.na(cd$is_frontline[1:5])))
  # New spots should be marked
  expect_true(all(cd$is_frontline[11:15] == "border"))
})

test_that("add_frontlines_spe handles border and core with same spots", {
  skip_if_not_installed("SpatialExperiment")
  skip_if_not_installed("SummarizedExperiment")
  
  data("visium_simulated_spe", package = utils::packageName())
  
  spot_ids <- colnames(visium_simulated_spe)
  
  # Create data with overlapping border/core - border should take precedence
  border_df <- data.frame(
    spot_id = spot_ids[1:5],
    interface = rep("1", 5),
    frontline_mode = rep("inner", 5)
  )
  
  core_df <- data.frame(
    spot_id = spot_ids[3:8],
    interface = rep("1", 6),
    frontline_mode = rep("inner", 6)
  )
  
  spe_result <- add_frontlines_spe(visium_simulated_spe, border = border_df, core = core_df)
  
  cd <- SummarizedExperiment::colData(spe_result)
  
  # Overlapping spots should be marked as border, not core
  expect_true(all(cd$is_frontline[3:5] == "border"))
  # Non-overlapping core spots should be core
  expect_true(all(cd$is_frontline[6:8] == "core"))
})

test_that("add_frontlines_spe returns a SpatialExperiment", {
  skip_if_not_installed("SpatialExperiment")
  skip_if_not_installed("SummarizedExperiment")
  
  data("visium_simulated_spe", package = utils::packageName())
  
  spot_ids <- colnames(visium_simulated_spe)
  border_df <- data.frame(
    spot_id = spot_ids[1:5],
    interface = rep("1", 5),
    frontline_mode = rep("inner", 5)
  )
  
  spe_result <- add_frontlines_spe(visium_simulated_spe, border = border_df)
  
  # Result should be a SpatialExperiment
  expect_s4_class(spe_result, "SpatialExperiment")
})

test_that("add_frontlines_spe handles empty border dataframe", {
  skip_if_not_installed("SpatialExperiment")
  skip_if_not_installed("SummarizedExperiment")
  
  data("visium_simulated_spe", package = utils::packageName())
  
  # Empty border dataframe with correct structure
  border_df <- data.frame(
    spot_id = character(0),
    interface = character(0),
    frontline_mode = character(0)
  )
  
  spe_result <- add_frontlines_spe(visium_simulated_spe, border = border_df)
  
  cd <- SummarizedExperiment::colData(spe_result)
  
  # All battlefield columns should be NA
  expect_true(all(is.na(cd$is_frontline)))
  expect_true(all(is.na(cd$frontline_mode)))
  expect_true(all(is.na(cd$interface)))
})

test_that("add_selections_to_spe handles border data without optional columns", {
  skip_if_not_installed("SpatialExperiment")
  skip_if_not_installed("SummarizedExperiment")
  
  data("visium_simulated_spe", package = utils::packageName())
  
  spot_ids <- colnames(visium_simulated_spe)
  
  # Border data with only spot_id column (missing interface and mode)
  border_df <- data.frame(
    spot_id = spot_ids[1:5]
  )
  
  spe_result <- add_selections_to_spe(visium_simulated_spe, border = border_df)
  
  cd <- SummarizedExperiment::colData(spe_result)
  
  # Border spots should be marked
  expect_true(all(cd$is_battlefield[1:5] == "border"))
  # interface and mode should be NA for these spots
  expect_true(all(is.na(cd$interface[1:5])))
  expect_true(all(is.na(cd$mode[1:5])))
})
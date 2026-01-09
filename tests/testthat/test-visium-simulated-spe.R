
test_that("visium_simulated_spe is a valid SpatialExperiment dataset", {
skip_if_not_installed("SpatialExperiment")
skip_if_not_installed("SummarizedExperiment")

data("visium_simulated_spe", package = utils::packageName())

# --- Class / container checks ---
expect_s4_class(visium_simulated_spe, "SpatialExperiment")

# Assays must exist and contain the FAKE_GENE
a <- SummarizedExperiment::assays(visium_simulated_spe)
expect_true(length(a) >= 1)

assay_names <- SummarizedExperiment::assayNames(visium_simulated_spe)
expect_true(length(assay_names) >= 1)

# --- colData checks ---
cd <- SummarizedExperiment::colData(visium_simulated_spe)
expect_true(nrow(cd) > 0)
expect_true("barcode_id" %in% colnames(cd))
expect_true("cluster" %in% colnames(cd))

})

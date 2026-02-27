#' Generate and visualize visiumHD_16um_simulated_spe dataset
#'
#' This script generates the hexagonal grid dataset for VisiumHD 16μm resolution
#' based on the exact original grid structure with 400 μm spacing.
#' Cluster assignment will be handled following the 8μm pattern.

library(SpatialExperiment)
library(dplyr)
library(ggplot2)
library(patchwork)

set.seed(42)

# Load original data for comparison
load('data/visiumHD_16um_simulated_spe.rda')

original_df <- data.frame(
  x = spatialCoords(visiumHD_16um_simulated_spe)[, 1],
  y = spatialCoords(visiumHD_16um_simulated_spe)[, 2],
  cluster = colData(visiumHD_16um_simulated_spe)$cluster,
  source = "Original"
)

# ============================================================================
# GENERATE GRID BASED ON ORIGINAL STRUCTURE (400 μm spacing, rectangular)
# ============================================================================

# Grid parameters for VisiumHD 16μm:
# - Regular spacing: 400 μm on both X and Y axes
# - X extent: 0-15600 (40 columns)
# - Y extent: 0-9600 (25 rows)
# - This creates a 40 × 25 = 1000 spot grid

spacing <- 400
x_max <- 15600
y_max <- 9600

# Create regular grid
hex_grid <- data.frame()

for (y in seq(0, y_max, spacing)) {
  for (x in seq(0, x_max, spacing)) {
    # Add small random jitter to mimic biological variation
    x_jit <- x + rnorm(1, 0, 5)
    y_jit <- y + rnorm(1, 0, 5)
    
    hex_grid <- rbind(hex_grid, data.frame(
      x = x_jit,
      y = y_jit
    ))
  }
}

# Verify we have exactly 1000 spots from 40 × 25 grid
stopifnot(nrow(hex_grid) == 1000)

# ============================================================================
# ASSIGN CLUSTERS (Same pattern as 8μm but scaled 2x)
# ============================================================================

generated_df <- hex_grid |>
  mutate(
    # Calculate distance from center for cluster 1
    dist_to_center = sqrt((x - 7800)^2 + (y - 4800)^2),
    
    cluster = case_when(
      x >= 1600 & x <= 3600 & y >= 6000 & y <= 8000 ~ 5,  # Cluster 5 (orange) - square
      x >= 9400 & y >= 4200 & y <= 5800 ~ 2,              # Cluster 2 (blue) - rectangle
      dist_to_center <= 2400 ~ 1,                          # Cluster 1 (red) - circular
      y > 5000 ~ 4,                                        # Cluster 4 (purple) - upper half
      y <= 5000 ~ 3                                        # Cluster 3 (green) - lower half
    ),
    cluster = factor(cluster),
    source = "Generated"
  ) |>
  select(-dist_to_center)

cat("Generated grid structure:\n")
cat("Total spots:", nrow(generated_df), "\n")
cat("X range: [", round(min(generated_df$x), 1), "-", round(max(generated_df$x), 1), "]\n", sep="")
cat("Y range: [", round(min(generated_df$y), 1), "-", round(max(generated_df$y), 1), "]\n")
cat("Grid dimensions: 40 columns × 25 rows (400 μm spacing)\n\n")

# ============================================================================
# CREATE COMPARISON PLOT
# ============================================================================

cluster_colors <- c("1" = "#E41A1C", "2" = "#377EB8", "3" = "#4DAF4A", 
                    "4" = "#984EA3", "5" = "#FF7F00")

p_orig <- ggplot(original_df, aes(x = x, y = y, color = cluster)) +
  geom_point(size = 2, alpha = 0.6) +
  scale_color_manual(values = cluster_colors) +
  coord_fixed() +
  theme_minimal() +
  labs(title = "ORIGINAL Dataset (VisiumHD 16μm)", x = "X (pixels)", y = "Y (pixels)") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 14))

p_gen <- ggplot(generated_df, aes(x = x, y = y, color = cluster)) +
  geom_point(size = 2, alpha = 0.6) +
  scale_color_manual(values = cluster_colors) +
  coord_fixed() +
  theme_minimal() +
  labs(title = "GENERATED Dataset (VisiumHD 16μm)", x = "X (pixels)", y = "Y (pixels)") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 14))

combined <- p_orig | p_gen

print(combined)

# Save plot
ggsave("inst/scripts/generated_visiumHD_16um_preview.png", combined, width = 14, height = 7, dpi = 100)
cat("\n✓ Plot saved to: inst/scripts/generated_visiumHD_16um_preview.png\n")

# ============================================================================
# SHOW SUMMARY STATISTICS
# ============================================================================

cat("\n=== SUMMARY ===\n")
cat("Original:\n")
print(table(original_df$cluster))
cat("\nGenerated:\n")
print(table(generated_df$cluster))
cat("\nX range: Original [", round(min(original_df$x), 1), "-", round(max(original_df$x), 1), 
    "] vs Generated [", round(min(generated_df$x), 1), "-", round(max(generated_df$x), 1), "]\n", sep="")
cat("Y range: Original [", round(min(original_df$y), 1), "-", round(max(original_df$y), 1),
    "] vs Generated [", round(min(generated_df$y), 1), "-", round(max(generated_df$y), 1), "]\n", sep="")

# ============================================================================
# CREATE SPATIALEXPERIMENT OBJECT
# ============================================================================

# Create spatial coordinates matrix (pxl_col_in_fullres, pxl_row_in_fullres)
spatial_coords <- as.matrix(generated_df[, c("x", "y")])
colnames(spatial_coords) <- c("pxl_col_in_fullres", "pxl_row_in_fullres")

# Create counts matrix (1 gene x 1000 spots)
n_spots <- nrow(generated_df)
counts_matrix <- matrix(
  sample(c(0, 1, 2, 3, 4, 5), size = n_spots, replace = TRUE, prob = c(0.6, 0.15, 0.1, 0.08, 0.05, 0.02)),
  nrow = 1,
  ncol = n_spots,
  dimnames = list(gene = "Gene1", spots = paste0("spot_", 1:n_spots))
)

# Create colData
col_data <- DataFrame(
  barcode_id = paste0("AAACCCAA-1_", formatC(1:n_spots, width = 5, flag = "0")),
  cluster = as.character(generated_df$cluster),
  sample_id = "visiumHD_16um_simulated"
)

# Create SpatialExperiment
visiumHD_16um_simulated_spe_generated <- SpatialExperiment(
  assays = list(counts = counts_matrix),
  spatialCoords = spatial_coords,
  colData = col_data
)

cat("\n✓ SpatialExperiment object created (NOT saved to /data)\n")
cat("  Dimensions:", dim(visiumHD_16um_simulated_spe_generated), "\n")
cat("  Use this object for testing, but /data/visiumHD_16um_simulated_spe.rda remains unchanged\n")

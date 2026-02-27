#' Generate and visualize visium_simulated_spe dataset
#'
#' This script generates the hexagonal grid dataset and immediately displays
#' a comparison plot with the original data for visual validation.

library(SpatialExperiment)
library(dplyr)
library(ggplot2)
library(patchwork)

set.seed(42)

# Load original data for comparison
load('data/visium_simulated_spe.rda')

original_df <- data.frame(
  x = spatialCoords(visium_simulated_spe)[, 1],
  y = spatialCoords(visium_simulated_spe)[, 2],
  cluster = colData(visium_simulated_spe)$cluster,
  source = "Original"
)

# ============================================================================
# GENERATE HEXAGONAL GRID DATASET
# ============================================================================

n_spots_target <- 1000
spot_spacing <- 60  # Reduced from 100 for denser points
vertical_spacing <- spot_spacing * sin(60 * pi/180)

# Create hexagonal grid directly
hex_grid <- data.frame()
spot_id <- 0

for (row in 0:30) {
  y_coord <- row * vertical_spacing
  if (y_coord > 1200) break
  
  x_offset <- ifelse(row %% 2 == 0, 0, spot_spacing/2)
  
  for (col in 0:40) {
    x_coord <- col * spot_spacing + x_offset
    if (x_coord > 2300) break
    
    # Smaller jitter for tighter clustering
    x_jit <- rnorm(1, 0, 1)
    y_jit <- rnorm(1, 0, 1)
    
    spot_id <- spot_id + 1
    
    hex_grid <- rbind(hex_grid, data.frame(
      x = x_coord + x_jit,
      y = y_coord + y_jit
    ))
    
    if (spot_id >= n_spots_target) break
  }
  
  if (spot_id >= n_spots_target) break
}

hex_grid <- hex_grid[1:n_spots_target, ]

# ============================================================================
# ASSIGN CLUSTERS TO HEXAGONAL GRID
# ============================================================================

# Cluster spatial assignments:
# C3 (green): top (y <= 523.9)
# C4 (purple): bottom (y > 523.9)
# C1 (red): circular region at center
# C2 (blue): horizontal lines from circle to right edge
# C5 (orange): square in middle-left of cluster 4

generated_df <- hex_grid |>
  mutate(
    # Calculate distance from center (for circular cluster 1)
    dist_to_center = sqrt((x - 1050)^2 + (y - 550)^2),
    
    cluster = case_when(
      # Cluster 5: square in middle-left of cluster 4 (orange)
      x >= 250 & x <= 600 & y >= 750 & y <= 1050 ~ 5,
      # Cluster 1: circular region at center (red) - priorité pour le cercle complet
      dist_to_center <= 260 ~ 1,
      # Cluster 2: horizontal lines from circle to right edge (blue)
      y > 475 & y < 650 & x > 1250 ~ 2,
      # Cluster 3: top region (green)
      y <= 523.9 ~ 3,
      # Cluster 4: bottom region (purple) - default
      TRUE ~ 4
    ),
    cluster = factor(cluster),
    source = "Generated"
  ) |>
  select(-dist_to_center)

cat("Generated cluster distribution:\n")
print(table(generated_df$cluster))

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
  labs(title = "ORIGINAL Dataset", x = "X (pixels)", y = "Y (pixels)") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 14))

p_gen <- ggplot(generated_df, aes(x = x, y = y, color = cluster)) +
  geom_point(size = 2, alpha = 0.6) +
  scale_color_manual(values = cluster_colors) +
  coord_fixed() +
  theme_minimal() +
  labs(title = "GENERATED Dataset", x = "X (pixels)", y = "Y (pixels)") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 14))

combined <- p_orig | p_gen

print(combined)

# Save plot
ggsave("inst/scripts/generated_dataset_preview.png", combined, width = 14, height = 7, dpi = 100)
cat("\n✓ Plot saved to: inst/scripts/generated_dataset_preview.png\n")

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
  sample_id = "visium_simulated"
)

# Create SpatialExperiment
visium_simulated_spe_generated <- SpatialExperiment(
  assays = list(counts = counts_matrix),
  spatialCoords = spatial_coords,
  colData = col_data
)

cat("\n✓ SpatialExperiment object created (NOT saved to /data)\n")
cat("  Dimensions:", dim(visium_simulated_spe_generated), "\n")
cat("  Use this object for testing, but /data/visium_simulated_spe.rda remains unchanged\n")

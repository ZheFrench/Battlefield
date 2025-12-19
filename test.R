#library(devtools)

#dir_package1 <- "/data2/villemin/Battlefield/"
#devtools::load_all(dir_package1)


# Packages
library(RANN)
library(SpatialExperiment)
library(SpatialExperimentIO)
library(readr)
library(ggplot2)

#devtools::install_github("thomasp85/ggforce")

file_tsv <- "visium_simulated_spots_1000.tsv"

# If you prefer TSV, comment the line above and use:
# df <- read_tsv(file_tsv, show_col_types = FALSE)

# --- Option A: clusters aléatoires (ex: 6 clusters) ---
set.seed(42)
k <- 6
df <- df %>%
  mutate(cluster = factor(sample(paste0("C", 1:k), n(), replace = TRUE)))

# --- Option B (alternative): clusters par bandes en X (plus "spatial")
# k <- 6
# breaks <- quantile(df$X, probs = seq(0, 1, length.out = k + 1))
# df <- df %>% mutate(cluster = cut(X, breaks = breaks, include.lowest = TRUE, labels = paste0("C", 1:k)))

# Plot
ggplot(df, aes(X, Y, color = cluster)) +
  geom_point(size = 1) +
  coord_equal() +
  theme_minimal(base_size = 12) +
  labs(title = "Simulated Visium spots colored by cluster", x = "X", y = "Y", color = "Cluster")
# Simulated Dataset Generation Scripts

This directory contains R scripts for generating reproducible simulated spatial transcriptomics datasets at different resolutions. These datasets are used throughout the Battlefield package for testing and demonstration.

## Files

### `generate_visium_simulated_spe.R`

**Purpose**: Generates the `visium_simulated_spe.rda` dataset (Visium 60μm resolution).

**What it does**:
- Creates a 1000-spot hexagonal grid with 60μm spot spacing
- Spatial extent: ~2173 × 1143 μm
- Assigns spots to 5 clusters with distinct spatial shapes:
  - **Cluster 1 (red)**: Central circular region
  - **Cluster 2 (blue)**: Horizontal rectangular barrier on the right
  - **Cluster 3 (green)**: Lower region
  - **Cluster 4 (purple)**: Upper region
  - **Cluster 5 (orange)**: Square patch in lower-left

**Target cluster distribution**:
- Cluster 1: ~100 spots
- Cluster 2: ~60-65 spots
- Cluster 3: ~430-440 spots
- Cluster 4: ~380-390 spots
- Cluster 5: ~20-25 spots

**Grid structure**: Regular hexagonal grid with small random jitter (±5 pixels)

---

### `generate_visiumHD_8um_simulated_spe.R`

**Purpose**: Generates the `visiumHD_8um_simulated_spe.rda` dataset (VisiumHD 8μm resolution).

**What it does**:
- Creates a 1000-spot rectangular grid with 200μm spacing on both axes
- Grid dimensions: 40 columns × 25 rows
- Spatial extent: 7800 × 4800 μm
- Assigns spots to 5 clusters using the same spatial strategy as Visium (scaled appropriately)

**Grid structure**: 
- Regular rectangular grid: (0, 200, 400, ..., 7800) × (0, 200, 400, ..., 4800)
- Small random jitter (±5 pixels) to mimic biological variation

**Cluster assignments**:
- **Cluster 1 (red)**: Circular region at center (x≈3900, y≈2400, radius≈1200)
- **Cluster 2 (blue)**: Rectangular barrier on the right (x≥4700, y∈[2100-2900])
- **Cluster 3 (green)**: Lower half (y≤2500)
- **Cluster 4 (purple)**: Upper half (y>2500)
- **Cluster 5 (orange)**: Square in left-middle of cluster 4 (x∈[800-1800], y∈[3000-4000])

---

### `generate_visiumHD_16um_simulated_spe.R`

**Purpose**: Generates the `visiumHD_16um_simulated_spe.rda` dataset (VisiumHD 16μm resolution).

**What it does**:
- Creates a 1000-spot rectangular grid with 400μm spacing on both axes (2× the 8μm dataset)
- Grid dimensions: 40 columns × 25 rows
- Spatial extent: 15600 × 9600 μm
- Assigns spots to 5 clusters using the same spatial strategy as VisiumHD 8μm (scaled 2×)

**Grid structure**:
- Regular rectangular grid: (0, 400, 800, ..., 15600) × (0, 400, 800, ..., 9600)
- Small random jitter (±5 pixels)

**Cluster assignments** (same structure as 8μm, scaled 2×):
- **Cluster 1 (red)**: Circular region at center (x≈7800, y≈4800, radius≈2400)
- **Cluster 2 (blue)**: Rectangular barrier on the right (x≥9400, y∈[4200-5800])
- **Cluster 3 (green)**: Lower half (y≤5000)
- **Cluster 4 (purple)**: Upper half (y>5000)
- **Cluster 5 (orange)**: Square in left-middle of cluster 4 (x∈[1600-3600], y∈[6000-8000])

---

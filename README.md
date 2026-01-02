# Battlefield <img  width="160" height="160" src="man/figures/logo.png" align="right" />
	
<br/>
<br/>

## Overview

**Battlefield**  is a Swiss-army toolkit originally developed to define and extract spatial spots from specific tissue regions—such as front regions, niche borders, invasive margins, and cluster interfaces—using spatial transcriptomics data or clustered tissue maps.

It has since been extended to support trajectory selection and layer inspection, and now provides a collection of low-level utilities for spatial transcriptomics analysis. These utilities are primarily intended to be reused within higher-level analytical packages.

Battlefield enables margin/trajectory/layer-aware spot selecIt is designed to work with sequencing-based platforms
such as Visium (classic 6.5 µm or 11 µm resolutions, and Visium HD).tion, allowing users to construct consistent comparison groups (e.g. front vs core or ordered spots along a trajectory).  

It is designed to work with sequencing-based platforms
such as Visium (classic 6.5 µm or 11 µm resolutions, and Visium HD).

For example, it can be integrated upstream to 
support insightful ligand–receptor analyses, 
using BulkSignalR, available from Bioconductor 
[here](https://www.bioconductor.org/packages/release/bioc/html/BulkSignalR.html).  

A number of visualization and data summary functions are proposed to
help defining spatial regions of interest on the top of previously defined 
clusters.

<!--#<img   src="man/figures/workflow.png" align="center" width="85%" height="85%" />-->
  

## Installation

``` R

# Battlefield directly from Bioconductor.
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("Battlefield")

# or development version via GitHub:
# install.packages("devtools")
devtools::install_github("ZheFrench/Battlefield",build_vignettes = TRUE)

# To read the vignette
# browseVignettes("Battlefield")

```

## Changes

Major changes are reported in the [NEWS file](https://github.com/ZheFrench/Battlefield/blob/master/NEWS), make sure to check it out if you want to follow the latest developments.

## Issues and bug reports
Please use https://github.com/ZheFrench/Battlefield/issues to submit issues, bug reports, and comments.


**Battlefield** has been successfully installed on Mac OS X, Linux, and Windows using R version 4.5.


The code in this repository is published with the [CeCILL](https://github.com/ZheFrench/Battlefield/blob/master/LICENSE.md) License.


<!-- badges: start -->
[![Generic badge](https://img.shields.io/badge/License-CeCILL-green.svg)](https://shields.io/)
<!-- badges: end -->




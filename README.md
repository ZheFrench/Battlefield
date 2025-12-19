# Battlefield <img  width="160" height="160" src="man/figures/logo.png" align="right" />
	
<br/>
<br/>

## Overview

**Battlefield**  is a Swiss-army toolkit to define and extract spatial 
"front" regions—niche borders, invasive margins, and cluster interfaces—from spatial transcriptomics or clustered tissue maps.  

It provides margin-aware spot selection to build consistent comparison groups (front vs core, front vs background, interface vs non-interface)
upstream of downstream analyses such as differential gene expression 
and cell–cell communication (ligand–receptor) inference.

For example, it can be integrated upstream to 
support insightful ligand–receptor analyses, 
using BulkSignalR, available from Bioconductor 
[here](https://www.bioconductor.org/packages/release/bioc/html/BulkSignalR.html).  

A number of visualization and data summary functions are proposed to
help navigating the predicted interactions.

<!--#<img   src="man/figures/workflow.png" align="center" width="85%" height="85%" />-->
  

## Installation

``` R

# Battlefield directly from Bioconductor.
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("Battlefield")

# or Installation goes via GitHub:
# install.packages("devtools")
devtools::install_github("ZheFrench/Battlefield",build_vignettes = TRUE)

# To read the vignette
# browseVignettes("Battlefield")

```

## Notes

For a version history/change logs, see the [NEWS file](https://github.com/ZheFrench/Battlefield/blob/master/NEWS).


**Battlefield** has been successfully installed on Mac OS X, Linux, and Windows using R version 4.5.


The code in this repository is published with the [CeCILL](https://github.com/ZheFrench/Battlefield/blob/master/LICENSE.md) License.


<!-- badges: start -->
[![Generic badge](https://img.shields.io/badge/License-CeCILL-green.svg)](https://shields.io/)
<!-- badges: end -->




# Palmtree
This package provides methods to interactively explore and visualize datasets with hierarchies. eg. single cells datasets with hierarchy over cells at different resolutions. 

For a walk through of various features of the app, checkout the vignette at 

https://hcbravolab.github.io/palmtree/Explore-PBMC.html

## Installation and requirements

`Palmtree` currently requires a development version of `epivizrData` which supports Sparse Measurements. We plan to push this to `bioc` soon. 

To install `Palmtree`, 

```{r}
BiocManager::install("epiviz/epivizrData")
BiocManager::install("HCBravoLab/palmtree")
```
***Note: you can also use devtools to install the package from github***

## Usage

To try out and explore the various features of the package, we provide the PBMC 3K dataset and the vignette `Explore-PBMC` walks you through the process. 

To checkout all the vignettes

```{r}
library(palmtree)
browseVignettes("palmtree")
```

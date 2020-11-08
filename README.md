# TreeViz
This package provides methods to interactively explore and visualize datasets with hierarchies. eg. single cells datasets with hierarchy over cells at different resolutions. 

## Installation and requiremenets

TreeViz currently requires a development version of `epivizrData` which supports Sparse Measurements. We plan to push this to `bioc` soon. 

To install TreeViz, 

```{r}
BiocManager::install("epiviz/epivizrData@SparseMeasurements")
BiocManager::install("HCBravoLab/TreeViz")
```
***Note: you can also use devtools to install the package from github***

## Usage

To try out and explore the various features of the package, we provide the PBMC 3K dataset and the vignette `Explore-PBMC` walks you through the process. 

To checkout all the vignettes

```{r}
library(TreeViz)
browseVignettes("TreeViz")
```
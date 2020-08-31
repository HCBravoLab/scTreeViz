# TreeViz
This package provides methods to interactively explore and visualize datasets with hierarchies. eg. microbiome datasets where features are mapped to a taxonomic hierarchy and single cells datasets with hierarchy over cells at different resolutions. 

## Installation and requiremenets

To install TreeViz, 

```{r}
BiocManager::install("HCBravoLab/TreeViz")
```

## Usage

To try out and explore the various features of the package, checkout the vignettes

```{r}
library(TreeViz)
browseVignettes("TreeViz")
```
## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----load-packages, message=FALSE, warning=FALSE------------------------------
library(palmtree)
library(metagenomeSeq)
library(msd16s)

## -----------------------------------------------------------------------------
data(mouseData)
counts <- MRcounts(mouseData)
hierarchy <- fData(mouseData)

## -----------------------------------------------------------------------------
tree <- TreeIndex(hierarchy)
tree

## -----------------------------------------------------------------------------
subset <- tree[1:1000,]
subset

## -----------------------------------------------------------------------------
tree_level <- splitAt(tree, selectedLevel = 2)
tree_level

## -----------------------------------------------------------------------------
tree_level <- splitAt(tree, selectedLevel = 2, format="TreeIndex")
tree_level


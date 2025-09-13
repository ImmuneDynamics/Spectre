# Spectre 1.3.0

* Added `utils.R` to keep internal reusable functions.
* Refactored `run.flowsom`:
    * Added `verbose` parameter to optionally print progress messages.
    * Overwrites `clust.name` and `meta.clust.name` columns in input data.table
    if they already exist.
    * Tidied up and simplified code.
* Updated `run.umap` to overwrite `umap.x.name` and `umap.y.name` columns in 
input data.table if they already exist.
* Updated `do.asinh`: 
  * Overwrite columns in input data.table if they already exist. 
  * Added support to specify a list of cofactors, so different values can be assigned to different markers.
* Refactored `make.colour.plot`:
  * Move some code out to internal functions for modularity.
  * Rewrote and simplified code.
  * Added more colour scheme options to support all options in viridis and RColorBrewer and jet.
  * Updated documentation.
  * Fixed issue #197 by using `geom_label_repel` to draw centroid labels. 
  * Added option to do fast plot using `scattermore`.
* Added deprecated message for `fast.colour.plot`.
* Fixed `package.check` for issue #186.
* Refactored `read.files` and remove any `setwd()` for issue #162.
* Add `add.label` to multi.plot for issue #202.
* Added a `NEWS.md` file to track changes to the package.
* Added `testthat` unittests.
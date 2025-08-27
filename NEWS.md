# Spectre 1.3.0

* Added `utils.R` to keep internal reusable functions.
* Refactored `run.flowsom`:
    * Added `verbose` parameter to optionally print progress messages.
    * Overwrites `clust.name` and `meta.clust.name` columns in input data.table
    if they already exist.
    * Tidied up and simplified code.
* Updated `run.umap` to overwrite `umap.x.name` and `umap.y.name` columns in 
input data.table if they already exist.
* Updated `do.asinh` to overwrite columns in input data.table if they already exist. Can now input string of cofactors, so different values can be assigned to different markers.
* Added a `NEWS.md` file to track changes to the package.
* Added `testthat` unittests.
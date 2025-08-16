# Spectre 1.3.0

* Added `utils.R` to keep internal reusable functions.
* Refactored `run.flowsom`:
    * Added `verbose` parameter to optionally print progress messages.
    * Overwrites `clust.name` and `meta.clust.name` columns in input data.table
    if they already exist.
    * Tidied up and simplified code.
* Added a `NEWS.md` file to track changes to the package.
* Added `testthat` unittests.
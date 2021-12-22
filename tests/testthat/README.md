# Unit tests

This folder contains all the unit tests for the Spectre package.

## Developer guideline

If you would like to contribute to the package, please make sure that all the existing unit tests are all passing.
This is to ensure that your changes has not broken anything.

In addition, if you are fixing a bug or introducing new features, please kindly update or add relevant unit tests.
It makes life easier for everyone! 

To avoid cluttering or overcrowding, please maintain (as much as possible) 1 unit test file per function, and name them sensibly.

## How to run all the tests

1. Install devtools and testthat.
2. Set your working directory to the package's root directory.
3. Run `devtools::test()`.

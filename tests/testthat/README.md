# Unit tests

This folder contains all the unit tests for the Spectre package.


## Developer guidelines

If you would like to contribute to the package, please make sure that all the existing unit tests are all passing.
This is the minimal step to ensure that your changes has not broken existing features.
If there are unit tests that are not relevant to your changes, but are failing, please raise an issue or if you are feeling inclined, fix them.

In addition, if you are fixing a bug or introducing new features, please kindly update or add and run the relevant unit tests.
It makes life easier for everyone! 

To avoid cluttering or overcrowding, please maintain (as much as possible) 1 unit test file per function, and name them sensibly.

## Preparation steps
These need to be done regardless of whether you are running only 1 unit test or all the unit tests:
1. Install devtools and testthat.
2. Set your working directory to the package's root directory.

Please do the above steps before proceeding.

## How to run all the unit tests

Run `devtools::test()`.

## How to run just one test file.
1. Open the test file.
2. Go to console and run `devtools::test_active_file()`.


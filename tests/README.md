# Test scripts

This folder contains all the unit and manual tests for the Spectre package.
Unit tests are stored in `testthat` folder while manual tests are stored in `manual_tests` folder.

**Please** do not change the `testthat` folder name as it is crucial for it to be named `testthat` for `devtools` to automatically run them.

Unit tests are well unit tests.
These are tests which validation and execution can be automated.
Manual tests are the ones which validation must be done manually, e.g. drawing colour plots requires manual inspection to see if certain features work.

## Developer guidelines

If you would like to contribute to the package, please make sure that at least all the existing unit tests are all passing.
This is the minimal step to ensure that your changes has not broken existing features.

In addition, if you are fixing a bug or introducing new features, please kindly update or add and run the relevant unit or manual tests.
It makes life easier for everyone! 
If you find any broken tests that are not relevant to your changes, either fix it (thank you!) or post an issue so someone can attend to it.

To avoid cluttering or overcrowding, please name your test file sensibly and provide a short description within the file as to what you are testing.

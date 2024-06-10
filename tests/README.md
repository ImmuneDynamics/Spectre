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
1. Install devtools, testthat, and usethis packages.
2. Set your working directory to the package's root directory.

Please do the above steps before proceeding.

## How to run all the unit tests

Run `devtools::test()`.

## How to run just one test file
1. Open the test file.
2. Go to console and run `devtools::test_active_file()`.

## How to create a new unit test file
Run `usethis::use_test(name='whatever_testname')`.
Change the `name` parameter to the name of the function you want to create test cases for. 
Avoid using `.` in the name. 
Replace it with `_` as required. 
For example, if the function is called `make.colour.plot`, the unit tests' name (fed to the `use_test` function) for the function should be `make_colour_plot`.

After running the function, there should be a new file in the `tests/testthat` folder with `test-<name>.R` filename.
For example, after running `usethis::use_test(name='make_colour_plot')`, there should be a new file called `test-make_colour_plot.R`.

## How to write a unit test

We use `testthat` package to mange unit tests. 

A unit test case can look something like the following:
```
test_that("make autograph works", {
    dat <- data.table::data.table(
        Samples = c("Mock_01", "Mock_02", "Mock_03", "WNV_01", "WNV_02", "WNV_03"), 
        Group = c(rep("Mock",3), rep("WNV", 3)), 
        Tcells = c(20, 40, 30, 60, 70, 80), 
        Bcells = c(90, 95, 70, 20, 15, 30), 
        Batch = c(1, 2, 1, 3, 2, 1))
    
    expect_no_error(
        make.autograph(
            dat = dat, 
            x.axis = "Group",
            y.axis = "Tcells", 
            colour.by = "Batch", 
            colours = c("Black", "Red","Blue"), 
            y.axis.label = "Proportion",
            save_to_disk = FALSE)
    )
    
})
```

Few things to note from the above:
* `test_that` function declares the unit test case.
* `make autograph works` is the description of the unit test case. It essentially describes what the unit test is testing.
* `expect_no_error` is a function in `testthat` package which make sure that a particular function (in this example `make.autograph`) can be run successfully (without producing any errors). There are various functions in `testthat` package you can use to test a variety of stuff, but `expect_no_error` is the simplest form of testing you can do to make sure the function can at least run properly.
* Ideally, each unit test only tests one aspect of the function. E.g., for `make.colour.plot`, if you want to test: (1) whether different colour schemes work, and (2) whether density plot can be produced, you should create 2 unit test cases (2 `test_that` functions), one for the colour scheme, the other for the density plot.

If you are unsure how to create unit tests, I highly recommend reading this [guide](https://r-pkgs.org/testing-basics.html) written by Hadley Wickham. 

# Contributing to Spectre

Feel inclined to contribute? Thank you! Your help is greatly appreciated!:tada:

So you have a fix/feature you want to add to Spectre, what do you do? 
Well, the easiest thing to do is fork the repository, make the changes in there, and create a pull request. 
One of maintainers will then review the pull request, then approve and merge the request if everything is well.

Now before you do that we just have a few guidelines that will make life much better for everyone if you follow them prior to submitting a pull request:
1. Make sure all the unit tests are passing. See the `tests/testthat` folder for the unit tests and how to run them. If you find tests that are failing but has nothing to do with your changes, raise an issue, or fix them (thanks!)
2. Fix up all the unit/manual tests relevant to your changes. If they don't exist, create one. This will certainly future proof further code changes.
3. Update the documentations. This is especially important. Read more below.

We're aware that there are functions which documentations need to be updated. 
We're constantly working on them.

## Updating documentations
We use Roxygen to generate documentations. 
If you have updated any, just run `devtools::document()`.
Devtools will automatically use Roxygen to take care of everything. 

Few caveats to ensure things are working properly:
1. Roxygen mixes up S3method and normal method if you use dots in your function name. To get around this, make sure you include the function name when you use `@export` tag. Something like `@export do.add.cols`.
2. In relation to the above, you also need to explicitly provide `@usage` tag. Otherwise Roxygen will, yet again, get it wrong.
3. If you are deprecating a function, simply add `.Deprecated(x)` in the beginning of the function. Fill x with the name of the function that replaces the deprecated function. Roxygen will take care of the rest (documentations, warning message, etc.).

Now, for the documentations, at the very least, please have:
1. Function name.
2. Short description of the function.
3. List of parameters and their description.
4. Usage (how to use it).
5. Return values (if it returns something).
6. Author.

It will be great if you can add examples as well. 

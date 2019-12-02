# Spectre
A computational toolkit for analysis of high-dimensional single-cell cytometry data. 

### Current version
v0.2.0

### About
The Spectre package is designed to allow for data analaysis workflows that simplify and streamline data manipulation and annotation, population identification (clustering, classificaiton), and dimensionality reduction (tSNE, UMAP) etc in high-dimensional cytometry data. Critically, the design of Spectre allows for a simple, clear, and modular design of analysis workflows, that can be utilised by data and laboratory scientists.

### How to use
In R, install and load the 'devtools' library.

```     
if(!require('devtools')) {install.packages('devtools')}
library('devtools')
```

Subsequently, use the 'install_github' function to install and load the Spectre package. By default this will load the 'master' branch, which is the same as the latest stable release version (listed at https://github.com/sydneycytometry/Spectre/releases). To install a specific release version, see https://cran.r-project.org/web/packages/githubinstall/vignettes/githubinstall.html.

```
if(!require('Spectre')) {install_github("sydneycytometry/spectre")}
```

You will see the following returned:
```
Downloading GitHub repo sydneycytometry/spectre@master
These packages have more recent versions available.
Which would you like to update?

 1: All                                 
 2: CRAN packages only                  
 3: None                                
 4: data.table (1.12.0 -> 1.12.2) [CRAN]
 ... etc
 ```
 
We suggest selecting 'none' (in this example, by entering '3' and pressing return) to avoid updating other packages.

If the package is sucessfully installed, you can load the library using:
```
library("Spectre")
```

Alternatively, you can go to releases (https://github.com/sydneycytometry/Spectre/releases) and download the latest stable release -- which can then be installed in R.

### Protocols and vignettes
Basic usage instructions, and vignettes are available from https://wiki.centenary.org.au/display/SPECTRE. More comprehensive protocols can be found by signing up to the Sydney Cytometry 'extranet': https://sydneycytometry.org.au/wiki-launch.

### Acknowledgements
Key contributors to the Spectre package are Thomas Ashhurst, Felix Marsh-Wakefield, and Givanna Putri. The Spectre package was constructed on the basis of the CAPX workflow in R (https://sydneycytometry.org.au/capx). Along with the various R packages used within Spectre, we would like to acknowledge the Seurat and cytofkit R packages from providing inspiration for elements of the package design.

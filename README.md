# Spectre
A computational toolkit for analysis of high-dimensional single-cell cytometry data. 

### How to use
You can install the package in R by using the following lines:
```     if(!require('devtools')) {install.packages('devtools')}
        library('devtools')

        if(!require('Spectre')) {install_github("sydneycytometry/spectre")}
        library("Spectre")
```

Alteranatively, you can go to releases (https://github.com/sydneycytometry/Spectre/releases) and download the latest stable release -- which can then be installed in R.

### About
The Spectre package is designed to allow for data analaysis workflows that simplify and streamline data manipulation and annotation, population identification (clustering, classificaiton), and dimensionality reduction (tSNE, UMAP) etc in high-dimensional cytometry data.

Critically, the design of Spectre allows for a simple, clear, and modular design of analysis workflows, that can be utilised by data and laboratory scientists.

### Protocols and vignettes
Protocols, usage instructions, and vignettes are available from https://wiki.centenary.org.au/display/SPECTRE.

### Acknowledgements
Key contributors to the Spectre package are Thomas Ashhurst, Felix Marsh-Wakefield, and Givanna Putri. The Spectre package was constructed on the basis of the CAPX workflow in R (https://sydneycytometry.org.au/capx). Along with the various R packages used within Spectre, we would like to acknowledge the Seurat and cytofkit R packages from providing inspiration for elements of the package design.

## Summary of components

### Data preparation


### Population identification


### Dimensionality reduction


### Plotting


### Summary Tables, heatmaps, and graphs

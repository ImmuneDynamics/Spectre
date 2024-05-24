# Spectre

A computational toolkit in R for the integration, exploration, and analysis of high-dimensional single-cell cytometry and imaging data.

<img src="https://raw.githubusercontent.com/tomashhurst/tomashhurst.github.io/master/images/Spectre.png"/>

**Current version**: [`v1.2.0`](https://github.com/ImmuneDynamics/Spectre/releases)

[![AppVeyor build status](https://ci.appveyor.com/api/projects/status/akhvb8wub6d6xhtd?svg=true)](https://ci.appveyor.com/project/tomashhurst/spectre)

<br/>

## About

Spectre is an R package that enables comprehensive end-to-end integration and analysis of high-dimensional cytometry data from different batches or experiments. Spectre streamlines the analytical stages of raw data pre-processing, batch alignment, data integration, clustering, dimensionality reduction, visualisation and population labelling, as well as quantitative and statistical analysis. To manage large cytometry datasets, Spectre was built on the data.table framework -- this simple table-like structure allows for fast and easy processing of large datasets in R. Critically, the design of Spectre allows for a simple, clear, and modular design of analysis workflows, that can be utilised by data and laboratory scientists. Recently we have extended the functionality of Spectre to support the analysis of Imaging Mass Cytometry (IMC) and scRNAseq data. For more information, please see our paper: [Ashhurst TM, Marsh-Wakefield F, Putri GH et al. (2021). Cytometry A. DOI: 10.1002/cyto.a.24350](https://doi.org/10.1002/cyto.a.24350).

Spectre was developed by [Thomas Ashhurst](https://immunedynamics.github.io/thomas-ashhurst/), [Felix Marsh-Wakefield](https://immunedynamics.github.io/felix-marsh-wakefield/), and [Givanna Putri](https://immunedynamics.github.io/givanna-putri/).

<br/>

## Citation

If you use Spectre in your work, please consider citing [Ashhurst TM, Marsh-Wakefield F, Putri GH et al. (2022). Cytometry A. DOI: 10.1002/cyto.a.24350](https://doi.org/10.1002/cyto.a.24350). To continue providing open-source tools such as Spectre, it helps us if we can demonstrate that our efforts are contributing to analysis efforts in the community. Please also consider citing the authors of the individual packages or tools (e.g. CytoNorm, FlowSOM, tSNE, UMAP, etc) that are critical elements of your analysis work.

<br/>

## Getting started

We recommend using Spectre with [R](https://cran.r-project.org/mirrors.html) and [RStudio](https://www.rstudio.com/products/rstudio/download/#download). If you are unfamiliar with using R and RStudio, check out our [R and Spectre basics guides](https://immunedynamics.io/spectre/install/#Basics_guide) for assistance. Once R and RStudio are installed, run the following to install the Spectre package.

```         
if(!require('remotes')) {install.packages('remotes')} # Installs the package 'remotes'
remotes::install_github(repo = "immunedynamics/spectre") # Install the Spectre package
```

For detailed installation instructions, and instructions for installing Spectre via Docker, see our [installation guide](https://immunedynamics.io/spectre/install/).

***In Spectre v1.1 and above we have removed the package dependencies `rgeos` and `rgdal` as these are no longer available on CRAN. The package should install fine without these dependencies, but some spatial functions may not work properly. If required, one can download the archived packages, unzip them, and then placed them in the R library location.***

<br />

## Workflows and protocols

Analysis workflows and protocols are available from <https://immunedynamics.github.io/spectre>.

<br/>

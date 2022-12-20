
<!-- README.md is generated from README.Rmd. Please edit that file -->

# BarleyEnvironmentGWS

## Description

This repository contains information and code for replicating the
analyses performed in the article below:

Article Title: *Accurate predictions of barley phenotypes using
genomewide markers and environmental covariates*  
Journal: *Crop Science*  
Authors: Jeffrey L. Neyhart, Kevin A.T. Silverstein, and Kevin P.
Smith  
Article doi:
<https://acsess.onlinelibrary.wiley.com/doi/10.1002/csc2.20782>

## Navigation

### Data

Data used in this study are available from the [Triticeae
Toolbox](https://triticeaetoolbox.org/barley). See [this
README](https://github.com/neyhartj/BarleyEnvironmentGWS/tree/master/Data)
for instructions on accessing this data.

### Code

Scripts used to complete the analysis and generate figures outlined in
the article above are available in the “Scripts” subfolder. See [this
README](https://github.com/neyhartj/BarleyEnvironmentGWS/tree/master/Scripts)
for information on the scripts and their intended execution order.

Three scripts in this directory are used by all other scripts:

1.  `source.R` - loads packages, creates directory links, and loads
    data.
2.  `source_MSI.R` - runs script \#1 with modifications for
    high-performance computing supported by the [Minnesota
    Supercomputing Institute](https://www.msi.umn.edu/).
3.  `source_functions.R` - loads additional functions into the
    environment.

The `figures.R` script in the “Scripts” subfolder produces the figures
found in the paper.

## Software/package versions

*R* version 4.0.2 was used to complete analyses in this project.

The following packages (with version information) were used to complete
this project:

| package   | version |
|:----------|:--------|
| APSIM     | 0.9.3   |
| broom     | 1.0.1   |
| car       | 3.0-10  |
| cowplot   | 1.1.1   |
| daymetr   | 1.5     |
| dbplyr    | 2.1.1   |
| dplyr     | 1.0.10  |
| ggdendro  | 0.1.22  |
| ggrepel   | 0.9.1   |
| ggsn      | 0.5.0   |
| glmnet    | 4.1-1   |
| lme4      | 1.1-27  |
| lmerTest  | 3.1-3   |
| lubridate | 1.8.0   |
| maps      | 3.3.0   |
| maptools  | 1.1-1   |
| modelr    | 0.1.8   |
| nasapower | 1.1.3   |
| paletteer | 1.3.0   |
| parallel  | 4.0.2   |
| patchwork | 1.1.1   |
| purrr     | 0.3.4   |
| raster    | 3.4-10  |
| readr     | 2.1.3   |
| rgeos     | 0.5-5   |
| rhwsd     | 1.0     |
| RSQLite   | 2.2.7   |
| rvest     | 1.0.3   |
| sp        | 1.4-5   |
| spdep     | 1.1-8   |
| stringr   | 1.4.0   |
| tibble    | 3.1.8   |
| tidyr     | 1.2.0   |

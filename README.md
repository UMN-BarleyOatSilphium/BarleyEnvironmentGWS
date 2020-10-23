
<!-- README.md is generated from README.Rmd. Please edit that file -->

# BarleyEnvironmentGWS

## Description

This repository contains information and code for replicating the
analyses performed in the article below:

Article Title: *Accurate genomewide predictions of barley phenotypes
across North American environments*  
Journal: *Science* (under review)  
Authors: Jeffrey L. Neyhart, Kevin A.T. Silverstein, and Kevin P.
Smith  
Link to article

## Navigation

### Data

Data used in this study are available from the [Triticeae
Toolbox](https://triticeaetoolbox.org/barley). See [this
README](https://github.com/neyhartj/S2MET_Predictions_Models/tree/master/Data)
for instructions on accessing this data.

### Code

Scripts used to complete the analysis and generate figures outlined in
the article above are available in the “Scripts” subfolder. See [this
README](https://github.com/neyhartj/S2MET_Predictions_Models/tree/master/Scripts)
for information on the scripts and their intended execution order.

Three scripts in this directory are used by all other scripts:

1.  `source.R` - loads packages, creates directory links, and loads
    data.
2.  `source_MSI.R` - runs script \#1 with modifications for
    high-performance computing supported by the [Minnesota
    Supercomputing Institute](https://www.msi.umn.edu/).
3.  `source_functions.R` - loads additional functions into the
    environment.

## Software/package versions

*R* version 3.5.3 was used to complete analyses in this project.

The following packages (with version information) were used to complete
this project:

| package   | version |
| :-------- | :------ |
| APSIM     | 0.9.2   |
| apsimr    | 1.2     |
| broom     | 0.5.2   |
| car       | 3.0-2   |
| cowplot   | 1.0.0   |
| daymetr   | 1.4     |
| dbplyr    | 1.3.0   |
| dplyr     | 0.8.0.1 |
| ggdendro  | 0.1-20  |
| ggrepel   | 0.8.0   |
| ggsn      | 0.5.0   |
| glmnet    | 2.0-18  |
| lme4      | 1.1-21  |
| lmerTest  | 3.1-0   |
| lubridate | 1.7.4   |
| maps      | 3.3.0   |
| maptools  | 0.9-5   |
| modelr    | 0.1.4   |
| nasapower | 1.1.0   |
| paletteer | 0.2.1   |
| parallel  | 3.5.3   |
| patchwork | 1.0.1   |
| purrr     | 0.3.4   |
| raster    | 2.8-19  |
| readr     | 1.3.1   |
| rgeos     | 0.4-2   |
| rhwsd     | 1.0     |
| RSQLite   | 2.1.1   |
| rvest     | 0.3.3   |
| sp        | 1.3-1   |
| spdep     | 1.1-2   |
| stringr   | 1.4.0   |
| tibble    | 2.1.1   |
| tidyr     | 0.8.3   |

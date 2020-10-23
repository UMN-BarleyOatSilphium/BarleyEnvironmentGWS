
<!-- README.md is generated from README.Rmd. Please edit that file -->

# BarleyEnvironmentGWS - Data

Below is information on the data in this subfolder, as well as
instructions for obtaining from the [Triticeae Toolbox
(T3)](https://triticeaetoolbox.org/barley) the experimental data used in
this study. These instructions valid as of 23 Oct 2020.

## Data in this subfolder

1.  `project_entries.xlsx` - metadata for genotypes used in this study,
    including breeding program name, population group, predigree,
    family, and other notes.
2.  `barley_base_crop_simulation.apsim` - a base APSIM simulation file
    that is modified to run the crop models described in this analysis.
3.  `barley_base_crop_simulation_factorial_cultivars.apsim` - a base
    APSIM simulation file that is modified to run the comparison of
    reference cultivars described in this analysis. .
4.  `trial_metadata.csv` - metadata for the experimental trials used in
    this study.

## Data from T3

### Genomewide marker data

The genomewide marker data used in this study are available from T3 and
(in a more readily usable form) from the Data Repository for the
University of Minnesota using the following persistent link:
<http://hdl.handle.net/11299/204785>.

### Phenotype data

1.  Go to <https://triticeaetoolbox.org/barley>.
2.  Under the “Select” tab, go to “Wizard (Lines, Traits, Trials)”
3.  Change the first drop-down box to “Experiment”
4.  Select the experiment “UMN Spring 2-row MET”
5.  Select all Years (2015-2017)
6.  Select all Trials
7.  Select the Traits: “grain yield”, “plant height”, “heading date”,
    “test weight”, and “grain protein”
8.  Select all Lines
9.  Click the “Save current selection” button.
10. Under the “Download” tab, go to “Genotype and Phenotype Data”.
11. Make sure only the “Phenotype” box is checked.
12. Click the “Create file” button with instructions for one column for
    each trait (not used by TASSEL).
13. Click the “Download Zip file of results” button.

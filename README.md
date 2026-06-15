# aSyn_spread

Code to reproduce network analyses and pathology progression plots in Vatsa, Brynildsen, Goralski, Kurgat et al. "Network analysis of α-synuclein pathology progression reveals p21-activated kinases as regulators of vulnerability". 

## Requirements:
  - Tested with R version 4.4.1 (2024-06-14)
    - Requisite packages are listed in `code/misc/packages.R`
   
## Data availability

The `data/` directory contains all pathology datasets required to reproduce the analyses in this manuscript. The PANGEA gene expression dataset is too large to distribute through this GitHub repository and should be downloaded separately from its Zenodo archive:

### PANGEA data: https://zenodo.org/records/20671427

After downloading, place the Pangea_data2.csv file in the `data/` directory. The code assumes that all pathology datasets and the PANGEA .csv are located within this directory.

## Tutorial

For a general overview of linear diffusion modeling in neurodegenerative disease, see tutorial in `example` folder here: https://github.com/ejcorn/tau-spread.

## Input specification

The file `aSyn_pipeline.R` is located in the main directory. This file will coordinate the sequential execution of all scripts within the code/ folder, generating all the figures in the paper and more. Custom specification of the following inputs at the top of `aSyn_pipeline.R` is required:
  - `basedir`:  path to the main directory containing the `code` and `data` folders 
  - `injection.site`: vector of character names of brain regions to inject pathology into
  - `measures`: character vector containing the names of the datasets being anayzed. For our dataset, these were measures of total, cell body, and neurite pathology.
  - `opdir`: here you can add some extra custom label for your output folder given a particular configuration at the top of the script

## Questions, suggestions, comments?

Please contact Kate Brynildsen (jbryn ~`AT`~ seas.upenn.edu) with any questions regarding network analysis and code, and contact Mike Henderson (Michael ~`DOT`~ Henderson ~`AT`~ vai.org) with any questions regarding experiments and data.

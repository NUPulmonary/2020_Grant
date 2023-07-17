
# Code for single-cell RNA-sequencing analysis

R virtual environment for the project was managed with `renv`, 
full list of dependencies and versions recorded in [`renv.lock`](renv.lock).
Python virtual environment was managed with `pipenv`, full list of dependencies with versions
recorded in [`Pipfile.lock`](Pipfile.lock).

This project uses [snakemake](https://snakemake.readthedocs.io/en/stable/) to generate data files.

### Contents:

* `01human-sars-ref`: snakefile to generate hybrid reference genome for 10x cellranger processing
* `02preprint`: code of preprint analysis
* `03raw-data`: snakefile to process all data from Illumina sequencing results to cellranger counts  
    and genotype demultiplexing with [`souporcell`](https://github.com/wheaton5/souporcell)
* `04dataset`: integration, analysis and plots of the data
* `05il6`: analysis of _IL6_ expression in the data
* `lib-r`: supplementary functions in R
* `lib`: supplementary functions in python

### Errata:

* July 17, 2023: metadata field `Patient` in `main.h5ad` and `supplement.h5ad` files deposited at  
    [GSE155249](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE155249) wrongly assigns all  
    cells in `Sample 16` to `Patient A`. Correct value of `Patient` field for cells from `Sample 16`  
    is deposited in [`2023-07-17-sample-16-patient.csv`](2023-07-17-sample-16-patient.csv).  
    GEO submission is being updated.


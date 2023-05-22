# README
## Article Information
This repository provides access to the data and source code used for the manuscript  
**Extinct and extant termites reveal the fidelity of behavior fossilization in amber**  
Nobuaki Mizumoto, Simon Hellemans, Michael S Engel, Thomas Bourguignon, Aleš Buček  

Preprint is available at bioRxiv.  
The all data will be uploaded in ZENODO upon acceptance:[![DOI](https://zenodo.org/badge/DOI/XXXDOIXXX.svg)](https://doi.org/XXXDOIXXX)

## Table of Contents
* [README](./README.md) - this file
* [scripts](./scripts)
  * [Preprocess.R](./scripts/Preprocess.R) - run before Output.R
  * [Output.R](./scripts/Output.R) - output all results (figures, statistics, tables)
  * [data_conversion.py](./scripts/data_conversion.py) - convert data from SLEAP for R (use before Preprocess.R)
* [img](./img) - all outputs are stored
* [data](./data)
  * [Cf-control](./data/Cf-control) - Coordinates generated by UMATracker (control)
  * [Cf-control-sleap](./data/Cf-control-sleap) - Coordinates generated by SLEAP (control)
  * [Cf-sticky](./data/Cf-sticky) - Coordinates generated by UMATracker (sticky)
  * [Cf-sticky-DLC](./data/Cf-sticky-DLC) - Coordinates generated by DLC (sticky)
  * [Fossil](./data/Fossil)
  * [bodysize_measurement](./data/bodysize_measurement)
  

## Session information
```
R version 4.0.1 (2020-06-06)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 19045)

Matrix products: default

locale:
[1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252    LC_MONETARY=English_United States.1252
[4] LC_NUMERIC=C                           LC_TIME=English_United States.1252    
system code page: 932

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] data.table_1.14.8

loaded via a namespace (and not attached):
 [1] Rcpp_1.0.10     lattice_0.20-41 here_1.0.1      png_0.1-8       withr_2.5.0     rprojroot_2.0.3 grid_4.0.1      jsonlite_1.8.4  Matrix_1.2-18  
[10] reticulate_1.28 tools_4.0.1     compiler_4.0.1 
```
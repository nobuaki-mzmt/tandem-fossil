# README
## Article Information
This repository provides access to the data and source code used for the manuscript  
### **Extinct and extant termites reveal the fidelity of behavior fossilization in amber**  
Nobuaki Mizumoto, Simon Hellemans, Michael S Engel, Thomas Bourguignon, Aleš Buček  

Preprint is available at bioRxiv. [![DOI:10.1101/2023.05.22.541647](http://img.shields.io/badge/DOI-10.1101/2023.05.22.541647-B31B1B.svg)](https://doi.org/10.1101/2023.05.22.541647)  
Zenodo link: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10557251.svg)](https://doi.org/10.5281/zenodo.10557251)

## Table of Contents
* [README](./README.md)
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
  * [Fossil](./data/Fossil) - data extracted from amber CT-scan
  * [bodysize_measurement](./data/bodysize_measurement) - body length measurement
  

## Session information
```
R version 4.3.1 (2023-06-16 ucrt)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 19045)

Matrix products: default


locale:
[1] LC_COLLATE=English_United States.utf8  LC_CTYPE=English_United States.utf8    LC_MONETARY=English_United States.utf8
[4] LC_NUMERIC=C                           LC_TIME=English_United States.utf8    

time zone: America/Chicago
tzcode source: internal

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] reticulate_1.30       dlcpr_0.1.0           stringr_1.5.0         data.table_1.14.8     extrafont_0.19       
 [6] car_3.1-2             carData_3.0-5         coxme_2.2-18.1        bdsmatrix_1.3-6       survival_3.5-5       
[11] survminer_0.4.9       ggpubr_0.6.0          exactRankTests_0.8-35 MASS_7.3-60           viridis_0.6.2        
[16] viridisLite_0.4.1     ggplot2_3.4.4        

loaded via a namespace (and not attached):
 [1] gtable_0.3.3      xfun_0.39         rstatix_0.7.2     lattice_0.21-8    vctrs_0.6.4       tools_4.3.1      
 [7] generics_0.1.3    tibble_3.2.1      fansi_1.0.4       pkgconfig_2.0.3   Matrix_1.5-4.1    lifecycle_1.0.3  
[13] compiler_4.3.1    farver_2.1.1      textshaping_0.3.6 munsell_0.5.0     Rttf2pt1_1.3.12   pillar_1.9.0     
[19] crayon_1.5.2      extrafontdb_1.0   tidyr_1.3.0       abind_1.4-5       nlme_3.1-162      km.ci_0.5-6      
[25] tidyselect_1.2.0  stringi_1.7.12    dplyr_1.1.4       purrr_1.0.1       labeling_0.4.2    splines_4.3.1    
[31] grid_4.3.1        colorspace_2.1-0  cli_3.6.1         magrittr_2.0.3    utf8_1.2.3        broom_1.0.4      
[37] withr_2.5.0       scales_1.2.1      backports_1.4.1   gridExtra_2.3     ggsignif_0.6.4    png_0.1-8        
[43] ragg_1.2.5        zoo_1.8-12        knitr_1.42        KMsurv_0.1-5      survMisc_0.5.6    rlang_1.1.0      
[49] Rcpp_1.0.10       isoband_0.2.7     xtable_1.8-4      glue_1.6.2        jsonlite_1.8.4    rstudioapi_0.14  
[55] R6_2.5.1          systemfonts_1.0.4
```

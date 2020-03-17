# R package to access DoRothEA's regulons

### Installation

```r
# install the development version from GitHub
# install.packages("devtools")
devtools::install_github("christianholland/dorothea")
```

### Session Info
```r
sessioninfo::session_info()
```

```r
─ Session info ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
 setting  value                       
 version  R version 3.6.2 (2019-12-12)
 os       macOS Mojave 10.14.5        
 system   x86_64, darwin15.6.0        
 ui       RStudio                     
 language (EN)                        
 collate  en_US.UTF-8                 
 ctype    en_US.UTF-8                 
 tz       Europe/Berlin               
 date     2020-03-17                  

─ Packages ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
 package       * version  date       lib source        
 AnnotationDbi   1.48.0   2019-10-29 [1] Bioconductor  
 askpass         1.1      2019-01-13 [1] CRAN (R 3.6.0)
 assertthat      0.2.1    2019-03-21 [1] CRAN (R 3.6.0)
 Biobase         2.46.0   2019-10-29 [1] Bioconductor  
 BiocFileCache   1.10.2   2019-11-08 [1] Bioconductor  
 BiocGenerics    0.32.0   2019-10-29 [1] Bioconductor  
 biomaRt       * 2.42.0   2019-10-29 [1] Bioconductor  
 bit             1.1-15.2 2020-02-10 [1] CRAN (R 3.6.0)
 bit64           0.9-7    2017-05-08 [1] CRAN (R 3.6.0)
 blob            1.2.1    2020-01-20 [1] CRAN (R 3.6.0)
 cli             2.0.2    2020-02-28 [1] CRAN (R 3.6.0)
 crayon          1.3.4    2017-09-16 [1] CRAN (R 3.6.0)
 curl            4.3      2019-12-02 [1] CRAN (R 3.6.0)
 DBI             1.1.0    2019-12-15 [1] CRAN (R 3.6.0)
 dbplyr          1.4.2    2019-06-17 [1] CRAN (R 3.6.0)
 digest          0.6.25   2020-02-23 [1] CRAN (R 3.6.0)
 dplyr         * 0.8.5    2020-03-07 [1] CRAN (R 3.6.0)
 evaluate        0.14     2019-05-28 [1] CRAN (R 3.6.0)
 fansi           0.4.1    2020-01-08 [1] CRAN (R 3.6.0)
 glue            1.3.1    2019-03-12 [1] CRAN (R 3.6.0)
 hms             0.5.3    2020-01-08 [1] CRAN (R 3.6.0)
 htmltools       0.4.0    2019-10-04 [1] CRAN (R 3.6.0)
 httr            1.4.1    2019-08-05 [1] CRAN (R 3.6.0)
 IRanges         2.20.2   2020-01-13 [1] Bioconductor  
 janitor       * 1.2.1    2020-01-22 [1] CRAN (R 3.6.0)
 knitr           1.28     2020-02-06 [1] CRAN (R 3.6.0)
 lifecycle       0.2.0    2020-03-06 [1] CRAN (R 3.6.0)
 magrittr        1.5      2014-11-22 [1] CRAN (R 3.6.0)
 memoise         1.1.0    2017-04-21 [1] CRAN (R 3.6.0)
 openssl         1.4.1    2019-07-18 [1] CRAN (R 3.6.0)
 packrat         0.5.0    2018-11-14 [1] CRAN (R 3.6.0)
 pillar          1.4.3    2019-12-20 [1] CRAN (R 3.6.0)
 pkgconfig       2.0.3    2019-09-22 [1] CRAN (R 3.6.0)
 prettyunits     1.1.1    2020-01-24 [1] CRAN (R 3.6.0)
 progress        1.2.2    2019-05-16 [1] CRAN (R 3.6.0)
 purrr         * 0.3.3    2019-10-18 [1] CRAN (R 3.6.0)
 R6              2.4.1    2019-11-12 [1] CRAN (R 3.6.0)
 rappdirs        0.3.1    2016-03-28 [1] CRAN (R 3.6.0)
 Rcpp            1.0.3    2019-11-08 [1] CRAN (R 3.6.0)
 readr         * 1.3.1    2018-12-21 [1] CRAN (R 3.6.0)
 rlang           0.4.5    2020-03-01 [1] CRAN (R 3.6.0)
 rmarkdown       2.1      2020-01-20 [1] CRAN (R 3.6.0)
 RSQLite         2.2.0    2020-01-07 [1] CRAN (R 3.6.0)
 S4Vectors       0.24.3   2020-01-18 [1] Bioconductor  
 sessioninfo     1.1.1    2018-11-05 [1] CRAN (R 3.6.0)
 stringi         1.4.6    2020-02-17 [1] CRAN (R 3.6.0)
 stringr       * 1.4.0    2019-02-10 [1] CRAN (R 3.6.0)
 tibble        * 2.1.3    2019-06-06 [1] CRAN (R 3.6.0)
 tidyr         * 1.0.2    2020-01-24 [1] CRAN (R 3.6.0)
 tidyselect      1.0.0    2020-01-27 [1] CRAN (R 3.6.0)
 vctrs           0.2.3    2020-02-20 [1] CRAN (R 3.6.0)
 withr           2.1.2    2018-03-15 [1] CRAN (R 3.6.0)
 xfun            0.12     2020-01-13 [1] CRAN (R 3.6.0)
 XML             3.99-0.3 2020-01-20 [1] CRAN (R 3.6.0)
```

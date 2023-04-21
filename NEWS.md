# dorothea 1.7.3 (2023-04-21)
* Added CollecTRI and deprecated run_viper function for `decoupleR`.

# dorothea 1.7.1 (2022-04-07)
* Removed outdated vignettes, now instead we point to decoupleR and decoupler-py's 
most up-to-date vignettes.

# dorothea 1.5.2 (2021-10-13)
* Daniel Dimitrov is assigned as the new maintainer

# dorothea 1.4.2 (2021-10-08)
* Fixed lazy data warning
* Improved test coverage

# dorothea 1.4.1 (2021-05-25)
* Rebuild all regulons
* Fixed ambiguously mode of regulation in mouse regulons

# dorothea 1.3.2 (2021-03-09)
* Added pancancer regulons for application in cancer.

# dorothea 1.2.1 (2021-02-08)
* Fixed bug in Seurat's related unit tests due to Seurats package update to version 4.0. `s@assays$dorothea@misc` is now `list()`, before it was `NULL`.

# dorothea 1.1.3 (2020-10-26)
* Fixed bug in single-cell vignette

# dorothea 1.1.2 (2020-10-08)
* Changed TF census from [TFclass](https://doi.org/10.1093/nar/gkx987) to the more recent version from [Lambert et al.](10.1016/j.cell.2018.01.029). Information of mode of regulation for each TF (activator, supressor, dual) is still taken from [Garcia-Alonso et al.](http://www.genome.org/cgi/doi/10.1101/gr.240663.118).
* Updated  gene symbols to their latest alias with the limma package (version 3.44.3).
* Shifted viper package from `suggest` to `depends` in the DESCRIPTION file.
* Added a further argument specifially for `run_viper.Seurat()` to select a specific assay name to extract the normalized gene expression values from.

# dorothea 1.1.1 (2020-09-02)
* Export `df2regulon()` function
* Improved documentation (added gh page URL to DESCRIPTION)

# dorothea 1.0.1 (2020-08-13)
* Improved package documentation
* Updated link to 10x genomics data set in single-cell vignette
* Fixed tests related to Seurat and SCE class

# dorothea 1.0.0 (2020-04-27)
* Official release in Bioconductor 3.11

# dorothea 0.99.9 (2020-04-22)
* Integration of `Travis CI`
* Integrated unit test coverage with `covr::codecov`
* Build github page via `pkgdown`

# dorothea 0.99.1 (2020-04-06)
* If the input to the viper wrapper is a Bioconductor class, this Bioconductor
class is retured with added TF activities at appropiate slots
* Expanded the viper wrapper to the Bioconductor class `SingleCellExperiment`

# dorothea 0.99.0 (2020-04-03)
* Initial submission to Bioconductor


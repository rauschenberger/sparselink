[![AppVeyor build status](https://ci.appveyor.com/api/projects/status/github/rauschenberger/sparselink?svg=true)](https://ci.appveyor.com/project/rauschenberger/sparselink)
[![R-CMD-check](https://github.com/rauschenberger/sparselink/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/rauschenberger/sparselink/actions/workflows/R-CMD-check.yaml)
[![Codecov test coverage](https://codecov.io/gh/rauschenberger/sparselink/graph/badge.svg)](https://app.codecov.io/gh/rauschenberger/sparselink)

# Sparse regression for related problems

Estimates sparse regression models (i.e., with few non-zero coefficients) in high-dimensional multi-task learning and transfer learning settings

## Installation

Install the current release from
[CRAN](https://CRAN.R-project.org/package=sparselink):

``` r
install.packages("sparselink")
```

or the latest development version from [GitHub](https://github.com/lcsb-bds/sparselink) or [GitLab](https://gitlab.lcsb.uni.lu/bds/sparselink):

``` r
#install.packages("remotes")
remotes::install_github("rauschenberger/sparselink")
remotes::install_gitlab("bds/sparselink",host="gitlab.lcsb.uni.lu") # mirror
```

## Reference

Armin Rauschenberger 
[![AR](https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png)](https://orcid.org/0000-0001-6498-4801),
Petr V. Nazarov[![PVN](https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png)](https://orcid.org/0000-0003-3443-0298),
and Enrico Glaab
[![EG](https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png)](https://orcid.org/0000-0003-3977-7469) (2023).
"Estimating sparse regression models in multi-task learning and transfer learning through adaptive penalisation".
*Manuscript in preparation*.

## Reproducibility

The code for reproducing the simulations and applications shown in the manuscript is available in a vignette (...). After installing the package with `remotes::install_github("rauschenberger/sparselink",build_vignettes=TRUE)` and restarting R, the vignette can also be loaded with `vignette(topic="...",package="sparselink")`.

<!--
[![CRAN version](https://www.r-pkg.org/badges/version/sparselink)](https://CRAN.R-project.org/package=sparselink)
[![CRAN RStudio mirror downloads](https://cranlogs.r-pkg.org/badges/sparselink)](https://CRAN.R-project.org/package=sparselink)
[![Total CRAN downloads](https://cranlogs.r-pkg.org/badges/grand-total/sparselink)](https://CRAN.R-project.org/package=sparselink)
-->

## Disclaimer

The R package `sparselink` implements sparse regression for related problems (Rauschenberger et al., 2024).

Copyright &copy; 2024 Armin Rauschenberger; Luxembourg Institute of Health (LIH), Department of Medical Informatics (DMI), Bioinformatics and Artificial Intelligence (BioAI); University of Luxembourg, Luxembourg Centre for Systems Biomedicine (LCSB), Biomedical Data Science (BDS)

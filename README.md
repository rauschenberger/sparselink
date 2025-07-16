
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

or the latest development version from [GitHub](https://github.com/rauschenberger/sparselink):

``` r
#install.packages("remotes")
remotes::install_github("rauschenberger/sparselink")
```

This repository is mirrored on two institutional GitLab instances
(see [LIH](https://git.lih.lu/bioinformatics-and-ai/sparselink) and 
[LCSB](https://gitlab.com/uniluxembourg/lcsb/bds/sparselink)).

## Reference

Armin Rauschenberger 
[![AR](https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png)](https://orcid.org/0000-0001-6498-4801),
Petr V. Nazarov
[![PVN](https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png)](https://orcid.org/0000-0003-3443-0298),
and Enrico Glaab
[![EG](https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png)](https://orcid.org/0000-0003-3977-7469) (2025).
"Estimating sparse regression models in multi-task learning and transfer learning through adaptive penalisation".
*Bioinformatics*. [doi: 10.1093/bioinformatics/btaf406](https://doi.org/10.1093/bioinformatics/btaf406)

## Reproducibility

The code for reproducing the simulations and applications shown in the manuscript is available in a vignette ([analysis](https://rauschenberger.github.io/sparselink/articles/analysis.html)). After installing the package with `remotes::install_github("rauschenberger/sparselink",build_vignettes=TRUE)` and restarting R, the vignette can also be loaded with `vignette(topic="analysis",package="sparselink")`.

<!--
[![CRAN version](https://www.r-pkg.org/badges/version/sparselink)](https://CRAN.R-project.org/package=sparselink)
[![CRAN RStudio mirror downloads](https://cranlogs.r-pkg.org/badges/sparselink)](https://CRAN.R-project.org/package=sparselink)
[![Total CRAN downloads](https://cranlogs.r-pkg.org/badges/grand-total/sparselink)](https://CRAN.R-project.org/package=sparselink)
-->

## Disclaimer

The R package `sparselink` implements sparse regression for related problems (Rauschenberger et al., 2025).

Copyright &copy; 2025 Armin Rauschenberger; Luxembourg Institute of Health (LIH), Department of Medical Informatics (DMI), Bioinformatics and Artificial Intelligence (BioAI); University of Luxembourg, Luxembourg Centre for Systems Biomedicine (LCSB), Biomedical Data Science (BDS)

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

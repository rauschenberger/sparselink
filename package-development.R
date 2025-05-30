
# --- code for package development ---

rm(list=ls(all.names=TRUE))
#install.packages(c("roxygen2","pkgdown","rcmdcheck","usethis","remotes","testthat","devtools"))
setwd("C:/Users/arauschenberger/Desktop/sparselink/package")
# usethis::use_mit_license()
roxygen2::roxygenise()
rcmdcheck::rcmdcheck()
pkgdown::check_pkgdown()
#devtools::build()
#devtools::submit_cran()
#pkgdown::build_site()
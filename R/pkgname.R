
#' @name sparselink-package
#' @keywords documentation
#' @docType package
#' 
#' @aliases sparselink-package
#' 
#' @title
#' Sparse regression for related problems
#' 
#' @description
#' The R package `sparselink` implements sparse
#' regression for related problems
#' (multi-task learning and transfer learning).
#' 
#' @details
#' Use function [sparselink()] for model fitting.
#' Type `library(sparselink)` and then `?sparselink` or
#' `help("sparselink")` to open its help file.
#' 
#' See the vignette for further examples.
#' Type `vignette("sparselink")` or `browseVignettes("sparselink")`
#' to open the vignette.
#' 
#' @seealso
#' First use \code{\link{sparselink}} to fit the models,
#' and then \code{\link[=coef.sparselink]{coef}} to extract coefficients
#' or \code{\link[=predict.sparselink]{predict}} to make predictions.
#' 
#' @references
#' \href{https://orcid.org/0000-0001-6498-4801}{Armin Rauschenberger},
#' \href{https://orcid.org/0000-0003-3443-0298}{Petr N. Nazarov}, and
#' \href{https://orcid.org/0000-0003-3977-7469}{Enrico Glaab}
#' (2025).
#' "Estimating sparse regression models in multi-task learning and transfer learning through adaptive penalisation".
#' \emph{Bioinformatics}. \doi{10.1093/bioinformatics/btaf406}
#' 
#' @examples
#' ?sparselink
#' ?coef.sparselink
#' ?predict.sparselink
#' 
"_PACKAGE"

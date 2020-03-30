#' Human DoRothEA
#'
#' A table reporting signed human TF-target interactions. This database covers
#'   in total 1395 TFs targeting 20,244 genes with 486,676 unique interactions.
#'   In addition, each TF is accompanied with an emperical confidence level that
#'   was derived from the number of supporting evidences for this
#'   TF/interaction. The range is from A (high quality) to E (low quality).
#'
#'
#' @format A table of human TF-target interactions:
#' \describe{
#'     \item{tf}{TF identifier as HGNC symbols}
#'     \item{confidence}{Confidence score classifying regulons based on their quality}
#'     \item{target}{target identifier as HGNC symbols}
#'     \item{mor}{mode of regulation indicating the effect of a TF on the target}
#' }
#'
#' @keywords datasets
#' @name dorothea_hs

#' @source \url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6673718/}
NULL

#' Mouse DoRothEA
#'
#' A table reporting signed human TF-target interactions. This database covers
#'   in total 1179 TFs targeting 17,410 genes with 410,455 unique interactions.
#'   In addition, each TF is accompanied with an emperical confidence level that
#'   was derived from the number of supporting evidences for this
#'   TF/interaction. The range is from A (high quality) to E (low quality).
#'
#'
#' @format A table of mouse TF-target interactions:
#' \describe{
#'     \item{tf}{TF identifier as MGI symbols}
#'     \item{confidence}{Confidence score classifying regulons based on their quality}
#'     \item{target}{target identifier as MGI symbols}
#'     \item{mor}{mode of regulation indicating the effect of a TF on the target}
#' }
#'
#' @keywords datasets
#' @name dorothea_mm
#' @example data(dorothea_mm, package = "dorothea")
#' @source \url{https://www.ncbi.nlm.nih.gov/pubmed/31525460}
NULL

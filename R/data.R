#' Human DoRothEA
#'
#' A table reporting signed human TF-target interactions. This database covers
#'   in total 1395 TFs targeting 20,244 genes with 486,676 unique interactions.
#'   In addition, each TF is accompanied with an empirical confidence level that
#'   was derived from the number of supporting evidences for this
#'   TF/interaction. The range is from A (high quality) to E (low quality).
#'
#'
#' @format A table of human TF-target interactions:
#' \describe{
#'     \item{tf}{TF identifier as HGNC symbols}
#'     \item{confidence}{Summary confidence score classifying regulons based on
#'         their quality}
#'     \item{target}{target identifier as HGNC symbols}
#'     \item{mor}{mode of regulation indicating the effect of a TF on the
#'         target}
#' }
#'
#' @keywords datasets
#' @name dorothea_hs
#' @source \url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6673718/}
NULL

#' Human DoRothEA (pancancer)
#'
#' A table reporting signed human TF-target interactions for cancer application. 
#'   The difference to "dorothea_hs" is that TGCA gene expression data was used 
#'   instead of GTEx to infer the networks with ARACNE. This database covers
#'   in total 1344 TFs targeting 20,582 genes with 213,230 unique interactions.
#'   In addition, each TF is accompanied with an empirical confidence level that
#'   was derived from the number of supporting evidences for this
#'   TF/interaction. The range is from A (high quality) to E (low quality).
#'
#'
#' @format A table of human TF-target interactions:
#' \describe{
#'     \item{tf}{TF identifier as HGNC symbols}
#'     \item{confidence}{Summary confidence score classifying regulons based on
#'         their quality}
#'     \item{target}{target identifier as HGNC symbols}
#'     \item{mor}{mode of regulation indicating the effect of a TF on the
#'         target}
#' }
#'
#' @keywords datasets
#' @name dorothea_hs_pancancer
#' @source \url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6673718/}
NULL

#' Mouse DoRothEA
#'
#' A table reporting signed human TF-target interactions. This database covers
#'   in total 1179 TFs targeting 17,410 genes with 410,455 unique interactions.
#'   In addition, each TF is accompanied with an empirical confidence level that
#'   was derived from the number of supporting evidences for this
#'   TF/interaction. The range is from A (high quality) to E (low quality).
#'
#'
#' @format A table of mouse TF-target interactions:
#' \describe{
#'     \item{tf}{TF identifier as MGI symbols}
#'     \item{confidence}{summary confidence score classifying regulons based on
#'         their quality}
#'     \item{target}{target identifier as MGI symbols}
#'     \item{mor}{mode of regulation indicating the effect of a TF on the
#'         target}
#' }
#'
#' @keywords datasets
#' @name dorothea_mm
#' @source \url{https://www.ncbi.nlm.nih.gov/pubmed/31525460}
NULL

#' Mouse DoRothEA (pancancer)
#'
#' A table reporting signed human TF-target interactions. The difference to 
#'   "dorothea_mm" is that TGCA gene expression data was used instead of GTEx to 
#'   infer the networks with ARACNE.This database covers in total 1096 TFs 
#'   targeting 17,695 genes with 187,955 unique interactions. In addition, each 
#'   TF is accompanied with an empirical confidence level that was derived from 
#'   the number of supporting evidences for this TF/interaction. The range is 
#'   from A (high quality) to E (low quality).
#'
#'
#' @format A table of mouse TF-target interactions:
#' \describe{
#'     \item{tf}{TF identifier as MGI symbols}
#'     \item{confidence}{summary confidence score classifying regulons based on
#'         their quality}
#'     \item{target}{target identifier as MGI symbols}
#'     \item{mor}{mode of regulation indicating the effect of a TF on the
#'         target}
#' }
#'
#' @keywords datasets
#' @name dorothea_mm_pancancer
#' @source \url{https://www.ncbi.nlm.nih.gov/pubmed/31525460}
NULL

#' Entire database with associated meta data
#'
#' This table lists all human TF-target interactions that were derived from the
#'   four lines of evidences. Each interaction is assigned a confidence score
#'   based on the number of supporting evidences. The table provides also all
#'   required information to trace back the origin of the interaction.
#'
#'
#' @keywords datasets
#' @name entire_database
#' @source \url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6673718/}
NULL

#' Construction of viper regulons
#'
#' This function converts DoRothEA's regulons that are stored in a table to the
#' format required by the \code{\link[=viper]{viper::viper()}} function.
#'
#' @param df A table from DoRothEA package.
#' @param confidence_levels Combination of DoRothEA's confidence levels.
#'   Corresponding regulons/interactions will be returned. This combination
#'   can contain any combination from the letters A-E and must be encoded as a
#'   single character string.
#' @param organism A character string indicating whether human DoRothEA (hs) or
#'   mouse DoRothEA (mm) is desired.
#'
#' @return Regulons in the \code{\link[=viper]{viper::viper()}} format.
#' @import dplyr
#' @export
df2regulon = function(df, confidence_levels = c("ABCDE"), organism = "hs") {

  regulon_list = df %>%
    filter(confidence %in% parse_confidence_levels(confidence_levels)) %>%
    split(.$tf)

  viper_regulons = lapply(regulon_list, function(regulon) {
    tfmode = setNames(regulon$mor, regulon$target)
    list(tfmode = tfmode, likelihood = rep(1, length(tfmode)))
  })

  return(viper_regulons)
}

#' Construction of gene sets
#'
#' This function coverts DoRothEA's regulons that are stored in a table to the
#' format required for GSEA methods (e.g. \code{\link[=fgsea]{fgsea::fgsea()}}).
#'
#' @param df A table from DoRothEA package.
#' @param confidence_levels Combination of DoRothEA's confidence levels.
#'   Corresponding regulons/interactions will be returned. This combination
#'   can contain any combination from the letters A-E and must be encoded as a
#'   single character string.
#' @param organism A character string indicating whether human DoRothEA (hs) or
#'   mouse DoRothEA (mm) is desired.
#'
#' @return Regulons in GSEA format (see \code{\link[=fgsea]{fgsea::fgsea()}}).
#'
#' @import dplyr
#' @export
df2geneset = function(df, confidence_levels = c("ABCDE"),
                      organism = "hs") {

  regulon_list = df %>%
    filter(confidence %in% parse_confidence_levels(confidence_levels)) %>%
    split(.$tf)

  gsea_regulons = lapply(regulon_list, function(regulon) {
    regulon$target
  })
}


#' Parsing confidende levels
#'
#' @param confidence_levels Combination of DoRothEA's confidence levels. This
#'   combination can contain any combination from the letters A-E and must be
#'   encoded as a single character string.
#'
#' @return Vector of confidence levels.
parse_confidence_levels = function(confidence_levels) {

  v = unlist(strsplit(x = confidence_levels, split = ""))

  if (!all(grepl(pattern = "[A|B|C|D|E]", v))) {
    stop("Confidence levels must contain only the letters A,B,C,D and E")
  }

  return(v)
}


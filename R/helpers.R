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

#' VIPER wrapper for single cell data
#'
#' This function is a convenient wrapper to run the
#' \code{\link[=viper]{viper::viper()}} function to apply it on single cell
#' data combined with Dorothea regulons
#'
#' @import viper
#' @export
#' @param InputObject An expression matrix with genes (HGNC/MGI symbols) in rows
#' and cells in columns. It also accepts \code{Seurat} object from which the
#' normalised expresion matrix of the selected assay will be extracted
#' @param regulon Object of class regulon. Check viper documentation for
#' further information.
#' @param assay_name it only applies if the input is a Seurat object. It selects
#' the name of the assay on which Viper will be run. Default to:
#' RNA, i.e. normalized expression values.
#' @param return_assay it only applies if the input is a Seurat object. A
#' logical value indicating whether to return Dorothea results as a new
#' assay called Dorothea in the Seurat object used as input.
#' Default to FALSE.
#' @param options A list of named options to pass to
#' \code{\link[=viper]{viper::viper()}} such as \code{minsize} or
#' \code{method}. These options should not include \code{eset} or
#' \code{regulon}.
#'
#' @return A matrix containing the activity of the different TFs provided in
#' the regulon object.

scViper = function(InputObject, regulon, assay_name = "RNA",
    return_assay = FALSE, options=list()) {

    requireNamespace("Seurat")

    if (class(InputObject) == "Seurat"){
        expr <- InputObject[[assay_name]]@data
    } else {
        expr <- InputObject
    }

    emat <- as.matrix(expr)
    viper_res  <- do.call(viper,c(list(eset = emat,regulon = regulon),options))

    if (return_assay){
        InputObject[['dorothea']] <-
            Seurat::CreateAssayObject(data = viper_res)
        Seurat::Key(object = InputObject[['dorothea']]) <- 'dorothea_'
        return(InputObject)
    } else {
        return(viper_res)
    }
}
>>>>>>> 04881b407527bb245ef1833511f155afecc0f718

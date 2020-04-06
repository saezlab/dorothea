#' VIPER wrapper
#'
#' This function is a convenient wrapper for the
#' \code{viper::\link[viper]{viper}} function using DoRothEA regulons.
#'
#' @param input An object containing a gene expression matrix with genes
#'   (HGNC/MGI symbols) in rows and samples in columns. The object can be a
#'   simple matrix/data frame or complexer objects such as
#'   \code{Biobase::\link[Biobase:ExpressionSet]{ExpressionSet}},
#'   \code{Seurat::\link[Seurat:Seurat]{Seurat}} or
#'   \code{SingleCellExperiment::\link[SingleCellExperiment:SingleCellExperiment]{SingleCellExperiment}} objects.
#' @param regulons \code{\link[=dorothea_hs]{DoRothEA}} regulons in table
#'   format.
#' @param options A list of named options to pass to
#'   \code{viper::\link[viper]{viper}} such as \code{minsize} or
#'   \code{method}. These options should not include, \code{eset} or
#'   \code{regulon}.
#' @param tidy Logical, whether computed tf activities scores should be returned
#'  in a tidy format.
#'
#' @return A matrix of normalized enrichment scores for each TF across all
#'  samples. Of note, if you provide Bioconductor objects as input the function
#'  will return this object with added tf activities at appropiate slots. e.g.
#'  Seurat object with a new assay called \code{dorothea}. For all
#'  other inputs the function will return a matrix. If \code{tidy} is
#'  \code{TRUE} the normalized enrichment scores are retured in a tidy format
#'  (not supported for Bioconductor objects).
#'
#' @export
#' @import dplyr
#' @import bcellViper
#'
#' @examples
#' # use example gene expression matrix from bcellViper package
#' library(bcellViper)
#' data(bcellViper, package = "bcellViper")
#' # acessing (human) dorothea regulons
#' # for mouse regulons: data(dorothea_mm, package = "dorothea")
#' data(dorothea_hs, package = "dorothea")
#' # run viper
#'tf_activities <- run_viper(dset, dorothea_hs,
#'                           options =  list(method = "scale", minsize = 4,
#'                           eset.filter = FALSE, cores = 1,
#'                           verbose = FALSE))
run_viper <- function(input, regulons, options = list(), tidy = FALSE) {
  UseMethod("run_viper")
}

#' @export
run_viper.ExpressionSet <- function(input, regulons, options = list(),
                                    tidy=FALSE) {

  if (tidy) {
    tidy <- FALSE
    warning("The argument 'tidy' cannot be TRUE for 'ExpressionSet' objects. ",
            "'tidy' is set to FALSE")
  }

  tf_activities <- run_viper(Biobase::exprs(input), regulons = regulons,
                             options = options, tidy = tidy)

  eset <- Biobase::ExpressionSet(assayData = tf_activities,
                                 phenoData = Biobase::phenoData(input))

  return(eset)
}

#' @export
run_viper.data.frame <- function(input, regulons, options = list(),
                                 tidy=FALSE) {
  run_viper(as.matrix(input), regulons = regulons, options = options,
            tidy = tidy)
}

#' @export
run_viper.SingleCellExperiment <- function(input, regulons, options = list(),
                                           tidy = FALSE) {
  if (tidy) {
    tidy <- FALSE
    warning("The argument 'tidy' cannot be TRUE for 'SingleCellExperiment' ",
            "objects. ","'tidy' is set to FALSE")
  }

  mat <- as.matrix(SingleCellExperiment::logcounts(input))

  tf_activities <- run_viper(mat, regulons = regulons, options = options,
                             tidy = FALSE)

  # include TF activities into SingleCellExperiment object
  dorothea_se <- SummarizedExperiment::SummarizedExperiment(tf_activities)
  SummarizedExperiment::assayNames(dorothea_se) <- "tf_activities"
  SingleCellExperiment::altExp(input, "dorothea") <- dorothea_se

  return(input)
}

#' @export
run_viper.Seurat <- function(input, regulons, options = list(), tidy = FALSE) {
  if (tidy) {
    tidy <- FALSE
    warning("The argument 'tidy' cannot be TRUE for 'Seurat' objects. ",
            "'tidy' is set to FALSE")
  }
  mat <- as.matrix(Seurat::GetAssayData(input, assay = "RNA", slot = "data"))

  tf_activities <- run_viper(mat, regulons = regulons, options = options,
                             tidy = FALSE)

  # include TF activities in Seurat object
  dorothea_assay <- Seurat::CreateAssayObject(data = tf_activities)
  Seurat::Key(dorothea_assay) <- "dorothea_"
  input[["dorothea"]] <- dorothea_assay

  return(input)
}

#' @export
run_viper.matrix <- function(input, regulons, options = list(), tidy=FALSE) {
  viper_res <- do.call(viper::viper,
                      c(list(eset = input,
                             regulon = .df2regulon(regulons)),
                        options))

  if (tidy) {
    metadata <- regulons %>%
      dplyr::select(-c("target", "mor")) %>%
      dplyr::distinct()

    tidy_viper_res <- viper_res %>%
      data.frame(check.names = FALSE, stringsAsFactors = FALSE) %>%
      tibble::rownames_to_column("tf") %>%
      tidyr::gather(sample, "activity", -c("tf")) %>%
      dplyr::inner_join(metadata, by="tf")

    return(tidy_viper_res)
  } else {
    return(viper_res)
  }
}

#' @export
run_viper.default <- function(input, regulons, options=list(), tidy = FALSE) {
  stop("Do not know how to access the data matrix from class ", class(input))
}


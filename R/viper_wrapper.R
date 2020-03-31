#' VIPER wrapper
#'
#' This function is a convenient wrapper for the
#' \code{\link[=viper]{viper()}} function using DoRothEA regulons.
#'
#' @param input An object containing a gene expression matrix with genes
#'   (HGNC/MGI symbols) in rows and samples in columns. The object can be a simple
#'   matrix/data frame or complexer objects such as ExpressionSet or Seurat
#'   objects.
#' @param regulons DoRothEA regulons in table format.
#' @param options A list of named options to pass to \code{viper} such as
#'   \code{minsize} or \code{method}. These options should not include,
#'   \code{eset} or \code{regulon}.
#' @param tidy Logical, whether computed viper scores should be returned in a
#' tidy format.
#'
#' @return A matrix of normalized enrichment scores for each TF across all
#'  samples. Of note, if you provide a Seurat object as input the function will
#'  return also a Seurat object with a new assay called \code{dorothea} For all
#'  other inputs the function will return a matrix. If \code{tidy} is
#'  \code{TRUE} the normalized enrichment scores are retured in a tidy format
#'  (not supported for Seurat objects).
#'
#' @export
#' @import dplyr
run_viper <- function(input, regulons, options = list(), tidy = F) {
  UseMethod("run_viper")
}

#' @export
run_viper.ExpressionSet <- function(input, regulons, options = list(), tidy=F) {
  run_viper(input@assayData$exprs, regulons = regulons, options = options, tidy = tidy)
}

#' @export
run_viper.data.frame <- function(input, regulons, options = list(), tidy=F) {
  run_viper(as.matrix(input), regulons = regulons, options = options,
            tidy = tidy)
}

#' @export
run_viper.Seurat <- function(input, regulons, options = list(), tidy = F) {
  if (tidy) {
    tidy <- F
    warning("The argument 'tidy' cannot be true for Seurat objects. tidy is set to FALSE")
  }
  mat <- as.matrix(Seurat::GetAssayData(input, assay = "RNA", slot = "data"))

  tf_activities <- run_viper(mat, regulons = regulons, options = options,
                             tidy = F)

  # include TF activities in Seurat object
  dorothea_assay <- Seurat::CreateAssayObject(data = tf_activities)
  Seurat::Key(dorothea_assay) <- "dorothea_"
  input[["dorothea"]] <- dorothea_assay

  return(input)
}

#' @export
run_viper.matrix <- function(input, regulons, options = list(), tidy=F) {
  viper_res <- do.call(viper::viper,
                      c(list(eset = input,
                             regulon = df2regulon(regulons)),
                        options))

  if (tidy) {
    metadata <- regulons %>%
      select(-c("target", "mor")) %>%
      distinct()

    tidy_viper_res <- viper_res %>%
      data.frame(check.names = F, stringsAsFactors = F) %>%
      tibble::rownames_to_column("tf") %>%
      tidyr::gather(sample, "activity", -c("tf")) %>%
      inner_join(metadata, by="tf")

    return(tidy_viper_res)
  } else {
    return(viper_res)
  }
}

#' @export
run_viper.default <- function(input, regulons, options=list(), tidy = F) {
  stop("Do not know how to access the data matrix from class ", class(input))
}


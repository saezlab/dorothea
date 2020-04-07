test_that("test run_viper with matrix as input", {
  m <- readRDS(
    system.file("testdata", "toy_matrix.rds", package = "dorothea")
  )
  r <- dplyr::filter(dorothea_hs, confidence %in% c("A", "B"))
  res <- run_viper(m, r, options =  list(method = "scale", minsize = 4,
                                         eset.filter = FALSE, verbose = FALSE))

  tidy_res <- run_viper(m, r, options =  list(method = "scale", minsize = 4,
                                              eset.filter = FALSE,
                                              verbose = FALSE),
                        tidy = TRUE)

  expected_res <- readRDS(
    system.file("testdata", "output_matrix.rds", package = "dorothea")
    )

  expected_res_tidy <- readRDS(
    system.file("testdata", "output_matrix_tidy.rds", package = "dorothea")
  )

  expect_equal(res, expected_res)
  expect_equal(tidy_res, expected_res_tidy)
})

test_that("test run_viper with eset as input", {
  m <- readRDS(
    system.file("testdata", "toy_eset.rds", package = "dorothea")
  )
  r <- dplyr::filter(dorothea_hs, confidence %in% c("A", "B"))

  res <- run_viper(m, r, options =  list(method = "scale", minsize = 4,
                                         eset.filter = FALSE, verbose = FALSE))
  suppressWarnings(
    tidy_res <- run_viper(m, r, options =  list(method = "scale", minsize = 4,
                                                eset.filter = FALSE,
                                                verbose = FALSE),
                          tidy = TRUE)
  )

  # we need to create the expected result on the fly because the created object 
  # is dependent on the used Rversion (see .@__classVersion__).
  expected_res <- viper::viper(m, .df2regulon(r), method = "scale", minsize=4,
                               eset.filter = F, verbose = F)

  expect_equal(res, expected_res)
  expect_equal(tidy_res, expected_res)

  # check raised warning when tidy is set to T
  expect_warning(
    run_viper(m, r, options = list(method = "scale", minsize = 4,
                                   eset.filter = FALSE,
                                   verbose = FALSE),
              tidy = TRUE),
    paste0("The argument 'tidy' cannot be TRUE for 'ExpressionSet' objects. ",
           "'tidy' is set to FALSE"))

})


test_that("test run_viper with seurat as input", {
  library(Seurat)
  m <- readRDS(
    system.file("testdata", "toy_seurat.rds", package = "dorothea")
  )
  r <- dplyr::filter(dorothea_hs, confidence %in% c("A", "B"))

  res <- run_viper(m, r, options =  list(method = "scale", minsize = 4,
                                         eset.filter = FALSE, verbose = FALSE))

  suppressWarnings(
    tidy_res <- run_viper(m, r, options =  list(method = "scale", minsize = 4,
                                                eset.filter = FALSE,
                                                verbose = FALSE),
                          tidy = TRUE)
  )

  expected_res <- readRDS(
    system.file("testdata", "output_seurat.rds", package = "dorothea")
  )

  expect_equal(res, expected_res)
  expect_equal(tidy_res, expected_res)

  # check raised warning when tidy is set to T
  expect_warning(
    run_viper(m, r, options = list(method = "scale", minsize = 4,
                                   eset.filter = FALSE,
                                   verbose = FALSE),
              tidy = TRUE),
    paste0("The argument 'tidy' cannot be TRUE for 'Seurat' objects. ",
           "'tidy' is set to FALSE"))

  # Check key of seurat assays
  expect_equal(unname(Seurat::Key(res)), c("rna_","dorothea_"))
})


test_that("test run_viper with sce as input", {
  m <- readRDS(
    system.file("testdata", "toy_sce.rds", package = "dorothea")
  )
  r <- dplyr::filter(dorothea_mm, confidence %in% c("A", "B"))

  res <- run_viper(m, r, options =  list(method = "scale", minsize = 4,
                                         eset.filter = FALSE, verbose = FALSE))

  suppressWarnings(
    tidy_res <- run_viper(m, r, options =  list(method = "scale", minsize = 4,
                                                eset.filter = FALSE,
                                                verbose = FALSE),
                          tidy = TRUE)
  )

  expected_res <- readRDS(
    system.file("testdata", "output_sce.rds", package = "dorothea")
  )

  expect_equal(res, expected_res)
  expect_equal(tidy_res, expected_res)

  # check raised warning when tidy is set to T
  expect_warning(
    run_viper(m, r, options = list(method = "scale", minsize = 4,
                                   eset.filter = FALSE,
                                   verbose = FALSE),
              tidy = TRUE),
    paste0("The argument 'tidy' cannot be TRUE for 'SingleCellExperiment' ",
           "objects. ", "'tidy' is set to FALSE"))

  # Check name of alternative experiment
  expect_equal(altExpNames(res), "dorothea")
})


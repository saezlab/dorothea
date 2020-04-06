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


  expected_res <- readRDS(
    system.file("testdata", "output_eset.rds", package = "dorothea")
  )

  expect_equal(res, expected_res)
  expect_equal(tidy_res, expected_res)

  # check raised warning when tidy is set to T
  expect_warning(
    run_viper(m, r, options = list(method = "scale", minsize = 4,
                                   eset.filter = FALSE,
                                   verbose = FALSE),
              tidy = TRUE),
    paste0("The argument 'tidy' cannot be TRUE for ExpressionSet objects. ",
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
    paste0("The argument 'tidy' cannot be TRUE for Seurat objects. ",
           "'tidy' is set to FALSE"))

  # Check key of seurat assays
  expect_equal(unname(Seurat::Key(res)), c("rna_","dorothea_"))

})


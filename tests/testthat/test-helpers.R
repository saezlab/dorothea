test_that("test df2regulon", {

  res = df2regulon(dorothea_hs)

  # test if a list is returned
  expect_equal(class(res), "list")

  # test number of tfs
  expect_equal(length(res), 1333)

  # test random position in set
  expect_equal(names(res$FOXF2$tfmode[10]), "ABCB8")
})

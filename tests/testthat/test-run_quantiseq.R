# Running on Racle dataset...

test_that("Racle runs", {

  expect_message(
    res_quantiseq_run <- quantiseqr::run_quantiseq(
      expression_data = dataset_racle$expr_mat,
      signature_matrix = "TIL10",
      is_arraydata = FALSE,
      is_tumordata = TRUE,
      scale_mRNA = TRUE
    )
  )

  res_quantiseq_run_df <- quantiseqr::run_quantiseq(
    expression_data = as.data.frame(dataset_racle$expr_mat),
    signature_matrix = "TIL10",
    is_arraydata = FALSE,
    is_tumordata = TRUE,
    scale_mRNA = TRUE
  )

  expect_identical(res_quantiseq_run, res_quantiseq_run_df)

  expect_is(res_quantiseq_run, "data.frame")
  expect_equal(nrow(res_quantiseq_run), 4)

  p <- quantiplot(res_quantiseq_run)
  expect_is(p, "gg")

  # forcing the TPM value to be rounded to integer (simulating counts...)
  racle_integered <- dataset_racle$expr_mat
  storage.mode(racle_integered) <- "integer"
  expect_warning(
    racle_test_integer <- quantiseqr::run_quantiseq(
      expression_data = racle_integered,
      signature_matrix = "TIL10",
      is_arraydata = FALSE,
      is_tumordata = TRUE,
      scale_mRNA = TRUE
    )
  )

  # providing wrong inputs
  expect_error({
    racle_strings <- dataset_racle$expr_mat
    storage.mode(racle_strings) <- "character"
    res_quantiseq_run <- quantiseqr::run_quantiseq(
      expression_data = racle_strings,
      signature_matrix = "TIL10",
      is_arraydata = FALSE,
      is_tumordata = TRUE,
      scale_mRNA = TRUE
    )
  })

  expect_error(
    quantiseqr::run_quantiseq(
      expression_data = dataset_racle$expr_mat,
      signature_matrix = "TIL2020",
      is_arraydata = FALSE,
      is_tumordata = TRUE,
      scale_mRNA = TRUE
    )
  )
})



# using a SummarizedExperiment object
test_that("Racle runs as SE", {

  expect_message(
    res_run_SE <- quantiseqr::run_quantiseq(
      expression_data = se_racle,
      signature_matrix = "TIL10",
      is_arraydata = FALSE,
      is_tumordata = TRUE,
      scale_mRNA = TRUE
    )
  )

  expect_equal(
    sum(grepl("quanTIseq_TIL10", colnames(colData(res_run_SE)))),
    11
  )


  expect_is(res_run_SE, "SummarizedExperiment")
  expect_equal(dim(res_run_SE), c(32467, 4))

  expect_message(
    p1 <- quantiplot(res_run_SE)
  )
  expect_is(p1, "gg")

  res_run_SE_mixedup <- res_run_SE
  colnames(colData(res_run_SE_mixedup))[5] <- "quanTIseq_TIL2020_Monocytes"
  expect_error(
    extract_ti_from_se(res_run_SE_mixedup)
  )
})

test_that("Racle runs as ExpressionSet", {

  expect_message(
    res_run_ESet <- quantiseqr::run_quantiseq(
      expression_data = es_racle,
      signature_matrix = "TIL10",
      is_arraydata = FALSE,
      is_tumordata = TRUE,
      scale_mRNA = TRUE
    )
  )

  expect_equal(dim(res_run_ESet), c(4, 12))

})



test_that("Conversions work...", {
  mat_racle <- eset_to_matrix(eset = es_racle, column = "gene_symbol")
  expect_is(mat_racle, "matrix")

  from_se_racle <- se_to_matrix(se_racle)
  expect_is(from_se_racle, "matrix")

  expect_error(se_to_matrix(se = es_racle))
  expect_error(extract_ti_from_se(se = es_racle))

  expect_error(se_to_matrix(se = se_racle, assay = "something_not_there"))

  se_racle_mod <- se_racle
  assays(se_racle_mod) <- List(
    counts = assay(se_racle_mod)
  )

  expect_warning(se_to_matrix(se = se_racle_mod, assay = "counts"))

  expect_warning(se_to_matrix(se = se_racle_fakeENS))


})




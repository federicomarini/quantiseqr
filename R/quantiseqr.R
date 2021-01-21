######    #' Perform an immune cell deconvolution on a dataset.
######    #'
######    #' @param gene_expression A gene expression matrix or a Biobase ExpressionSet.
######    #'   Either: A numeric matrix or data.frame with HGNC gene symbols as rownames and sample ######    identifiers as colnames.
######    #'   Or: A Biobase ExpressionSet with HGNC symbols in an fData column (see `column` ######    parameter)
######    #'   In both cases, data must be on non-log scale.
######    #' @param column Only relevant if `gene_expression` is an ExpressionSet. Defines in which ######    column
######    #'   of fData the HGNC symbol can be found.
######    #' @param method a string specifying the method.
######    #'   Supported methods are `xcell`, `...`
######    #' @param indications a character vector with one indication per
######    #'   sample for TIMER. Argument is ignored for all other methods.
######    #' @param tumor use a signature matrix/procedure optimized for tumor samples,
######    #'   if supported by the method. Currently affects EPIC and quanTIseq.
######    #' @param arrays Runs methods in a mode optimized for microarray data.
######    #'   Currently affects quanTIseq and CIBERSORT.
######    #' @param rmgenes a character vector of gene symbols. Exclude these genes from the analysis.
######    #'   Use this to exclude e.g. noisy genes.
######    #' @param scale_mrna logical. If FALSE, disable correction for mRNA content of different ######    cell types.
######    #'   This is supported by methods that compute an absolute score (EPIC and quanTIseq)
######    #' @param expected_cell_types Limit the analysis to the cell types given in this list. If ######    the cell
######    #'   types present in the sample are known *a priori*, setting this can improve results for
######    #'   xCell (see https://github.com/grst/immunedeconv/issues/1).
######    #' @param ... arguments passed to the respective method
######    #' @return `data.frame` with `cell_type` as first column and a column with the
######    #'     calculated cell fractions for each sample.
######    #'
######    #' @examples
######    #' # Not run: deconvolute(gene_expression_matrix, "xcell")
######    #' @name deconvolute
######    #' @export deconvolute
######    deconvolute <- function(gene_expression, method = deconvolution_methods,
######                            indications = NULL, tumor = TRUE,
######                            arrays = FALSE, column = "gene_symbol",
######                            rmgenes = NULL, scale_mrna = TRUE,
######                            expected_cell_types = NULL,
######                            ...) {
######      message(paste0("\n", ">>> Running ", method))
######
######      # convert expression set to matrix, if required.
######      if (is(gene_expression, "ExpressionSet")) {
######        gene_expression <- gene_expression %>% eset_to_matrix(column)
######      }
######
######      if (!is.null(rmgenes)) {
######        gene_expression <- gene_expression[!rownames(gene_expression) %in% rmgenes, ]
######      }
######
######      # run selected method
######      res <- switch(method,
######                    xcell = deconvolute_xcell(gene_expression, arrays = arrays, ######    expected_cell_types = expected_cell_types, ...),
######                    mcp_counter = deconvolute_mcp_counter(gene_expression, ...),
######                    epic = deconvolute_epic(gene_expression, tumor = tumor, scale_mrna = ######    scale_mrna, ...),
######                    quantiseq = deconvolute_quantiseq(gene_expression,
######                                                      tumor = tumor, arrays = arrays, scale_mrna ######    = scale_mrna, ...
######                    ),
######                    cibersort = deconvolute_cibersort(gene_expression,
######                                                      absolute = FALSE,
######                                                      arrays = arrays, ...
######                    ),
######                    cibersort_abs = deconvolute_cibersort(gene_expression,
######                                                          absolute = TRUE,
######                                                          arrays = arrays, ...
######                    ),
######                    timer = deconvolute_timer(gene_expression, indications = indications, ...)
######      )
######
######      # convert to tibble and annotate unified cell_type names
######      res <- res %>%
######        as_tibble(rownames = "method_cell_type") %>%
######        annotate_cell_type(method = method)
######
######      return(res)
######    }






#' Deconvolute using quanTIseq (only needed for immunedeconv, not for quantiseqr -> can be removed)
#'
#' @param gene_expression_matrix a m x n matrix with m genes and n samples. Mandatory.
#' @param tumor Set to TRUE if dealing with a tumor samples. if TRUE, signature genes with
#'   high expression in tumor samples are removed. Default: FALSE
#' @param arrays Set to TRUE if working with Microarray data instead of RNA-seq. Default: FALSE
#' @param scale_mrna Set to FALSE to disable correction for cell type-specific differences
#'  in mRNA content. Default: TRUE
#' @param ... passed through to original quanTIseq method. A native argument takes precedence
#'   over an immunedeconv argument (e.g. `mRNAscale` takes precedence over `scale_mrna`)
#'   See `deconvolute_quantiseq.default()`.
#'
#' @export
deconvolute_quantiseq <- function(gene_expression_matrix,
                                  tumor,
                                  arrays,
                                  scale_mrna,
                                  ...) {

  # TODO: we might not even need this way of calling the function
  arguments <- rlang::dots_list(gene_expression_matrix, tumor = tumor, arrays = arrays, mRNAscale = scale_mrna, ..., .homonyms = "last")
  call <- rlang::call2(run_quantiseq, !!!arguments)
  res <- eval(call)

  sample_names <- res$Sample
  res_mat <- res %>%
    as_tibble() %>%
    select(-Sample) %>%
    as.matrix()
  rownames(res_mat) <- sample_names

  t(res_mat)
}










#' Use quanTIseq to deconvolute a gene expression matrix.
#'
#' Source code from https://github.com/FFinotello/quanTIseq - TODO: will need to update link?
#'
#' F. Finotello, C. Mayer, C. Plattner, G. Laschober, D. Rieder,
#' H. Hackl, A. Krogsdam, Z. Loncova, W. Posch, D. Wilflingseder, S. Sopper, M. Jsselsteijn,
#' T. P. Brouwer, D. Johnsons, Y. Xu, Y. Wang, M. E. Sanders, M. V. Estrada, P. Ericsson-Gonzalez,
#' P. Charoentong, J. Balko, N. F. d. C. C. de Miranda, Z. Trajanoski.
#' "Molecular and pharmacological modulators of the tumor immune contexture revealed by deconvolution of RNA-seq data".
#' Genome Medicine 2019;11(1):34. doi: 10.1186/s13073-019-0638-6.
#'
#' @param mix.mat table with the gene TPM or microarray expression values for all samples to be deconvoluted
#'     (Gene symbols on the first column and sample IDs on the first row). Expression data should be on non-log scale
#' @param arrays Set to TRUE if the expression data are from microarrays instead of RNA-seq. If TRUE, the "rmgenes" parameter is set to "none". Default: FALSE
#' @param signame name of the signature matrix. Currently only `TIL10` is available.
#' @param tumor 	Set to TRUE if the expression data is from tumor samples. Default: FALSE
#' @param mRNAscale Set to FALSE to disable the correction of cell-type-specific mRNA content bias. Default: TRUE
#' @param method deconvolution method to be used: "lsei" for constrained least squares regression, "hampel", "huber", or "bisquare" for robust regression with Huber, Hampel, or Tukey bisquare estimators. Default: "lsei".
#' @param column Character, specifies which column contains the information of the gene symbol identifiers
#' @param rmgenes Specifies which genes have to be removed from the deconvolution analysis.
#'  Can be a vector of gene symbols or a string among "none" (no genes are removed) and "default" (a list of genes with noisy expression RNA-seq data is removed as explained in quanTIseq paper). Default: "default" for RNA-seq data, "none" for microarrrays.
#'
#' @export
run_quantiseq <- function(mix.mat,
                          arrays = FALSE,
                          signame = "TIL10",
                          tumor = FALSE,
                          mRNAscale = TRUE,
                          method = "lsei",
                          column = "gene_symbol",
                          rmgenes = "unassigned") {




  # convert expression set to matrix, if required.
  if (is(mix.mat, "ExpressionSet")) {
    mix.mat <- mix.mat %>% eset_to_matrix(column)
  }






  if (!is.null(rmgenes)) {
    mix.mat <- mix.mat[!rownames(mix.mat) %in% rmgenes, ]
  }






  message("\nRunning quanTIseq deconvolution module\n")

  # List of genes to be discarded
  if (rmgenes == "unassigned" && arrays == TRUE) { # For Microarrays
    rmgenes <- "none"
  } else if (rmgenes == "unassigned" && arrays == FALSE) { # For RNA-seq
    rmgenes <- "default"
  }

  # Files
  listsig <- c("TIL10")
  if (signame %in% listsig) {
    sig.mat.file <- system.file("extdata", paste0(signame, "_signature.txt"),
      package = "quantiseqr", mustWork = TRUE
    )

    mRNA.file <- system.file("extdata", paste0(signame, "_mRNA_scaling.txt"),
      package = "quantiseqr", mustWork = TRUE
    )

    fileab <- system.file("extdata", paste0(signame, "_TCGA_aberrant_immune_genes.txt"),
      package = "quantiseqr", mustWork = TRUE
    )

    if (rmgenes == "default") {
      filerm <- system.file("extdata", paste0(signame, "_rmgenes.txt"),
        package = "quantiseqr", mustWork = TRUE
      )
    } else if (rmgenes == "path") {
      filerm <- system.file("extdata", paste0(signame, "rmgenes.txt"),
        package = "quantiseqr", mustWork = TRUE
      )
    }
  } else {
    sig.mat.file <- paste0(signame, "_signature.txt")
    mRNA.file <- paste0(signame, "_mRNA_scaling.txt")
  }

  if (is.numeric(mix.mat[[1, 1]]) != TRUE) {
    stop("Wrong input format for the mixture matrix! Please follow the instructions of the documentation.")
  }

  # Load signature
  sig.mat <- read.table(sig.mat.file, header = TRUE, sep = "\t", row.names = 1)

  # Load normalization factors (set all to 1 if mRNAscale==FALSE)
  if (mRNAscale) {
    mRNA <- read.table(mRNA.file,
      sep = "\t",
      header = FALSE,
      stringsAsFactors = FALSE
    )
    colnames(mRNA) <- c("celltype", "scaling")
    mRNA <- as.vector(as.matrix(mRNA$scaling[match(colnames(sig.mat), mRNA$celltype)]))
  } else {
    mRNA <- rep(1, ncol(sig.mat))
  }

  # Preprocess mixture matrix
  message(paste0(
    "Gene expression normalization and re-annotation (arrays: ",
    arrays, ")\n"
  ))
  mix.mat <- fixMixture(mix.mat, arrays = arrays)

  # Remove noisy genes
  if (rmgenes != "none") {
    if (signame %in% listsig) {
      lrmgenes <- as.vector(read.table(filerm, header = FALSE, sep = "\t")[, 1])
      n1 <- nrow(sig.mat)
      sig.mat <- sig.mat[!rownames(sig.mat) %in% lrmgenes, , drop = FALSE]
      n2 <- nrow(sig.mat)
      message(paste0("Removing ", n1 - n2, " noisy genes\n"))
    }
  }

  # Fix tumor data
  if (tumor) {
    if (signame %in% listsig) {
      abgenes <- as.vector(read.table(fileab, header = FALSE, sep = "\t")[, 1])
      n1 <- nrow(sig.mat)
      sig.mat <- sig.mat[!rownames(sig.mat) %in% abgenes, , drop = FALSE]
      n2 <- nrow(sig.mat)
      message(paste0("Removing ", n1 - n2, " genes with high expression in tumors\n"))
    }
  }

  # Signature genes present in the mixture
  ns <- nrow(sig.mat)
  us <- length(intersect(rownames(sig.mat), rownames(mix.mat)))
  perc <- round(us * 100 / ns, digits = 2)
  message(paste0(
    "Signature genes found in data set: ",
    us, "/", ns, " (", perc, "%)\n"
  ))

  # Run deconvolution
  message(paste0("Mixture deconvolution (method: ", method, ")\n"))
  results1 <- quanTIseq(sig.mat,
    mix.mat,
    scaling = mRNA,
    method = method
  )
  if ("Tregs" %in% colnames(sig.mat) && "T.cells.CD4" %in% colnames(sig.mat) && method %in% c("lsei")) {
    minTregs <- 0.02
    i <- which(colnames(sig.mat) == "T.cells.CD4")
    results2 <- quanTIseq(sig.mat[, -i],
      mix.mat,
      scaling = mRNA[-i],
      method = method
    )

    ind <- which(results1[, "Tregs"] < minTregs)
    if (length(ind) > 0) {
      results1[ind, "Tregs"] <- (results2[ind, "Tregs"] + results1[ind, "Tregs"]) / 2
      results1[ind, "T.cells.CD4"] <- pmax(0, results1[ind, "T.cells.CD4"] - (results2[ind, "Tregs"] + results1[ind, "Tregs"]) / 2)
    }
  }
  results <- results1
  results <- results / apply(results, 1, sum)

  results = data.frame(results)
  results<-cbind(rownames(results), results)
  colnames(results)[1]<-"Sample"

  message("Deconvolution sucessful!")

  return(results)
}

#' Run the quanTIseq algorithm
#'
#' Use quanTIseq to deconvolute a gene expression matrix.
#'
#' @param expression_data The gene expression information, containing the TPM
#' values for the measured features.
#' Can be provided as
#' - a simple gene expression matrix, or a data frame (with HGNC gene symbols as
#' row names and sample identifiers as column names)
#' - an `ExpressionSet` object (from the Biobase package), where the HGNC gene symbols
#' are provided in a column of the `fData` slot - that is specified by the `column`
#' parameter below
#' - a `SummarizedExperiment` object, or any of the derivative classes (e.g. DESeq2's
#' `DESeqDataSet`), in which the assay (default: "abundance") is containing the TPMs
#' as expected. Internally, `quantiseqr` handles the conversion to an object which
#' is used in the deconvolution procedure.
#' @param signature_matrix Character string, specifying the name of the signature matrix.
#' At the moment, only the original `TIL10` signature can be selected.
#' @param is_arraydata Logical value. Should be set to TRUE if the expression data
#' are originating from microarray data. For RNA-seq data, this has to be FALSE
#' (default value). If set to TRUE, the `rmgenes` parameter (see below) is set
#' to "none". 
#' @param is_tumordata Logical value. Should be set to TRUE if the expression data
#' is from tumor samples. Default: FALSE (e.g. for RNA-seq from blood samples)
#' @param scale_mRNA Logical value. If set to FALSE, it disables the correction
#' of cell-type-specific mRNA content bias. Default: TRUE
#' @param method Character string, defining the deconvolution method to be used:
#' `lsei` for constrained least squares regression, `hampel`, `huber`, or `bisquare`
#' for robust regression with Huber, Hampel, or Tukey bisquare estimators,
#' respectively. Default: `lsei`.
#' @param column Character, specifies which column in the `fData` slot (for the
#' ExpressionSet object) contains the information of the HGNC gene symbol
#' identifiers
#' @param rm_genes Character vector, specifying which genes have to be excluded
#' from the deconvolution analysis. It can be provided as
#' - a vector of gene symbols (contained in the `expression_data`)
#' - a single string among the choices of "none" (no genes are removed) and "default"
#'   (a list of genes with noisy expression RNA-seq data is removed, as explained
#'   in the quanTIseq paper).
#' Default: "default" for RNA-seq data, "none" for microarrays. 
#' @param return_se Logical value, controls the format of how the quantification
#' is returned. If providing a `SummarizedExperiment`, it can simply extend its
#' `colData` component, without the need to create a separate data frame as output.
#'
#' @details The values contained in the `expression_data` need to be provided as
#' TPM values, as this is the format also used to store the `TIL10` signature, upon
#' which quanTIseq builds to perform the immune cell type deconvolution.
#' Expression data should _not_ be provided in logarithmic scale.
#'
#' If providing the `expression_data` as a `SummarizedExperiment`/`DESeqDataSet`
#' object, it might be beneficial that this has been created via `tximport` -
#' if this is the case, the assay named "abundance" will be automatically
#' created upon importing the transcript quantification results.
#'
#' @return A data.frame containing the quantifications of the cell type proportions,
#' or alternatively, if providing `expression_data` as `SummarizedExperiment` and
#' setting `return_se` to TRUE, a `SummarizedExperiment` with the quantifications
#' included by expanding the `colData` slot of the original object
#'
#' @references
#' F. Finotello, C. Mayer, C. Plattner, G. Laschober, D. Rieder,
#' H. Hackl, A. Krogsdam, Z. Loncova, W. Posch, D. Wilflingseder, S. Sopper,
#' M. Jsselsteijn, T. P. Brouwer, D. Johnsons, Y. Xu, Y. Wang, M. E. Sanders,
#' M. V. Estrada, P. Ericsson-Gonzalez, P. Charoentong, J. Balko,
#' N. F. d. C. C. de Miranda, Z. Trajanoski.
#' "Molecular and pharmacological modulators of the tumor immune contexture
#' revealed by deconvolution of RNA-seq data".
#' Genome Medicine 2019;11(1):34. doi: 10.1186/s13073-019-0638-6.
#'
#' C. Plattner, F. Finotello, D. Rieder.
#' "Chapter Ten - Deconvoluting tumor-infiltrating immune cells from RNA-seq
#' data using quanTIseq".
#' Methods in Enzymology, 2020. doi: 10.1016/bs.mie.2019.05.056.
#'
#' @export
#'
#' @examples
#'
#' data(dataset_racle)
#' dim(dataset_racle$expr_mat)
#' res_quantiseq_run <- quantiseqr::run_quantiseq(
#'   expression_data = dataset_racle$expr_mat,
#'   signature_matrix = "TIL10",
#'   is_arraydata = FALSE,
#'   is_tumordata = TRUE,
#'   scale_mRNA = TRUE
#' )
#'
#' # using a SummarizedExperiment object
#' library("SummarizedExperiment")
#' se_racle <- SummarizedExperiment(
#'   assays = List(
#'     abundance = dataset_racle$expr_mat
#'   ),
#'   colData = DataFrame(
#'     SampleName = colnames(dataset_racle$expr_mat)
#'   )
#' )
#'
#' res_run_SE <- quantiseqr::run_quantiseq(
#'     expression_data = se_racle,
#'     signature_matrix = "TIL10",
#'     is_arraydata = FALSE,
#'     is_tumordata = TRUE,
#'     scale_mRNA = TRUE
#' )
#'
run_quantiseq <- function(expression_data,
                          signature_matrix = "TIL10",
                          is_arraydata = FALSE,
                          is_tumordata = FALSE,
                          scale_mRNA = TRUE,
                          method = "lsei",
                          column = "gene_symbol",
                          rm_genes = NULL,
                          return_se = is(expression_data, "SummarizedExperiment")) {

  stopifnot(is.logical(is_arraydata))
  stopifnot(is.logical(is_tumordata))
  stopifnot(is.logical(scale_mRNA))

  # automatically handles that a right option is passed
  stopifnot(is.character(method))
  method <- match.arg(method, c("lsei", "hampel", "huber", "bisquare"))

  ## handle case of SummarizedExperiment
  if (is(expression_data, "SummarizedExperiment")) {
    mix.mat <- se_to_matrix(expression_data)
  }

  # convert expression set to matrix, if required.
  if (is(expression_data, "ExpressionSet")) {
    mix.mat <- eset_to_matrix(expression_data, column)
  }

  if (is(expression_data, "matrix")) {
    mix.mat <- expression_data
  }

  if (is(expression_data, "data.frame")) {
    mix.mat <- as.matrix(expression_data)
  }

  if (is.integer(mix.mat)) {
    warning("Discrete values detected in the expression data, please keep in mind ",
            "that quanTIseq required the values formatted as TPM!")
  }

  # TODO: probably move the check above
  if (!is.numeric(mix.mat)) {
    stop("Expecting a matrix/an object with numeric values to be provided")
  }


  # TODO - slightly trickier to see how rmgenes should be structured
  # should be a character, at least

  message("\nRunning quanTIseq deconvolution module\n")

  # TODO: a thought: after the checks, what about printing out all
  # options and explain in brief what they should do?
  # OK for the printing, not sure about explanation (the user should be aware of them before)

  # List of genes to be discarded
  if (is.null(rm_genes) && is_arraydata == TRUE) { # For Microarrays
    rm_genes <- "none"
  } else if (is.null(rm_genes) && is_arraydata == FALSE) { # For RNA-seq
    rm_genes <- "default"
  }

  # Files
  listsig <- c("TIL10")
  if (signature_matrix %in% listsig) {
    sig.mat.file <- system.file("extdata", paste0(signature_matrix, "_signature.txt"),
      package = "quantiseqr", mustWork = TRUE
    )

    mRNA.file <- system.file("extdata", paste0(signature_matrix, "_mRNA_scaling.txt"),
      package = "quantiseqr", mustWork = TRUE
    )

    fileab <- system.file("extdata", paste0(signature_matrix, "_TCGA_aberrant_immune_genes.txt"),
      package = "quantiseqr", mustWork = TRUE
    )
    abgenes <- as.vector(read.table(fileab, header = FALSE, sep = "\t")[, 1])

    lrmgenes <- NULL
    if (length(rm_genes) == 1) {

      if (rm_genes == "default") {
        filerm <- system.file("extdata", paste0(signature_matrix, "_rmgenes.txt"),
                              package = "quantiseqr", mustWork = TRUE
        )
        lrmgenes <- as.vector(read.table(filerm, header = FALSE, sep = "\t")[, 1])

      } else if (rm_genes == "none") {
        lrmgenes <- c()

      }

    } else {
      lrmgenes <- rm_genes

    }

  } else {

    sig.mat.file <- paste0(signature_matrix, "_signature.txt")
    mRNA.file <- paste0(signature_matrix, "_mRNA_scaling.txt")

    if (!file.exists(sig.mat.file)) {
      stop("Signature matrix file not found! ",
           "quantiseqr is expecting to find a file called ", sig.mat.file)
    }
    if (!file.exists(mRNA.file)) {
      stop("Scaling info file not found! ",
           "quantiseqr is expecting to find a file called ", mRNA.file)
    }
    # TODO: would need to check that these files, if user-specified signatures are allowed
      ## exist
      ## are formatted as expected?
      ## I would can skip this for the first release
  }


  # Load signature
  sig.mat <- read.table(sig.mat.file, header = TRUE, sep = "\t", row.names = 1)
  check_signature(sig.mat, mix.mat)

  # Load normalization factors (set all to 1 if mRNAscale==FALSE)
  if (scale_mRNA) {
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
    is_arraydata, ")\n"
  ))
  mix.mat <- fixMixture(mix.mat, arrays = is_arraydata)

  # Remove noisy genes
  n1 <- nrow(sig.mat)
  sig.mat <- sig.mat[!rownames(sig.mat) %in% lrmgenes, , drop = FALSE]
  n2 <- nrow(sig.mat)
  message(paste0("Removing ", n1 - n2, " noisy genes\n"))

  # Fix tumor data
  if (is_tumordata) {
    n1 <- nrow(sig.mat)
    sig.mat <- sig.mat[!rownames(sig.mat) %in% abgenes, , drop = FALSE]
    n2 <- nrow(sig.mat)
    message(paste0("Removing ", n1 - n2, " genes with high expression in tumors\n"))
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
  results <- results / rowSums(results)
  results <- data.frame(Sample = rownames(results), results)

  message("Deconvolution successful!")

  if (is(expression_data, "SummarizedExperiment") & return_se) {
    colnames(results) <- paste0("quanTIseq_",
                                signature_matrix, "_",
                                colnames(results))
    quantiseq_coldata <- cbind(
      SummarizedExperiment::colData(expression_data),
      results[, -1]
    )
    SummarizedExperiment::colData(expression_data) <- quantiseq_coldata
    return(expression_data)
  }

  return(results)
}

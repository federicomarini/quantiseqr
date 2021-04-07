#' Helper functions for quanTIseq
#'
#' @name quantiseq_helper
#'
NULL





#' Convert a `Biobase::ExpressionSet` to a gene-expression matrix.
#'
#' @param eset `ExpressionSet`
#' @param column column name of the `fData()` table, which contains the HGNC gene symbols.
#' @return A matrix with gene symbols as rownames and sample identifiers as colnames.
#'
#' @export
#'
#' @examples
#'
#' data(dataset_racle)
#' dim(dataset_racle$expr_mat)
#'
#' library("Biobase")
#' es_racle <- ExpressionSet(assayData = dataset_racle$expr_mat)
#' featureData(es_racle)$gene_symbol <- rownames(dataset_racle$expr_mat)
#'
#' es_racle
#'
#' head(eset_to_matrix(es_racle, "gene_symbol"))
#'
eset_to_matrix <- function(eset, column) {
  expr_mat <- exprs(eset)
  rownames(expr_mat) <- fData(eset)[[column]]
  return(expr_mat)
}

#' SummarizedExperiment to matrix
#'
#' @param se A `SummarizedExperiment` object, or any of its derivates, which
#' contains the information on the TPM expression values, which are stored in a
#' specified assay slot.
#' @param assay A character string, specifying the name of the `assays` component
#' of the `se` object. Defaults to "abundance", as this is the common convention
#' used e.g. by the `tximport` package to store the values imported from the
#' transcript level quantifications
#'
#' @return A matrix object, containing the TPM values, ready to be used in the
#' framework of `quantiseqr`
#' @export
#'
#' @examples
#' library("SummarizedExperiment")
#' library("macrophage")
#' data("gse", package = "macrophage")
#' se <- gse
#'
#' # If using ENSEMBL or Gencode gene annotation, you might want to convert the row names
#' ## in this case, the gene symbols are provided as rowData information
#' rownames(se) <- rowData(se)$SYMBOL
#'
#' tpm_matrix <- se_to_matrix(se, assay = "abundance")
#'
#' ## otherwise, you can map the identifiers via
#' library("org.Hs.eg.db")
#' library("AnnotationDbi")
#' se <- gse
#' # keep the parts before the '.', used in the Gencode annotation
#' rownames(se) <- substr(rownames(se), 1, 15)
#' gene_names <- mapIds(org.Hs.eg.db,
#'                      keys = rownames(se),
#'                      column = "SYMBOL",
#'                      keytype = "ENSEMBL")
#' rownames(se) <- gene_names
#'
#' # If you require to convert the counts to TPMs by hand, you need a vector of
#' # gene lengths as well, and then run this simple function on the count matrix
#' counts_to_tpm <- function(counts, lengths) {
#'   ratio <- counts / lengths
#'   mytpm <- ratio / sum(ratio) * 1e6
#'   return(mytpm)
#' }
#' # then run via
#' # tpmdata <- counts_to_tpm(count_matrix, genelength_vector)
#'
se_to_matrix <- function(se,
                         assay = "abundance") {
  if (!is(se, "SummarizedExperiment"))
    stop("Please provide a SummarizedExperiment as input, or a derivate of this class")

  if (!(assay %in% names(assays(se))))
    stop("The specified name of the assay could not be found in the provided `se` object")

  if (assay == "counts")
    warning("Please consider quanTIseq expects expression values formatted as TPM values, ",
            "you selected an assay name which is commonly used for raw count values. \n",
            "Consider using a different name, or a function that converts the counts to TPMs - ",
            "see the example section for a simple implementation (counts and lengths are needed)")

  exp_mat <- assays(se)[[assay]]

  guessed_ensembl_gencode <- grepl(pattern = "^ENS", rownames(exp_mat))
  if (mean(guessed_ensembl_gencode) > 0.8)
    warning("Found a large majority of row names starting with 'ENS', mostly common ",
            "in the ENSEMBL/Gencode annotation schemes. quanTIseq requires you to ",
            "provide the row names as HGNC gene symbols, please see the example section ",
            "for an example of how to use annotation packages to convert to this format.") # TODO: or even throw an error? A lot would just probably not work...

  return(exp_mat)
}


#' Extract tumor immune quantifications
#'
#' Extract tumor immune quantifications from a SummarizedExperiment object,
#' previously processed with `run_quantiseqr()`
#'
#' @param se A `SummarizedExperiment` object, or any of its derivates, which
#' contains the quantifications extracted via `quantiseqr` in its `colData` slot.
#'
#' @return A data.frame, formatted as required by downstream functions
#' @export
#'
#' @examples
#' data(dataset_racle)
#' dim(dataset_racle$expr_mat)
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
#' extract_ti_from_se(res_run_SE)
#'
extract_ti_from_se <- function(se) {
  if (!is(se, "SummarizedExperiment"))
    stop("Please provide a SummarizedExperiment as input, or a derivate of this class")

  ti_cols <-
    colData(se)[, grepl(pattern = "quanTIseq", colnames(colData(se)))]

  ti_celltypes <- unlist(
    lapply(strsplit(colnames(ti_cols), split = "_"), function(arg) arg[[3]])
  )

  quanti_sig <- unlist(
    lapply(strsplit(colnames(ti_cols), split = "_"), function(arg) arg[[2]])
  )

  if (length(unique(quanti_sig)) == 1) {
    quanti_sig <- quanti_sig[1]
    message("Found quantifications for the ", quanti_sig, " signature...")
  } else {
    stop("Found mixed-up information for the signature in use, please check the ",
         "colData slot of the provided `se` object")
  }

  colnames(ti_cols) <- ti_celltypes
  ti_quant <- data.frame(
    Sample = colnames(se),
    ti_cols
  )

  return(ti_quant)
}



#' Plot the information on the tumor immune contexture
#'
#' Plot the information on the tumor immune contexture, as extracted with
#' `run_quantiseqr()`
#'
#' @param obj An object, either
#' - a `SummarizedExperiment` where the quantifications are stored
#' - a simple data.frame object, as obtained by `run_quantiseqr()`
#'
#' @return A ggplot object
#' @export
#'
#' @examples
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
#' quantiplot(res_quantiseq_run)
#' # equivalent to...
#' quantiplot(res_run_SE)
#'
quantiplot <- function(obj) {
  # if providing SE...
  if (is(obj, "SummarizedExperiment")) {
    ti_quant <- extract_ti_from_se(obj)
  } else if (is.data.frame(obj)) {
    ti_quant <- obj
  }

  # checks on the columns

  ti_mat <- t(ti_quant[, -1])
  ti_df <- as.data.frame(ti_mat)
  ti_df$cell_type <- rownames(ti_df)

  ti_df_long <- gather(ti_df, key = "sample",
                       value = "fraction", -"cell_type")

  # plot as stacked bar chart
  p <- ggplot(ti_df_long, aes_string(x = "sample", y = "fraction", fill = "cell_type")) +
    geom_bar(stat = "identity") +
    coord_flip() +
    scale_fill_brewer(palette = "Paired") +
    scale_x_discrete(limits = rev(colnames(ti_mat))) +
    theme_bw()
  return(p)
}





#' Check the signature matrix
#'
#' Checks requirements for the signature matrix, with respect to the expression
#' matrix data provided (the one on which the deconvolution algorithm needs to
#' be run)
#'
#' Performs a number of checks to ensure the compatibility of the provided
#' signature matrix in `quantiseqr`, referring also to the content of the `mix_mat`
#' mixture matrix, to be deconvoluted
#'
#' @param signature_matrix A data.frame or a matrix object, containing the
#' signature matrix
#' @param mix_mat Mixture matrix, storing the information provided as `expression_data`
#' to the main function, `run_quantiseqr()`
#'
#' @return Invisible NULL
check_signature <- function(signature_matrix, mix_mat) {
  if (!(is.data.frame(signature_matrix) | is.matrix(signature_matrix)))
    stop("The signature matrix provided is not formatted as a matrix or a data.frame")
  if (is.matrix(signature_matrix)) {
    if (!is.numeric(signature_matrix))
      stop("Signature Matrix provided not in numeric format")
  }
  if (is.data.frame(signature_matrix)) {
    if (!all(unlist(lapply(signature_matrix, is.numeric))))
      stop("Signature Data frame provided not in numeric format")
  }

  found_in_sig <- length(intersect(rownames(mix_mat), rownames(signature_matrix)))
  if (found_in_sig == 0)
    stop("No match found between signature genes and identifiers provided in the expression data!")
  if ((found_in_sig / nrow(signature_matrix)) < 0.1)
    warning("Found less than 10% of the signature genes in the provided expression data!")

  return(invisible(NULL))
}



#' Title
#'
#' @param mix.mat Matrix or data.frame with RNA-seq gene TPM or microarray
#' expression values for all samples to be deconvoluted, with gene
#' symbols as row names and sample IDs as column names. Expression
#' levels should be on non-log scale.
#' @param arrays Logical value. Should be set to TRUE if the expression data
#' are from microarrays. For RNA-seq data, this has to be FALSE (default value).
#'
#' @return The input matrix transformed to the natural scale (if needed), 
#' with fixed gene names on the rows, and TPM (for RNA-seq) or quantile (for microarrays) 
#' normalized. 
#'
#' @examples
#' 
#' data(dataset_racle)
#' mixture.fix <- fixMixture(dataset_racle$expr_mat)
fixMixture <- function(mix.mat, arrays = FALSE) {

  # Map gene names
  mix.mat <- mapGenes(mix.mat)

  # Un-log data in log2 base
  if (max(mix.mat) < 50) {
    mix.mat <- 2^mix.mat
  }

  # Quantile normalization
  if (arrays) mix.mat <- makeQN(mix.mat)

  # TPM normalization
  mix.mat <- t(t(mix.mat) * 1e6 / apply(mix.mat, 2, sum))

  return(mix.mat)
}

#' Title
#'
#' @param mix.mat Matrix or data.frame with microarray
#' gene expression values for all samples to be deconvoluted,
#' with gene symbols as row names and sample IDs as column names.
#' Expression levels should be on non-log scale.
#'
#' @return The input matrix transfromed with quantile normalization. 
#'
#' @examples
#' 
#' data(dataset_racle)
#' mixture.quantile <- makeQN(dataset_racle$expr_mat)
makeQN <- function(mix.mat) {
  cnames <- colnames(mix.mat)
  rnames <- rownames(mix.mat)
  mix.mat <- normalize.quantiles(as.matrix(mix.mat))
  colnames(mix.mat) <- cnames
  rownames(mix.mat) <- rnames
  return(mix.mat)
}



#' mapGenes
#'
#' @param mydata Matrix or data.frame with RNA-seq gene TPM or microarray
#' gene expression values for all samples to be deconvoluted,
#' with gene symbols as row names and sample IDs as column names.
#'
#' @return The input matrix with updated gene names on the rows.
#'
#' @examples
#' 
#' data(dataset_racle)
#' mixture.fixgenes <- mapGenes(dataset_racle$expr_mat)
mapGenes <- function(mydata) {
  HGNC <- read.csv(system.file("extdata", "HGNC_genenames_20170418.txt.gz", package = "quantiseqr", mustWork = TRUE),
    header = TRUE, sep = "\t"
  )

  curgenes <- rownames(mydata)
  newgenes <- rep(NA, length(curgenes))
  newgenes2 <- rep(NA, length(curgenes))
  ind <- match(curgenes, HGNC$ApprovedSymbol)

  # Current symbols and withdrawn ones
  genes.ind.notNA <- which(!is.na(ind))
  for (i in genes.ind.notNA) {
    genei <- curgenes[i]
    if (HGNC$Status[ind[i]] == "Approved") {
      newgenes[i] <- curgenes[i]
    } else if (HGNC$Status[ind[i]] == "EntryWithdrawn") {
      next
    } else {
      Wstring <- "symbolwithdrawn,see"
      newsymbol <- gsub(Wstring, "", HGNC$ApprovedName[ind[i]])
      newgenes2[i] <- newsymbol
    }
  }

  # Not found as symbols
  genes.ind.NA <- which(is.na(ind))
  for (i in genes.ind.NA) {
    genei <- curgenes[i]

    # Previous symbol?
    ind1 <- grep(genei, HGNC$PreviousSymbols)
    for (i1 in ind1) {
      array1 <- unlist(strsplit(as.character(HGNC$PreviousSymbols[i1]), ","))
      flag1 <- length(which(array1 == genei)) > 0
      if (flag1) {
        newsymbol <- as.character(HGNC$ApprovedSymbol[i1])
        newgenes2[i] <- newsymbol
      }
    }
    # Synonym?
    ind2 <- grep(genei, HGNC$Synonyms)
    for (i2 in ind2) {
      array2 <- unlist(strsplit(as.character(HGNC$Synonyms[i2]), ","))
      flag2 <- length(which(array2 == genei)) > 0
      if (flag2) {
        newsymbol <- as.character(HGNC$ApprovedSymbol[i2])
        newgenes2[i] <- newsymbol
      }
    }
  }

  newgenes2[which(newgenes2 %in% setdiff(newgenes, NA))] <- NA
  ind <- intersect(
    which(is.na(newgenes)),
    which(!is.na(newgenes2))
  )
  newgenes[ind] <- newgenes2[ind]

  mydata <- mydata[which(!is.na(newgenes)), ]
  newgenes <- newgenes[which(!is.na(newgenes))]

  # Take the median if duplicates are present
  if (any(duplicated(rownames(mydata)))) {
    message("dupe dupes, might take a little longer: TODO, faster options?")
    outdata <- aggregate(mydata, by = list(newgenes), FUN = median)
    rownames(outdata) <- outdata[, 1]
    outdata <- outdata[, -1, drop = FALSE]
    outdata <- as.data.frame(outdata)
  } else {
    message("no dupes")
    outdata <- as.data.frame(mydata)
  }

  return(outdata)
}

#' Title
#'
#' @param currsig Signature matrix to be used for deconvolution (format: genes by cell types).
#' @param currmix Mixture matrix to be deconvoluted (format: genes by samples).
#' @param scaling Logical value. If set to FALSE, it disables the correction
#' of cell-type-specific mRNA content bias. Default: TRUE
#' @param method Character string, defining the deconvolution method to be used:
#' `lsei` for constrained least squares regression, `hampel`, `huber`, or `bisquare`
#' for robust regression with Huber, Hampel, or Tukey bisquare estimators,
#' respectively. Default: `lsei`.
#'
#' @return A data.frame of cell fractions (TODO format)
#'
#' @examples
#' 
#' data(dataset_racle)
#' mixture <- dataset_racle$expr_mat
#' cellfrac <- quanTIseq(mixture, signature)
#' TODO: how to load TIL10 as signature in the example?
quanTIseq <- function(currsig, currmix, scaling = TRUE, method = "lsei") {
  method <- match.arg(method, c("lsei", "hampel", "huber", "bisquare"))

  cgenes <- intersect(rownames(currsig), rownames(currmix))
  currsig <- as.matrix(currsig[cgenes, ])
  currmix <- as.matrix(currmix[cgenes, ])

  if (method == "lsei") {

    # Run deconvolution with constrained least squares
    G <- matrix(0, ncol = ncol(currsig), nrow = ncol(currsig))
    diag(G) <- 1
    G <- rbind(G, rep(-1, ncol(G)))
    H <- c(rep(0, ncol(currsig)), -1)

    results <- apply(currmix, 2, DClsei,
      A = currsig, G = G, H = H,
      scaling = scaling
    )
  } else {

    # Run deconvolution with robust regression
    results <- apply(currmix, 2, DCrr,
      A = currsig,
      method = method, scaling = scaling
    )
  }

  if (!identical(rownames(results),colnames(currmix))) results<-t(results)
  # TODO do extensive testing

  return(results)
}

#' Title
#'
#' @param b Numeric vector containing the right-hand side of the quadratic function to be minimised.
#' @param A Numeric matrix containing the coefficients of the quadratic function to be minimised.
#' @param G Numeric matrix containing the coefficients of the inequality constraints.
#' @param H Numeric vector containing the right-hand side of the inequality constraints.
#' @param scaling TODO
#'
#' @return TODO
#'
#' @details The [limsolve::lsei()] function is used as underlying framework. Please
#' refer to that function for more details.
#'
#' @examples
#' # TODO
DClsei <- function(b, A, G, H, scaling) {
  sc <- norm(A, "2")
  A <- A / sc
  b <- b / sc

  res <- lsei(
    A = A,
    B = b,
    G = G,
    H = H,
    verbose = FALSE
  )
  est <- res$X

  est.sum <- sum(est)
  est <- est / scaling
  est <- est / sum(est) * est.sum
  est <- c(est, pmax(0, 1 - sum(est)))
  names(est)[length(est)] <- "Other"

  return(est)
}

#' Title
#'
#' @param b TODO
#' @param A TODO
#' @param method TODO
#' @param scaling TODO
#'
#' @return TODO
#'
#' @examples
#' # TODO
DCrr <- function(b, A, method, scaling) {

  # Robust regression
  m <- paste0("psi.", method)
  if (m == "psi.hampel") {
    bres <- rlm(b ~ A, psi = m, a = 1.5, b = 3.5, c = 8, maxit = 1e3)
  } else {
    bres <- rlm(b ~ A, psi = m, maxit = 1e3)
  }
  est <- bres$coefficients

  # Remove intercept
  est <- est[-1]

  # Set negative values to 0
  est[est < 0] <- 0

  # Normalize total to 1=100%
  est <- est / sum(est)

  # Scale by mRNA content
  est.sum <- sum(est)
  est <- est / scaling
  est <- est / sum(est) * est.sum

  names(est) <- gsub("^A", "", names(est))

  return(est)
}

#' Title
#'
#' @param DCres Data.frame of deconvoluted cell fractions computed with the "run_quantiseq" function, with sample identifiers as row names.
#' @param density_info Named numeric vector of total cell densities per sample. The vector names should match the sample identifiers specified in DCres.
#' derived from imaging data
#'
#' @return TODO
#'
#' @examples
#' # TODO
get_densities <- function(DCres,
                          density_info) {

  # checks:
  ## DCres should be formatted as output coming from quantiseqr
  ## density info should have the same samples included
    ## think of which format this needs to be provided

  ## TODO: check if intersection is empty?
  csbj <- intersect(rownames(DCres), names(density_info))
  density_info <- density_info[csbj]
  DCres <- DCres[csbj, , drop = FALSE]

  celldens <- data.frame(DCres, row.names=1)
  for (i in seq_len(nrow(celldens))) {
    celldens[i, ] <- celldens[i, ] * density_info
  }

  return(celldens)

}

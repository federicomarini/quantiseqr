#' Helper functions for quanTIseq
#'
#' Source code from https://github.com/FFinotello/quanTIseq
#'
#' @name quantiseq_helper
#'
NULL





#' Convert a `Biobase::ExpressionSet` to a gene-expression matrix.
#'
#' @param eset `ExpressionSet`
#' @param column column name of the `fData()` table, which contains the HGNC gene symbols.
#' @return matrix with gene symbols as rownames and sample identifiers as colnames.
#'
#' @export
eset_to_matrix <- function(eset, column) {
  expr_mat <- exprs(eset)
  # TODO: could be done also without dplyr code to reduce dependencies to minimum!
  rownames(expr_mat) <- fData(eset) %>% pull(!!column)
  return(expr_mat)
}


## TODO: potential arguments: which assay should I pick if more are available?
se_to_matrix <- function(se) {
  # TODO: define all behavior and edge cases
}


#' TODO: something similar, but for the SummarizedExperiment class?





#' TODO: even if unexported, it would not hurt to have the functions a little more
#' documented



########




#' Title
#'
#' @param mix.mat TODO
#' @param arrays TODO
#'
#' @return TODO
#'
#' @examples
#' # TODO
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
#' @param mix.mat TODO
#'
#' @return TODO
#'
#' @examples
#' # TODO
makeQN <- function(mix.mat) {
  cnames <- colnames(mix.mat)
  rnames <- rownames(mix.mat)
  mix.mat <- normalize.quantiles(as.matrix(mix.mat))
  colnames(mix.mat) <- cnames
  rownames(mix.mat) <- rnames
  return(mix.mat)
}



## TODO: maybe this can be cleverly handled with orgDb packages?
#' Title TODO
#'
#' @param mydata TODO
#'
#' @return TODO
#'
#' @examples
#' # TODO
mapGenes <- function(mydata) {
  HGNC <- read.csv(system.file("extdata", "HGNC_genenames_20170418.txt", package = "quantiseqr", mustWork = TRUE),
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

    # Previos symbol?
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
  outdata <- aggregate(mydata, by = list(newgenes), FUN = median)
  rownames(outdata) <- outdata[, 1]
  outdata <- outdata[, -1, drop = FALSE]
  outdata <- as.data.frame(outdata)


  return(outdata)
}

#' Title
#'
#' @param currsig TODO
#' @param currmix TODO
#' @param scaling TODO
#' @param method TODO
#'
#' @return TODO
#'
#' @examples
#' # TODO
quanTIseq <- function(currsig, currmix, scaling, method) {
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

  # if (nrow(results)!=ncol(currmix))
  results <- t(results)

  return(results)
}

#' Title
#'
#' @param b TODO
#' @param A TODO
#' @param G TODO
#' @param H TODO
#' @param scaling TODO
#'
#' @return TODO
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
#' @param DCres TODO
#' @param density_info TODO
#'
#' @return TODO
#'
#' @examples
#' # TODO
celldensities <- function(DCres,
                          density_info) {

  # checks:
  ## DCres should be formatted as output coming from quantiseqr
  ## density info should have the same samples included
    ## think of which format this needs to be provided



  # TODO: this one is not available?
  ## TODO: if we provide this with the package, it is good for showing which format is expected
  imageinfo <- system.file("extdata", "quantiseq", "totalcells.txt", package = "immunedeconv", mustWork = TRUE)
  csbj <- intersect(rownames(DCres), imageinfo[, 1])
  imageinfo <- imageinfo[imageinfo[, 1] %in% csbj, , drop = FALSE]
  DCres <- DCres[csbj, , drop = FALSE]

  celldens <- DCres

  for (i in 1:nrow(celldens)) {
    celldens[i, ] <- celldens[i, ] * imageinfo[i, 2]
  }

  return(celldens)
}

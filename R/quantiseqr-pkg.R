#' quantiseqr package
#'
#' pkg description
#'
#'
#' @importFrom utils read.table read.csv
#' @importFrom MASS rlm
#' @importFrom limSolve lsei
#' @importFrom preprocessCore normalize.quantiles
#' @importFrom methods is
#' @importFrom stats aggregate median
#' @importFrom Biobase exprs fData
#'
#' @name quantiseqr-pkg
#' @docType package
NULL



#' An exemplary dataset
#'
#' @details TODO
#'
#' @references TODO
#'
#' @name dataset_racle
#' @docType data
NULL


#' quanTIseq output for the simulation data of 1700 mixtures for RNA-seq data
#'
#' @details quanTIseq output for the simulation data of 1700 mixtures for RNA-seq data,
#' stored as a data.frame with 1700 rows (all the single instances of the different
#' mixtures) as returned by `run_quantiseq()`. Column names, accordingly, contain
#' the names of the component of the TIL10 signature, namely "B.cells", "Macrophages.M1",
#' "Macrophages.M2", "Monocytes", "Neutrophils", "NK.cells", "T.cells.CD4", "T.cells.CD8",
#' "Tregs", "Dendritic.cells", and "Other" (indicating for example a proxy for the
#' amount of tumor tissue).
#'
#' This can be compared (see Vignette for an example) to the ground truth information
#' on the components of the mixtures.
#'
#' @references Finotello, F., Mayer, C., Plattner, C. et al. Correction to: Molecular
#' and pharmacological modulators of the tumor immune contexture revealed by deconvolution
#' of RNA-seq data. Genome Med 11, 50 (2019). https://doi.org/10.1186/s13073-019-0655-5
#'
#' @name ti_quant_sim1700mixtures
#' @docType data
NULL

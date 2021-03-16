library("quantiseqr")

data(dataset_racle)

library("SummarizedExperiment")
se_racle <- SummarizedExperiment(
  assays = List(
    abundance = dataset_racle$expr_mat
  ),
  colData = DataFrame(
    SampleName = colnames(dataset_racle$expr_mat)
  )
)
se_racle

se_racle_fakeENS <- se_racle
rownames(se_racle_fakeENS) <- paste0("ENSG00000_fakenames_", rownames(se_racle))

library("Biobase")
es_racle <- ExpressionSet(assayData = dataset_racle$expr_mat)
featureData(es_racle)$gene_symbol <- rownames(dataset_racle$expr_mat)
es_racle



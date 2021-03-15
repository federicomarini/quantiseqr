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

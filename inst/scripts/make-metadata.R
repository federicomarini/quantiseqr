## validate with `ExperimentHubData::makeExperimentHubMetadata()`
## (above pkg directory)

Song2017_MAPKi_treatment <- data.frame(
  stringsAsFactors = FALSE,
  Title = "Song2017_MAPKi_treatment",
  Description = paste(
    "RNA-seq data, quantified as TPM, for tumor samples ",
    "pre and on MAPKi treatment.",
    "Represented as a SummarizedExperiment;",
    "derived from GEO accession number GSE75299"),
  BiocVersion = "3.13",
  Genome = NA,
  SourceType = "tar.gz",
  SourceUrl = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE75299",
  SourceVersion = "May 15 2019",
  Species = "Homo sapiens",
  TaxonomyId = "9606",
  Coordinate_1_based  = NA,
  DataProvider = "GEO",
  Maintainer = "Federico Marini <marinif@uni-mainz.de>",
  RDataClass = "SummarizedExperiment",
  DispatchClass = "Rda",
  RDataPath = "quantiseqr/Song2017_MAPKi_treatment.Rda")

# write to .csv
write.csv(Song2017_MAPKi_treatment, file = "../extdata/metadata.csv", row.names = FALSE)

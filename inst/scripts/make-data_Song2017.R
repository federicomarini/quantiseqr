# download files


# after processing with quanTIseq full pipeline, we will work with the TPM object

# this is the tabular format as plain text file, compressed
tpmdata <- read.table("tpm_GSE75299.txt.gz", header = TRUE)

dim(tpmdata)
colnames(tpmdata)

# extracting the gene symbols from the first column
tpm_genesymbols <- tpmdata$GENE
tpmdata <- as.matrix(tpmdata[, -1])
rownames(tpmdata) <- tpm_genesymbols

# constructing the SummarizedExperiment object
library(SummarizedExperiment)

# preparing something more on the colData to have the information available directly
coldata_Song2017 <- DataFrame(
  sra_id = colnames(tpmdata)
)

# would be used like...
ti_quant <- quantiseqr::run_quantiseq(
  expression_data = tpmdata,
  signature_matrix = "TIL10",
  is_arraydata = FALSE,
  is_tumordata = TRUE,
  scale_mRNA = TRUE
)


# annotating the samples
library("GEOquery")
gds <- getGEO("GSE75299")
pdata <- pData(gds[[1]])

# from the run table, also made available in quantiseqr directly
sraInfo <- read.csv(
  system.file("extdata", "SRP066571runtable.txt", package = "quantiseqr"),
  header = TRUE, sep = ",", stringsAsFactors = FALSE)

coldata_Song2017$geo_accession <-
  sraInfo[, "SampleName"][match(coldata_Song2017$sra_id, sraInfo$Run)]

# keeping the GEO info only for matched patient datasets
dim(pdata)

is_patient <- grep("^Pt", pdata[, "title"])

coldata_Song2017_full <- merge(
  x = coldata_Song2017,
  y = pdata,
  by = "geo_accession"
)

coldata_Song2017_full$is_patient <- grepl("^Pt", coldata_Song2017_full[, "title"])


Song2017_MAPKi_treatment <- SummarizedExperiment(
  assays = SimpleList(
    abundance = tpmdata
  ),
  colData = coldata_Song2017_full
)

save(Song2017_MAPKi_treatment, file = "Song2017_MAPKi_treatment.Rda")

# so that one can use it in quantiseqr with...
Song2017_MAPKi_treatment_tiquant <- quantiseqr::run_quantiseq(
  expression_data = Song2017_MAPKi_treatment,
  signature_matrix = "TIL10",
  is_arraydata = FALSE,
  is_tumordata = TRUE,
  scale_mRNA = TRUE
)

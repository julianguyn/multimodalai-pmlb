# load libraries
suppressPackageStartupMessages({
    library(data.table)
    library(readxl)
    library(maftools)
    library(factoextra)
    library(ggplot2)
    library(ggpubr)
})

source("workflow/scripts/utils.R")
source("workflow/scripts/palettes.R")
set.seed(101)

###########################################################
# Load in data
###########################################################

RAWDATA <- "data/rawdata/omics/"

# read in omics data
atc <- fread(paste0(RAWDATA, "atac/ATACseq_PMLB_multimodal_20260310.tsv"), data.table = FALSE)
rna <- fread(paste0(RAWDATA, "rnaseq/RNAseq_TPM_PMLB_multimodal_20260310.tsv"), data.table = FALSE)
cnv <- fread(paste0(RAWDATA, "cnv/CNV_PMLB_multimodal_20260310.tsv"), data.table = FALSE)
maf <- read.maf(paste0(RAWDATA, "mutation/mutation_PMLB_multimodal_20260310.maf"))

# read in atac peaks
peaks <- fread(paste0(RAWDATA, "atac/All_samples_consensus_peak.bed"))

# read in metadata
meta <- read_excel("metadata/PMLB_panOrganoid_multimodal_metadata_20260310.xlsx", sheet = 1) |> 
    as.data.frame()
meta <- meta[complete.cases(meta$doubling_rate),]
meta$sex <- ifelse(meta$sex == "F", "Female", "Male")

# load in snf clusters (from James)
snf2 <- read.csv("data/rawdata/snf_clusters_2.csv")
snf6 <- read.csv("data/rawdata/snf_clusters_6.csv")

###########################################################
# Format omics matrices
###########################################################

# format ATAC matrix
rownames(atc) <- paste(peaks$V1, peaks$V2, peaks$V3, sep = ":")
atc$peak <- NULL

# format RNA matrix
rownames(rna) <- rna$gene_id
rna$gene_id <- NULL

# format CNV matrix
rownames(cnv) <- cnv$Hugo_Symbol
cnv$Hugo_Symbol <- NULL

# format mutation matrix
mut <- mutCountMatrix(
  maf,
  includeSyn = FALSE,
  countOnly = NULL,
  removeNonMutated = TRUE
) |> as.data.frame()

# remove sample that failed QC
to_remove <- "PHLC0362.TXO"
meta <- meta[meta$PMLB_organoidID != "PHLC0362.TXO", ]

# keep only 65 samples of interest
atc <- atc[,match(meta$PMLB_organoidID, colnames(atc))]
rna <- rna[,match(meta$PMLB_organoidID, colnames(rna))]
cnv <- cnv[,match(meta$PMLB_organoidID, colnames(cnv))]
mut <- mut[,match(meta$PMLB_organoidID, colnames(mut))]
snf2 <- snf2[match(meta$PMLB_organoidID, snf2$Organoid.ID),]
snf6 <- snf6[match(meta$PMLB_organoidID, snf6$Organoid.ID),]

###########################################################
# Save data inputs for modeling
###########################################################

# R session save
save(meta, atc, rna, cnv, mut, file = "data/procdata/RData/omics_inputs.RData")

# files for model inputs
write.csv(atc, file = "data/procdata/files/atac.csv", quote = FALSE, row.names = TRUE)
write.csv(rna, file = "data/procdata/files/rna.csv", quote = FALSE, row.names = TRUE)
write.csv(cnv, file = "data/procdata/files/cnv.csv", quote = FALSE, row.names = TRUE)
write.csv(mut, file = "data/procdata/files/mut.csv", quote = FALSE, row.names = TRUE)
write.csv(meta, file = "data/procdata/files/meta.csv", quote = FALSE, row.names = FALSE)

###########################################################
# Data exploration: PCA Omics
###########################################################

# perform PCA
plot_panel(atc, meta, "ATAC")
plot_panel(rna, meta, "RNA")
plot_panel(cnv, meta, "CNV")
plot_panel(mut, meta, "MUT")

###########################################################
# Data exploration: Doubling Time
###########################################################

# add snf cluster assignments
meta$SNF2 <- factor(snf2$Cluster[match(meta$PMLB_organoidID, snf2$Organoid.ID)])
meta$SNF6 <- factor(snf6$Cluster[match(meta$PMLB_organoidID, snf6$Organoid.ID)])

plot_doubling("sex", 5)
plot_doubling("organoid_sample_class", 5)
plot_doubling("primary_tumor_site", 6)
plot_doubling("SNF2", 5)
plot_doubling("SNF6", 6)

# plot doubling time by snf cluster coloured by cancer type
plot_snf_clusters(snf2, meta, "SNF2")
plot_snf_clusters(snf6, meta, "SNF6")

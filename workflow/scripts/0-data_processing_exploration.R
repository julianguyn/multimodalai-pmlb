# load libraries
suppressPackageStartupMessages({
    library(data.table)
    library(dplyr)
    library(readxl)
    library(maftools)
    library(factoextra)
    library(ggplot2)
    library(ggpubr)
    library(GenomicRanges)
    library(ChIPseeker)
    library(TxDb.Hsapiens.UCSC.hg38.knownGene)
    library(biomaRt)
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
# Map to common genes
###########################################################

# get peak coordinates
coords <- do.call(rbind, strsplit(rownames(atc), ":")) |> as.data.frame()
gr <- GRanges(seqnames = coords$V1, ranges = IRanges(as.numeric(coords$V2), as.numeric(coords$V3)))

# annotate peaks
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
anno <- annotatePeak(gr, TxDb=txdb, annoDb="org.Hs.eg.db")

p <- ggplot(anno@annoStat, aes(fill=Feature, y=Frequency, x="PDOs")) + 
    geom_bar(position="fill", stat="identity", color = "black") +
    scale_fill_manual(values = genefeat_pal) +
    theme_minimal() +
    labs(y = "Percentage (%)", x = "")
ggsave("data/results/figures/0-DataExploration/atac_peakanno.png", plot=p, width = 4, height = 5)

# keep promoter regions
anno_promoter <- as.data.frame(anno)[grep("Promoter", as.data.frame(anno)$annotation), ]
anno_promoter$promoter_peaks <- paste(anno_promoter$seqnames, anno_promoter$start, anno_promoter$end, sep = ":")
atc_promoter <- atc[rownames(atc) %in% anno_promoter$promoter_peaks,] #68929 peaks

# map to genes
gene_map <- anno_promoter$SYMBOL[match(rownames(atc_promoter), anno_promoter$promoter_peaks)]

#rownames(atc_promoter) <- anno_promoter$SYMBOL[match(rownames(atc_promoter), anno_promoter$promoter_peaks)]
atc_gene <- as.data.frame(atc_promoter) %>%
  mutate(Gene = gene_map) %>%
  group_by(Gene) %>%
  summarise(across(everything(), mean, na.rm = TRUE))
atc_gene <- as.data.frame(atc_gene)
atc_gene <- atc_gene[!is.na(atc_gene$Gene), ]
rownames(atc_gene) <- atc_gene$Gene
atc_gene$Gene <- NULL
write.csv(atc_gene, file = "data/procdata/files/atac_gene.csv", quote = FALSE, row.names = TRUE)

# map RNA ensembl to gene names
mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
mapping <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = sub("\\..*", "", rownames(rna)),
  mart = mart
)
tt <- as.data.frame(rna)
tt$ensembl_gene_id <- sub("\\..*", "", rownames(tt))
rna_annot <- left_join(tt, mapping, by = "ensembl_gene_id")
rna_annot <- rna_annot %>% filter(hgnc_symbol != "")

# average expression of overlapping gene symbols
rna_gene <- rna_annot %>%
  group_by(hgnc_symbol) %>%
  summarise(across(where(is.numeric), mean)) %>%
  as.data.frame()
rownames(rna_gene) <- avg$hgnc_symbol
rna_gene$hgnc_symbol <- NULL
write.csv(rna_gene, file = "data/procdata/files/rna_gene.csv", quote = FALSE, row.names = TRUE)

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

plot_doubling("sex", 5)
plot_doubling("organoid_sample_class", 5)
plot_doubling("primary_tumor_site", 6)

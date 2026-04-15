# check RNA batches of TCGA RNA counts

# dirs to create:
# data/rawdata/TCGA-rna
# data/results/figures/TCGA-rna

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(dplyr)
})

INDIR = "data/rawdata/TCGA-rna/"
OUTDIR = "data/results/figures/TCGA-rna/"

# load data
files <- list.files(INDIR)
rna <- fread(paste0(INDIR, files[1]))
genes <- rna$attrib_name
rna <- rna[,-1]
rna <- as.data.frame(t(rna))
colnames(rna) <- genes
cohort <- c(rep("BRCA", nrow(rna)))

for (file in files[-1]) {
    cohort_name <- sub("Human__TCGA_", "", sub("__UNC__RNAseq__HiSeq_RNA__01_28_2016__BI__Gene__Firehose_RSEM_log2.cct", "", file))
    cat("---Starting", cohort_name, "\n")

    df <- fread(paste0(INDIR, file))
    print(table(df$attrib_name %in% tt))
    genes <- df$attrib_name
    df <- df[,-1]
    df <- as.data.frame(t(df))
    colnames(df) <- genes
    cohort <- c(cohort, rep(cohort_name, nrow(df)))
    rna <- bind_rows(rna, df)
}
cat("done\n\n")

# perform pca
cat("---Converting NAs to 0\n")
rna[is.na(rna)] <- 0

cat("---Starting PCA\n")
pca_res <- prcomp(rna)

save(pca_res, cohort, file = paste0(OUTDIR, "pca_res.RData"))

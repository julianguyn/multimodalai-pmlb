# archived script for lasso feature seleection plots

# load libraries
suppressPackageStartupMessages({
    library(ggplot2)
    library(ggpubr)
    library(patchwork)
})

source("workflow/scripts/utils.R")
source("workflow/scripts/palettes.R")
set.seed(101)

###########################################################
# Early fusion
###########################################################

OUTDIR <- "5-CommonGene2"
MODEL <- "rf"

# load in unimodal-1
pattern <- "_folds.csv"
path <- paste0("data/results/data/", OUTDIR, "/unimodal-1")
unimodal1 <- list.files(path, pattern = pattern, full.names = TRUE)

model_res <- data.frame(matrix(nrow=0, ncol=8))

# unimodal model outputs
for (file in unimodal1) {
    anno <- sub(paste0(path, "/"), "", sub(pattern, "", file))
    modality <- sub("_.*", "", anno)
    model <- sub(".*_", "", anno)
    df <- read.csv(file)
    model_res <- rbind(
        model_res,
        data.frame(
            approach = "Unimodal-1",
            modality = modality,
            model = model,
            spearman = df$spearman,
            pearson = df$pearson,
            rmse = df$rmse,
            mae = df$mae
        )
    )
}

# load in unimodal-2
path <- paste0("data/results/data/", OUTDIR, "/unimodal-2")
unimodal2 <- list.files(path, pattern = pattern, full.names = TRUE)
for (file in unimodal2) {
    anno <- sub(paste0(path, "/"), "", sub(pattern, "", file))
    modality <- sub("_.*", "", anno)
    model <- sub(".*_", "", anno)
    df <- read.csv(file)
    model_res <- rbind(
        model_res,
        data.frame(
            approach = "Unimodal-2",
            modality = modality,
            model = model,
            spearman = df$spearman,
            pearson = df$pearson,
            rmse = df$rmse,
            mae = df$mae
        )
    )
}

# load in EF1
SUBDIR <- "ef-1"
filepath <- paste0("data/results/data/", OUTDIR, "/", SUBDIR, "/", MODEL, "_folds.csv")
ef1 <- read.csv(filepath)[,c("spearman", "pearson", "rmse", "mae")]
ef1$approach <- "EF-1"
ef1$modality <- "EF-1"
ef1$model <- MODEL
ef1 <- ef1[,match(colnames(model_res), colnames(ef1))]

# load in EF2
SUBDIR <- "ef-2"
filepath <- paste0("data/results/data/", OUTDIR, "/", SUBDIR, "/", MODEL, "_folds.csv")
ef2 <- read.csv(filepath)[,c("spearman", "pearson", "rmse", "mae")]
ef2$approach <- "EF-2"
ef2$modality <- "EF-2"
ef2$model <- MODEL
ef2 <- ef2[,match(colnames(model_res), colnames(ef2))]

model_res <- rbind(model_res, ef1, ef2)

# load in late fusion (unimodal-1)
SUBDIR <- "unimodal-1"
filepath <- paste0("data/results/data/", OUTDIR, "/", SUBDIR, "/", MODEL, "_late_fusion_average.csv")
df <- read.csv(filepath)[,c("spearman", "pearson", "rmse", "mae")]
df$approach <- "LF-1-average"
df$modality <- "LF-1-average"
df$model <- MODEL
df <- df[,match(colnames(model_res), colnames(df))]
model_res <- rbind(model_res, df)

filepath <- paste0("data/results/data/", OUTDIR, "/", SUBDIR, "/", MODEL, "_late_fusion_weighted.csv")
df <- read.csv(filepath)[,c("spearman", "pearson", "rmse", "mae")]
df$approach <- "LF-1-weighted"
df$modality <- "LF-1-weighted"
df$model <- MODEL
df <- df[,match(colnames(model_res), colnames(df))]
model_res <- rbind(model_res, df)

# load in late fusion (unimodal-2)
SUBDIR <- "unimodal-2"
filepath <- paste0("data/results/data/", OUTDIR, "/", SUBDIR, "/", MODEL, "_late_fusion_average.csv")
df <- read.csv(filepath)[,c("spearman", "pearson", "rmse", "mae")]
df$approach <- "LF-2-average"
df$modality <- "LF-2-average"
df$model <- MODEL
df <- df[,match(colnames(model_res), colnames(df))]
model_res <- rbind(model_res, df)

filepath <- paste0("data/results/data/", OUTDIR, "/", SUBDIR, "/", MODEL, "_late_fusion_weighted.csv")
df <- read.csv(filepath)[,c("spearman", "pearson", "rmse", "mae")]
df$approach <- "LF-2-weighted"
df$modality <- "LF-2-weighted"
df$model <- MODEL
df <- df[,match(colnames(model_res), colnames(df))]
model_res <- rbind(model_res, df)


###########################################################
# Plot model outputs
###########################################################

p <- ggplot(model_res, aes(x = modality, y = spearman)) +
    geom_boxplot() +
    geom_jitter(shape = 21, alpha = 0.8, position = position_jitterdodge(jitter.width = 0.2)) +
    #scale_fill_manual("Model", values = model_pal) +
    theme_minimal() +
    theme(
        panel.border = element_rect(),
        legend.key.size = unit(0.7, 'cm')
    ) +
    labs(x = "Modalities", y = "Spearman Correlation")

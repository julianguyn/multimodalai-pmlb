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
# Load in model final data
###########################################################

pattern <- "_folds.csv"
uni_data <- "data/results/data/1-Unimodal-OmicOnly"
multi_data <- "data/results/data/2-Multimodal-OmicOnly"

unimodal <- list.files(
    uni_data,
    pattern = pattern,
    full.names = TRUE
)

multimodal <- list.files(
    multi_data,
    pattern = pattern,
    full.names = TRUE
)

###########################################################
# Parse model final outputs
###########################################################

model_res <- data.frame(matrix(nrow=0, ncol=8))

# unimodal model outputs
for (file in unimodal) {
    anno <- sub(paste0(uni_data, "/"), "", sub(pattern, "", file))
    modality <- sub("_.*", "", anno)
    model <- sub(".*_", "", anno)
    df <- read.csv(file)
    model_res <- rbind(
        model_res,
        data.frame(
            approach = "Unimodal",
            modality = modality,
            model = model,
            fold = df$fold,
            spearman = df$spearman,
            pearson = df$pearson,
            rmse = df$rmse,
            mae = df$mae
        )
    )
}

# multimodal model outputs
for (file in multimodal) {
    anno <- sub(paste0(multi_data, "/"), "", sub(pattern, "", file))
    model <- sub(".*_", "", anno)
    modality <- sub(paste0("_", model), "", anno)
    df <- read.csv(file)
    model_res <- rbind(
        model_res,
        data.frame(
            approach = "Multimodal",
            modality = modality,
            model = model,
            fold = df$fold,
            spearman = df$spearman,
            pearson = df$pearson,
            rmse = df$rmse,
            mae = df$mae
        )
    )
}

###########################################################
# Format levels
###########################################################

# modality levels
modality_levels <- c(
    "atac\nrna\ncnv\nmut",
    "atac\ncnv\nmut", "atac\nrna\ncnv", "atac\nrna\nmut", "rna\ncnv\nmut",
    "atac\nrna", "atac\ncnv", "atac\nmut", "rna\ncnv", "rna\nmut", "cnv\nmut",
    "atac", "rna", "cnv", "mut"
)

# model levels
model_labs <- c("Elastic Net", "Random Forest", "LASSO")

# format modalities and models
model_res$modality[model_res$modality == "all_modalities"] <- "atac_rna_cnv_mut"
model_res$modality <- gsub("_", "\n", model_res$modality)
model_res$modality <- factor(model_res$modality, levels = modality_levels)
model_res$model <- factor(model_res$model, levels = c("en", "rf", "lasso"), labels = model_labs)

###########################################################
# Model performance comparison
###########################################################

toPlot <- model_res[model_res$modality %in% c("atac\nrna\ncnv\nmut", "atac\nrna", "atac", "rna", "cnv", "mut"), ]

p <- ggplot(toPlot, aes(x = modality, y = spearman, fill = model)) +
    geom_boxplot() +
    geom_jitter(shape = 21, alpha = 0.8, position = position_jitterdodge(jitter.width = 0.2)) +
    scale_fill_manual("Model", values = model_pal) +
    theme_minimal() +
    theme(
        panel.border = element_rect(),
        legend.key.size = unit(0.7, 'cm')
    ) +
    labs(x = "Modalities", y = "Spearman Correlation")

filename <- paste0("data/results/figures/12-OmicOnly/model_comparisons.png")
cat("Saving figure to", filename, "\n")
png(filename, width = 8, height = 4, res = 600, units = "in")
print(p)
dev.off()

###########################################################
# Format data for plotting
###########################################################

# build long dataframe for the tile plot
modality_df <- expand.grid(
    modality = modality_levels,
    omic = c("atac", "rna", "cnv", "mut"),
    stringsAsFactors = FALSE
) %>%
  mutate(
    present = mapply(function(mod, omic) grepl(omic, mod), modality, omic),
    modality = factor(modality, levels = modality_levels),
    omic = factor(omic, levels = rev(c("atac", "rna", "cnv", "mut")))
  )

###########################################################
# Plot model outputs
###########################################################

# modality tile plot
p1 <- ggplot(modality_df, aes(x = modality, y = omic, fill = present)) +
    geom_tile(color = "black") +
    scale_fill_manual(values = c("TRUE" = col, "FALSE" = "white")) +
    theme_minimal() +
    theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none"
    ) +
    labs(y = NULL, x = NULL)

# spearman
p2 <- ggplot(model_res, aes(x = modality, y = spearman, fill = approach)) +
    geom_boxplot() +
    geom_jitter(shape = 21, alpha = 0.8, width = 0.2) +
    scale_fill_manual(values = c(pale_blue, dark_blue)) +
    theme_minimal() +
    theme(
        panel.border = element_rect(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none"

    ) +
    labs(x = "Modalities", y = "Spearman Correlation")

# pearson
p3 <- ggplot(model_res, aes(x = modality, y = pearson, fill = approach)) +
    geom_boxplot() +
    geom_jitter(shape = 21, alpha = 0.8, width = 0.2) +
    scale_fill_manual(values = c(pale_blue, dark_blue)) +
    theme_minimal() +
    theme(
        panel.border = element_rect(),
        legend.position = "none") +
    labs(x = "Modalities", y = "Pearson Correlation")

# rmse
p4 <- ggplot(model_res, aes(x = modality, y = rmse, fill = approach)) +
    geom_boxplot() +
    geom_jitter(shape = 21, alpha = 0.8, width = 0.2) +
    scale_fill_manual(values = c(pale_blue, dark_blue)) +
    theme_minimal() +
    theme(
        panel.border = element_rect(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none"

    ) +
    labs(x = "Modalities", y = "Root Mean Square Error")

# mae
p5 <- ggplot(model_res, aes(x = modality, y = mae, fill = approach)) +
    geom_boxplot() +
    geom_jitter(shape = 21, alpha = 0.8, width = 0.2) +
    scale_fill_manual(values = c(pale_blue, dark_blue)) +
    theme_minimal() +
    theme(
        panel.border = element_rect(),
        legend.position = "none"
    ) +
    labs(x = "Modalities", y = "Mean Absolute Error")

# correlations plot
filename <- paste0("data/results/figures/12-OmicOnly/correlations_", model, ".png")
cat("Saving figure to", filename, "\n")
png(filename, width = 7, height = 5.5, res = 600, units = "in")
p1 / p2 / p3  + plot_layout(heights = c(1, 2.5, 2.75))
dev.off()

# errors plot
filename <- paste0("data/results/figures/12-OmicOnly/errors_", model, ".png")
cat("Saving figure to", filename, "\n")
png(filename, width = 7, height = 5.5, res = 600, units = "in")
p1 / p4 / p5  + plot_layout(heights = c(1, 2.5, 2.75))
dev.off()

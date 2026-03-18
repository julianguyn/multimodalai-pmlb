#' Helper function to plot PCs
#' 
plot_PCs <- function(res.pca, meta, var) {

  # format dataframe
  toPlot <- as.data.frame(res.pca$x[,1:2]) # first two PCs
  toPlot$var <- meta[[var]][match(rownames(toPlot), meta$PMLB_organoidID)]

  # get label
  label <- switch(
    var,
    sex = "Sex                  ",
    primary_tumor_site = "Tissue Type",
    organoid_sample_class = "Organoid Class"
  )

  # make plot
  p <- ggplot(toPlot, aes(x = PC1, y = PC2, fill = var)) +
    geom_point(shape = 21, size = 3) +
    theme_minimal() +
    theme(
      panel.border = element_rect()
    ) +
    labs(fill = label)

  return(p)
}

#' Helper function to compile plots
#' 
plot_panel <- function(mat, meta, label) {

  # run PCA
  mat <- t(mat) |> as.data.frame()
  if (label != "MUT") {
    res.pca <- prcomp(mat, scale = TRUE)
  } else {
    rm <- colSums(mat)[which(colSums(mat) == 0)]
    mat <- mat[,-which(colnames(mat) %in% names(rm))]
    res.pca <- prcomp(mat, scale = FALSE)
  }

  # get eigenvalues
  p1 <- fviz_eig(res.pca)

  # plot PCs
  p2 <- plot_PCs(res.pca, meta, "sex") + scale_fill_manual(values = sex_pal)
  p3 <- plot_PCs(res.pca, meta, "primary_tumor_site") + scale_fill_manual(values = tissue_pal)
  p4 <- plot_PCs(res.pca, meta, "organoid_sample_class") + scale_fill_manual(values = orgclass_pal)

  filename <- paste0("data/results/figures/0-DataExploration/pca_", label, ".png")
  cat("Saving figure to", filename, "\n")
  png(filename, width = 9, height = 5.5, res = 600, units = "in")
  print(ggarrange(p1, p2, p3, p4, nrow=2, ncol=2))
  dev.off()

}

#' Function to plot doubling rate vs variables
#'
plot_doubling <- function(var, w) { 

  # get label
  label <- switch(
    var,
    sex = "Sex",
    primary_tumor_site = "Tissue Type",
    organoid_sample_class = "Organoid Class",
    SNF2 = "SNF2",
    SNF6 = "SNF6"
  )

  p <- ggplot(meta, aes(x = .data[[var]], y = doubling_rate)) +
    geom_boxplot(fill = "#9AA899") + 
    geom_jitter(width = 0.2, alpha = 0.6) +
    theme_minimal() +
    theme(
      panel.border = element_rect()
    ) +
    labs(y = "Doubling Rate", x = label)
  
  filename <- paste0("data/results/figures/0-DataExploration/DT_", sub(" ", "_", label), ".png")
  cat("Saving figure to", filename, "\n")
  png(filename, width = w, height = 4, res = 600, units = "in")
  print(p)
  dev.off()

}


#' Plot doubling rate vs snf clusters + cancer type
#' 
plot_snf_clusters <- function(snf, meta, label) {

  w <- switch(label, SNF2 = 7, SNF6 = 9)
  h <- switch(label, SNF2 = 3, SNF6 = 5.5)

  snf$doubling_rate <- meta$doubling_rate[match(snf$Organoid.ID, meta$PMLB_organoidID)]
  snf$Cluster <- factor(snf$Cluster, levels = unique(snf$Cluster))

  p <- ggplot(snf, aes(x = Primary.Tissue, y = doubling_rate, fill = Primary.Tissue)) +
    geom_boxplot() +
    geom_jitter(
      aes(fill = Primary.Tissue),
      position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.75),
      shape = 21
    ) +
    facet_wrap(~Cluster, scales = "free_x") +
    scale_fill_manual("Primary Tissue", values = tissue_pal) +
    theme_minimal() +
    theme(
      panel.border = element_rect(),
      legend.key.size = unit(0.7, 'cm'),
      axis.text.x = element_blank()
    ) +
    labs(y = "Doubling Rate", x = "")

  filename <- paste0("data/results/figures/0-DataExploration/snf_", label, ".png")
  cat("Saving figure to", filename, "\n")
  png(filename, width = w, height = h, res = 600, units = "in")
  print(p)
  dev.off()
}
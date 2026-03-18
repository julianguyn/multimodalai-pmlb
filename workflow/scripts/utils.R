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
  p2 <- plot_PCs(res.pca, meta, "sex")
  p3 <- plot_PCs(res.pca, meta, "primary_tumor_site")
  p4 <- plot_PCs(res.pca, meta, "organoid_sample_class")

  filename <- paste0("data/results/figures/0-DataExploration/pca_", label, ".png")
  cat("Saving figure to", filename, "\n")
  png(filename, width = 9, height = 5.5, res = 600, units = "in")
  print(ggarrange(p1, p2, p3, p4, nrow=2, ncol=2))
  dev.off()

}

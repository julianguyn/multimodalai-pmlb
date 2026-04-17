# correlation matrix

suppressPackageStartupMessages({
    library(data.table)
    library(tidyverse)
    library(ComplexHeatmap)
    library(circlize)
})


toPlot <- read.csv("data/rawdata/correlations.csv")
toPlot$View.2[toPlot$View.2 == "Methyl."] <- "Methylation"
toPlot$View.1[toPlot$View.1 == "Methyl."] <- "Methylation"

toPlot <- toPlot %>%
  select(View.1, View.2, Correlation) %>%
  pivot_wider(names_from = View.2, values_from = Correlation) %>%
  column_to_rownames("View.1") %>%
  as.matrix()

toPlot[upper.tri(toPlot)] <- NA

Heatmap(
  toPlot,
  name = "Correlation",
  na_col = "white",
  col = colorRamp2(c(0.7, 0.85, 1), c("#66c2a5", "#fc8d62", "#8da0cb")),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  cell_fun = function(j, i, x, y, width, height, fill) {
    val <- toPlot[i, j]
    if (!is.na(val)) {
        grid.text(sprintf("%.2f", val), x, y, gp = gpar(fontsize = 8))
    }
    }
)

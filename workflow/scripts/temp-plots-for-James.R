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

# get rid of diagonal
for (i in 1:5) {
  toPlot[i,i] <- NA
}


Heatmap(
  toPlot,
  name = "Correlation",
  na_col = "white",
  col = colorRamp2(c(0.75, 0.85, 1), c("#66c2a5", "#fc8d62", "#8da0cb")),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
)

# color option 1: blue white pueple
Heatmap(
  toPlot,
  name = "Correlation",
  na_col = "white",
  col = colorRamp2(c(0.75, 0.85, 1), c("#7F7CAF", "#9FB4C7", "#28587B")),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
)

# color option 2: red orange yellow
Heatmap(
  toPlot,
  name = "Correlation",
  na_col = "white",
  col = colorRamp2(c(0.75, 0.85, 1), c("#F7B538", "#DB7C26", "#C32F27")),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
)

# color option 3: orange and browns
Heatmap(
  toPlot,
  name = "Correlation",
  na_col = "white",
  col = colorRamp2(c(0.75, 0.85, 1), c("#A57F60", "#E8AE68", "#DB5A42")),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
)

# color option 4: pink and blue
Heatmap(
  toPlot,
  name = "Correlation",
  na_col = "white",
  col = colorRamp2(c(0.75, 0.85, 1), c("#93BFB2", "#C5E0D8", "#CEABB1")),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
)

# color option 5: teal and taupe
Heatmap(
  toPlot,
  name = "Correlation",
  na_col = "white",
  col = colorRamp2(c(0.75, 0.85, 1), c("#A18276", "#FCDFA6", "#AAC0AA")),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
)

# color option 5: teal and taupe again ish
Heatmap(
  toPlot,
  name = "Correlation",
  na_col = "white",
  col = colorRamp2(c(0.75, 0.85, 1), c("#6D5959", "#9DCBBA", "#ABEBD2")),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
)

# color option 6: YELLOW blue
Heatmap(
  toPlot,
  name = "Correlation",
  na_col = "white",
  col = colorRamp2(c(0.75, 0.85, 1), c("#294C60", "#ADB6C4", "#FFC49B")),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
)

# color option 6: YELLOW blue
Heatmap(
  toPlot,
  name = "Correlation",
  na_col = "white",
  col = colorRamp2(c(0.75, 0.85, 1), c("#565676", "#AEADF0", "#A76571")),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
)

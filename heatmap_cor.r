#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

option_list <- list(
  make_option(c("-i", "--input"), type = "character",
              help = "Input correlation matrix TSV file", metavar = "FILE"),
  make_option(c("-o", "--output"), type = "character",
              help = "Output PDF file", metavar = "FILE")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$input) || is.null(opt$output)) {
  stop("Usage: Rscript plot_germline_qc_heatmap.R -i export_corr_matrix.tsv -o germline_qc.pdf")
}

## read in data
mat <- read.table(opt$input, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
row_names <- rownames(mat)
col_names <- colnames(mat)
mat <- as.matrix(mat)

## plot heatmap
pdf(opt$output, width = 20, height = 20, useDingbats = FALSE)

p <- Heatmap(
  mat,
  col = colorRamp2(c(0, 0.5, 1), c("dodgerblue", "white", "red3")),
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 8),
  heatmap_legend_param = list(
    title = "correlation",
    title_position = "topleft",
    legend_height = unit(5, "cm"),
    legend_direction = "vertical"
  ),
  row_order = row_names,
  column_order = col_names
)

draw(p)
dev.off()

cat("Heatmap saved to:", opt$output, "\n")

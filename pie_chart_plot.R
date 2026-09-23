#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(ggplot2))

frac <- read.csv("results/03_estimated_fractions.csv", check.names = FALSE)

if (any(c("celltype", "cell_type") %in% names(frac))) {
  # long format: one row per sample x celltype
  ct_col  <- intersect(c("celltype", "cell_type"), names(frac))[1]
  id_cols <- c("sample_id", "case_id", "patient_id")
  val_col <- setdiff(names(frac), c(ct_col, id_cols))[1]
  df <- aggregate(frac[[val_col]], by = list(celltype = frac[[ct_col]]), FUN = mean, na.rm = TRUE)
  names(df) <- c("celltype", "prop")
} else {
  # wide format: one column per cell type
  id_cols   <- intersect(c("sample_id", "case_id", "patient_id", "X"), names(frac))
  type_cols <- setdiff(names(frac), id_cols)
  mean_props <- colMeans(frac[, type_cols], na.rm = TRUE)
  df <- data.frame(celltype = names(mean_props), prop = as.numeric(mean_props))
}

df$prop <- df$prop / sum(df$prop)   # renormalise to sum to 1
write.csv(df, "results/tcga_mean_celltype_proportions.csv", row.names = FALSE)
print(df)

p <- ggplot(df, aes(x = "", y = prop, fill = celltype)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  geom_text(aes(label = paste0(round(prop * 100, 1), "%")),
            position = position_stack(vjust = 0.5), size = 4) +
  labs(title = "Mean estimated cell-type composition — TCGA-OV bulk (deconvolved)",
       fill = "Cell type", x = NULL, y = NULL) +
  theme_void(base_size = 13)

ggsave("figs/tcga_celltype_pie.png", p, width = 6, height = 6, dpi = 300)
message("Saved figs/tcga_celltype_pie.png and results/tcga_mean_celltype_proportions.csv")
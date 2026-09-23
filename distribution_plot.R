#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(ggplot2))

frac <- read.csv("results/03_estimated_fractions.csv", check.names = FALSE)

if (any(c("celltype", "cell_type") %in% names(frac))) {
  # long format: one row per sample x celltype
  ct_col  <- intersect(c("celltype", "cell_type"), names(frac))[1]
  id_cols <- c("sample_id", "case_id", "patient_id")
  val_col <- setdiff(names(frac), c(ct_col, id_cols))[1]
  df <- data.frame(celltype = frac[[ct_col]], prop = frac[[val_col]])
} else {
  # wide format: one column per cell type -> reshape to long
  id_cols   <- intersect(c("sample_id", "case_id", "patient_id", "X"), names(frac))
  type_cols <- setdiff(names(frac), id_cols)
  df <- do.call(rbind, lapply(type_cols, function(ct) {
    data.frame(celltype = ct, prop = frac[[ct]])
  }))
}

# quick summary table alongside the plot
summary_df <- aggregate(prop ~ celltype, df, function(x)
  c(mean = mean(x, na.rm = TRUE), median = median(x, na.rm = TRUE),
    sd = sd(x, na.rm = TRUE), n = sum(!is.na(x))))
write.csv(df, "results/tcga_celltype_scores_long.csv", row.names = FALSE)
print(summary_df)

p <- ggplot(df, aes(x = celltype, y = prop, fill = celltype)) +
  geom_boxplot(outlier.alpha = 0.4) +
  geom_jitter(width = 0.15, alpha = 0.15, size = 0.8) +
  labs(title = "Distribution of estimated cell-type fractions across TCGA-OV samples (n=425)",
       x = "Cell type", y = "Estimated fraction") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none")

ggsave("figs/tcga_celltype_score_distributions.png", p, width = 8, height = 5.5, dpi = 300)
message("Saved figs/tcga_celltype_score_distributions.png and results/tcga_celltype_scores_long.csv")
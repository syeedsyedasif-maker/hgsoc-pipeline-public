# ============================================================================
# run_all.R — run the whole pipeline end to end and print the PASS/FAIL summary.
# Usage (from inside sim_pipeline/):
#     Rscript run_all.R            # full default profile (~minutes)
#     Rscript run_all.R tiny       # fast smoke test
# ============================================================================
args <- commandArgs(trailingOnly = TRUE)
profile <- if (length(args) >= 1) args[1] else Sys.getenv("SIM_PROFILE", "default")
Sys.setenv(SIM_PROFILE = profile)
cat("Running pipeline with profile:", profile, "\n")

if (file.exists("results/validation_log.tsv")) file.remove("results/validation_log.tsv")
source("R/utils.R")

scripts <- c("00_simulate.R", "01_qc_cluster.R", "02_dge.R",
             "03_deconvolution.R", "04_spatial.R")
t0 <- Sys.time()
for (s in scripts) {
  cat("\n>>>>>>>>>>>>>>>>>>>> running", s, ">>>>>>>>>>>>>>>>>>>>\n")
  source(s)
}
cat(sprintf("\nTotal pipeline time: %.2f min\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))

cat("\n==================== VALIDATION SUMMARY ====================\n")
log <- read.delim("results/validation_log.tsv")
print(log[, c("stage", "status", "label")], row.names = FALSE)
cat(sprintf("\n%d/%d checks PASSED\n", sum(log$status == "PASS"), nrow(log)))

save_session_info("sessionInfo.txt")

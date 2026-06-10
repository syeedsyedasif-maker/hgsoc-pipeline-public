# ============================================================================
# R/utils.R — shared plumbing used by every script.
#   - load_config(): read config.yml, merge the chosen profile over `default`
#   - logging + the one-line PASS/FAIL printer that each stage emits
#   - small helpers (seeding, saving inspectable intermediates, sessionInfo)
# Nothing here is biology; it is the scaffolding that keeps the 5 scripts tidy.
# ============================================================================

suppressPackageStartupMessages({
  library(yaml)
})

# ---- recursive list merge: override wins, nested lists merged key-by-key -----
merge_lists <- function(base, override) {
  for (k in names(override)) {
    if (is.list(base[[k]]) && is.list(override[[k]])) {
      base[[k]] <- merge_lists(base[[k]], override[[k]])
    } else {
      base[[k]] <- override[[k]]
    }
  }
  base
}

# ---- load config.yml and merge profile (env var SIM_PROFILE, default "default")
load_config <- function(path = "config.yml", profile = Sys.getenv("SIM_PROFILE", "default")) {
  raw <- yaml::read_yaml(path)
  cfg <- raw[["default"]]
  if (profile != "default") {
    if (is.null(raw[[profile]])) stop("Unknown profile in config.yml: ", profile)
    cfg <- merge_lists(cfg, raw[[profile]])
  }
  cfg$.profile <- profile
  cfg
}

# ---- canonical cell-type / group order (single source of truth) -------------
# The 7 splatter groups, in the order their proportions appear in config.
group_levels <- function(cfg) names(cfg$group_prob)

# Map a group name to its higher-level cell type (collapses the two malignant).
celltype_of_group <- function(g) {
  ifelse(grepl("^malignant", g), "malignant", g)
}

# ---- threads: single-threaded by default for reproducibility -----------------
set_threads <- function(cfg) {
  n <- as.integer(cfg$threads %||% 1L)
  options(mc.cores = n)
  Sys.setenv(OMP_NUM_THREADS = n, OPENBLAS_NUM_THREADS = n)
  invisible(n)
}
`%||%` <- function(a, b) if (is.null(a)) b else a

# ---- per-script init: print header, seed, set threads -----------------------
init_script <- function(name, cfg) {
  cat(strrep("=", 74), "\n", sep = "")
  cat("  ", name, "   [profile: ", cfg$.profile, " | seed: ", cfg$seed,
      " | threads: ", cfg$threads, "]\n", sep = "")
  cat(strrep("=", 74), "\n", sep = "")
  set.seed(cfg$seed)
  set_threads(cfg)
  invisible(NULL)
}

log_msg <- function(...) cat("  -", ..., "\n")

# ---- the PASS/FAIL line each stage prints, also appended to results/ ---------
# Returns the logical so callers can stop()/branch if they wish.
pass_fail <- function(stage, label, condition, detail = "",
                      logfile = "results/validation_log.tsv") {
  status <- if (isTRUE(condition)) "PASS" else "FAIL"
  cat(sprintf("\n[%s] %s :: %s %s\n", status, stage, label,
              if (nzchar(detail)) paste0("(", detail, ")") else ""))
  if (!dir.exists(dirname(logfile))) dir.create(dirname(logfile), recursive = TRUE)
  line <- data.frame(time = as.character(Sys.time()), stage = stage,
                     label = label, status = status, detail = detail)
  write.table(line, logfile, sep = "\t", row.names = FALSE,
              col.names = !file.exists(logfile), append = file.exists(logfile),
              quote = FALSE)
  invisible(isTRUE(condition))
}

# ---- save a small intermediate the user can open and eyeball -----------------
save_peek <- function(obj, file, n = 8) {
  utils::write.csv(utils::head(as.data.frame(obj), n), file, row.names = TRUE)
  log_msg("peek ->", file)
}

# ---- write sessionInfo for reproducibility ----------------------------------
save_session_info <- function(file = "sessionInfo.txt") {
  writeLines(capture.output(sessionInfo()), file)
  log_msg("sessionInfo ->", file)
}

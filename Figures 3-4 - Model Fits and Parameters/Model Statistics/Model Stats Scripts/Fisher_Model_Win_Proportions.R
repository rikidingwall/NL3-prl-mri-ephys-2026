# R script: fisher_one_stage_with_posthoc.R
# Single-stage analysis. Set WT_col and G1_col (exact names), run the stage.
# If initial Fisher p < alpha, run per-model 2x2 Fisher post-hoc tests and BH-adjust p-values.

# ----------------- USER SETTINGS -----------------
input_csv <- "model_segment_wins_by_stage.csv"  # path to your CSV
WT_col    <- "PRL_Probe_3_BIC_WT"        # EXACT column name in CSV for WT for this stage
G1_col    <- "PRL_Probe_3_BIC_NL3" # EXACT column name in CSV for Genotype1 for this stage
B         <- 1e5                # Monte-Carlo replicates for fisher.test; e.g. 1e5 or 1e6
alpha     <- 0.05               # significance threshold for running post-hoc tests
verbose   <- TRUE
# -------------------------------------------------

# --- helpers ---
compute_chi2_stat <- function(tbl) {
  tbl <- as.matrix(tbl)
  n <- sum(tbl)
  if (n == 0) return(list(chi2 = NA_real_, expected = matrix(0, nrow=nrow(tbl), ncol=ncol(tbl))))
  rs <- rowSums(tbl); cs <- colSums(tbl)
  expected <- outer(rs, cs, FUN = function(r, c) r * c / n)
  valid <- expected > 0
  chi2 <- NA_real_
  if (any(valid)) chi2 <- sum(((tbl - expected)^2 / expected)[valid])
  list(chi2 = as.numeric(chi2), expected = expected)
}
cramers_v <- function(chi2, n, r, c) {
  if (is.na(chi2) || n <= 0) return(NA_real_)
  denom <- n * (min(r, c) - 1)
  if (denom <= 0) return(NA_real_)
  sqrt(chi2 / denom)
}

# --- read data and checks ---
df <- read.csv(input_csv, stringsAsFactors = FALSE, check.names = FALSE)
if (! "model" %in% names(df)) stop("CSV must contain a 'model' column with model names.")
if (! WT_col %in% names(df)) stop(sprintf("WT_col '%s' not found in CSV columns.", WT_col))
if (! G1_col %in% names(df)) stop(sprintf("G1_col '%s' not found in CSV columns.", G1_col))

# --- build contingency table (rows = models, cols = WT / Genotype1) ---
tab_df <- data.frame(
  WT = as.integer(df[[WT_col]]),
  Genotype1 = as.integer(df[[G1_col]]),
  stringsAsFactors = FALSE
)
rownames(tab_df) <- df$model

if (any(is.na(tab_df))) stop("NA values present in the selected columns.")
if (any(tab_df < 0)) stop("Negative counts present in the selected columns.")

tbl <- as.matrix(tab_df)
total_n <- sum(tbl)

cat("\n========================\n")
cat(sprintf("Testing stage: WT_col='%s'  vs  G1_col='%s'\n", WT_col, G1_col))
cat("========================\n\n")

cat("Contingency table (rows = models, columns = WT / Genotype1):\n")
print(tbl)
cat("\nRow sums (models):\n"); print(rowSums(tbl))
cat("Column sums (genotypes):\n"); print(colSums(tbl))
cat(sprintf("\nTotal segments (n): %d\n", total_n))

if (total_n == 0) {
  cat("Total count is 0 — nothing to test.\n")
  quit(save="no")
}

# --- initial Fisher test (simulated) ---
cat(sprintf("\nRunning fisher.test(simulate.p.value = TRUE, B = %d) for overall 2xK table...\n", as.integer(B)))
fisher_overall <- tryCatch(
  fisher.test(tbl, simulate.p.value = TRUE, B = as.integer(B)),
  error = function(e) e
)
if (inherits(fisher_overall, "error")) {
  stop(sprintf("Overall fisher.test failed: %s", fisher_overall$message))
}

p_overall <- fisher_overall$p.value
method_overall <- fisher_overall$method

chi_info <- compute_chi2_stat(tbl)
chi2_stat <- chi_info$chi2
expected_mat <- chi_info$expected
small_expected_cells <- if (is.matrix(expected_mat)) sum(expected_mat < 5) else NA_integer_
v <- cramers_v(chi2_stat, total_n, nrow(tbl), ncol(tbl))

cat("\nStage summary:\n")
cat(sprintf("  Fisher p-value (simulated): %s\n", format(p_overall, digits=6)))
cat(sprintf("  Method used: %s\n", method_overall))
cat(sprintf("  Chi-square (effect size only): %s\n", ifelse(is.na(chi2_stat), "NA", format(chi2_stat, digits=6))))
cat(sprintf("  Cramér's V: %s\n", ifelse(is.na(v), "NA", format(v, digits=6))))
cat(sprintf("  Expected cells < 5: %d\n", small_expected_cells))
cat("----------------------------------------------------\n")

# --- conditional post-hoc per-model 2x2 Fisher tests ---
if (!is.na(p_overall) && p_overall < alpha) {
  cat("\nOverall test significant (p < alpha). Running per-model 2x2 Fisher post-hoc tests...\n")
  models <- rownames(tbl)
  posthoc_list <- vector("list", length(models))
  for (i in seq_along(models)) {
    m <- models[i]
    a <- tbl[m, "WT"]
    b <- tbl[m, "Genotype1"]
    others_wt <- sum(tbl[-i, "WT"])
    others_g1 <- sum(tbl[-i, "Genotype1"])
    mat <- matrix(c(a, b, others_wt, others_g1), nrow = 2, byrow = TRUE)
    # rows: [model; others], cols: [WT, Genotype1]
    # fisher.test usually returns estimate = odds ratio for 2x2, unless degenerate.
    ft <- tryCatch(
      fisher.test(mat, simulate.p.value = FALSE),
      error = function(e) e
    )
    if (inherits(ft, "error")) {
      # fallback to simulated p-value if exact fails
      ft <- tryCatch(fisher.test(mat, simulate.p.value = TRUE, B = as.integer(B)), error = function(e) e)
      if (inherits(ft, "error")) {
        # ultimate fallback: record NA and continue
        posthoc_list[[i]] <- list(model = m, a_WT = a, b_G1 = b,
                                  OR = NA_real_, conf_low = NA_real_, conf_high = NA_real_,
                                  p_value = NA_real_, method = "fisher_failed")
        next
      } else {
        method_used <- "fisher_simulated"
      }
    } else {
      method_used <- "fisher_exact"
    }
    
    OR <- if (!is.null(ft$estimate)) as.numeric(ft$estimate) else NA_real_
    CIs <- if (!is.null(ft$conf.int)) ft$conf.int else c(NA_real_, NA_real_)
    pval <- ft$p.value
    posthoc_list[[i]] <- list(model = m, a_WT = a, b_G1 = b,
                              OR = OR, conf_low = CIs[1], conf_high = CIs[2],
                              p_value = pval, method = method_used)
  }
  
  # combine into data.frame and BH adjust
  post_df <- do.call(rbind, lapply(posthoc_list, as.data.frame))
  post_df$p_value <- as.numeric(as.character(post_df$p_value))
  # Adjust only on non-NA p-values
  valid_idx <- which(!is.na(post_df$p_value))
  post_df$p_adj_BH <- NA_real_
  if (length(valid_idx) > 0) {
    post_df$p_adj_BH[valid_idx] <- p.adjust(post_df$p_value[valid_idx], method = "BH")
  }
  
  # Interpret OR direction:
  # matrix layout was: rows = (model, others), cols = (WT, Genotype1)
  # OR = (a * others_g1) / (b * others_wt)
  # Therefore: OR > 1 => model is relatively enriched in WT (vs Genotype1).
  post_df$interpretation <- ifelse(is.na(post_df$OR), NA_character_,
                                   ifelse(post_df$OR > 1, "enriched_in_WT", "enriched_in_Genotype1"))
  
  # order by adjusted p then raw p
  post_df <- post_df[order(post_df$p_adj_BH, post_df$p_value), ]
  
  cat("\nPer-model post-hoc 2x2 Fisher results (rows = model, cols = WT / Genotype1):\n")
  print(post_df, digits = 6, right = TRUE, row.names = FALSE)
  cat("\nNotes: OR > 1 means the model is relatively more likely in WT (vs Genotype1).\n")
  cat("P-values: raw and BH-adjusted (across the models tested).\n")
} else {
  cat("\nOverall test NOT significant (p >= alpha). Skipping per-model post-hoc Fisher tests.\n")
}
cat("\nDone.\n")

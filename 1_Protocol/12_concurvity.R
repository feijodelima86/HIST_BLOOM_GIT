# ============================================================================
# 12_concurvity.R
# UCFR Filamentous Algae Project
# Concurvity diagnostics for M1 — supplementary table
#
# Input:   3_models/bloom_model_M1.rds
# Outputs: 2_incremental/concurvity_worst.csv   (pairwise worst-case matrix)
#          2_incremental/concurvity_summary.csv  (per-term observed summary)
#          console: scorecard
#
# Concurvity is the GAM analogue of multicollinearity: the degree to which
# one smooth term can be approximated by a combination of the other smooths.
# Three measures returned by mgcv::concurvity():
#
#   worst    — upper bound; assumes worst-case linear combination of others.
#              Most conservative; this is what referees will ask about.
#   observed — concurvity with the actual fitted smooths (not worst-case).
#              More realistic but model-specific.
#   estimate — based on the estimated smooths only (ignores uncertainty).
#
# We report:
#   (1) Pairwise worst-case matrix — identifies which predictor pairs
#       share the most information. Flags values > 0.8 (concern) and
#       > 0.9 (serious concern).
#   (2) Per-term summary — each term's worst and observed concurvity
#       vs the rest of the model combined. The single number for the
#       supplement.
#
# Ecological context for known overlaps:
#   anomaly × logQ     — both hydrology-derived; anomaly = f(peak Q,
#                        baseflow Q), logQ = summer observed Q. Structural
#                        overlap possible but they represent different
#                        temporal windows (spring vs summer).
#   Days_Since_Freshet × logQ — both increase as season progresses and
#                        flow recedes. Most likely high-concurvity pair.
#   lag_y × logTP      — prior bloom depletes P; high lag_y may correlate
#                        with low current TP. Moderate overlap expected.
# ============================================================================

library(mgcv)

# ============================================================================
# 1. LOAD MODEL
# ============================================================================
m_M1     <- readRDS("3_models/bloom_model_M1.rds")
RESPONSE <- as.character(formula(m_M1)[[2]])

# ============================================================================
# 2. COMPUTE CONCURVITY
# ============================================================================
# full=TRUE: each term vs rest of model combined.
# full=FALSE: pairwise term x term matrices.
#
# mgcv returns concurvity(full=TRUE) as a matrix with rows = measures
# (worst, observed, estimate) and cols = terms — NOT a named list.
# Extract by row name.

cc_full <- concurvity(m_M1, full = TRUE)
cc_pair <- concurvity(m_M1, full = FALSE)

# full=TRUE: rows are measures, columns are terms
worst_vec    <- cc_full["worst",    ]
observed_vec <- cc_full["observed", ]

# full=FALSE: list with named matrices (worst, observed, estimate)
worst_mat    <- cc_pair$worst

# ============================================================================
# 3. FORMAT PAIRWISE WORST-CASE MATRIX
# ============================================================================
# Diagonal = 1 by definition; set to NA for readability
diag(worst_mat) <- NA

clean_names <- function(x) {
  x <- gsub('^s\\(',      '',        x)
  x <- gsub('\\)$',       '',        x)
  x <- gsub(', bs="re"',  ' [RE]',   x)
  x
}
rownames(worst_mat) <- clean_names(rownames(worst_mat))
colnames(worst_mat) <- clean_names(colnames(worst_mat))

worst_df <- as.data.frame(round(worst_mat, 4))

# ============================================================================
# 4. FORMAT PER-TERM SUMMARY
# ============================================================================
term_names_clean <- clean_names(names(worst_vec))

summary_df <- data.frame(
  Term          = term_names_clean,
  Worst_vs_rest = round(worst_vec,    4),
  Obs_vs_rest   = round(observed_vec, 4),
  Flag          = ifelse(worst_vec > 0.9, "!!",
                         ifelse(worst_vec > 0.8, "!",  "")),
  stringsAsFactors = FALSE
)
row.names(summary_df) <- NULL

# ============================================================================
# 5. SAVE OUTPUTS
# ============================================================================
if (!dir.exists("2_incremental")) dir.create("2_incremental", recursive = TRUE)

write.csv(worst_df,   "2_incremental/concurvity_worst.csv")
write.csv(summary_df, "2_incremental/concurvity_summary.csv",
          row.names = FALSE)

# ============================================================================
# 6. SCORECARD
# ============================================================================
cat("\n")
cat("============================================================\n")
cat("CONCURVITY DIAGNOSTICS — M1\n")
cat("Response:", RESPONSE, "\n")
cat("============================================================\n\n")

cat("--- Per-term: worst-case and observed concurvity vs rest of model ---\n")
cat("(Values > 0.8 flagged '!'; > 0.9 flagged '!!')\n\n")
print(summary_df, row.names = FALSE)
cat("\n")

cat("--- Pairwise worst-case concurvity matrix ---\n")
cat("(Row = term being approximated; col = approximating term)\n")
cat("(Diagonal = NA; values > 0.8 warrant discussion)\n\n")
print(worst_df)
cat("\n")

# Identify flagged pairs explicitly
high_pairs <- which(worst_mat > 0.8, arr.ind = TRUE)
if (nrow(high_pairs) > 0) {
  cat("--- Flagged pairs (worst-case > 0.8) ---\n")
  pair_df <- data.frame(
    Term_A = rownames(worst_mat)[high_pairs[, 1]],
    Term_B = colnames(worst_mat)[high_pairs[, 2]],
    Worst  = round(worst_mat[high_pairs], 4),
    stringsAsFactors = FALSE
  )
  pair_df <- pair_df[pair_df$Term_A < pair_df$Term_B, ]
  if (nrow(pair_df) > 0) {
    print(pair_df, row.names = FALSE)
  } else {
    cat("None after removing symmetric duplicates.\n")
  }
} else {
  cat("--- No pairwise worst-case values exceed 0.8 ---\n")
  cat("Concurvity is not a concern for M1.\n")
}
cat("\n")

cat("--- Interpretation guide ---\n")
cat("< 0.8  : acceptable\n")
cat("0.8-0.9: moderate concern — interpret affected terms cautiously\n")
cat("> 0.9  : high concern — smooth estimates may be unreliable\n\n")

cat("--- Output files ---\n")
sc_files <- data.frame(
  File = c("2_incremental/concurvity_worst.csv",
           "2_incremental/concurvity_summary.csv"),
  Contents = c("Pairwise worst-case concurvity matrix.",
               "Per-term worst and observed concurvity vs rest of model."),
  stringsAsFactors = FALSE
)
print(sc_files, row.names = FALSE)
cat("\n")
cat("============================================================\n")
cat("Done. Paste summary_df into supplementary Table S-X.\n")
cat("============================================================\n")
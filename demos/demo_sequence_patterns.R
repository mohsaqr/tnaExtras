library(tnaExtras)

# ==============================================================================
# DEMO: Unified Sequence Pattern Analysis
# ==============================================================================
# This demo shows how to use the unified explore_patterns() function
# for n-gram extraction, abstract patterns, and full sequence analysis.
# ==============================================================================

# Load sample data
load("data/seqdata.rda")
seq_data <- seqdata[, -1]  # Remove group column if present

cat("======================================================\n")
cat("SEQUENCE PATTERN ANALYSIS DEMO\n")
cat("======================================================\n\n")

# ------------------------------------------------------------------------------
# 1. N-GRAM ANALYSIS
# ------------------------------------------------------------------------------
cat("=== 1. N-GRAM ANALYSIS ===\n\n")

# Basic n-gram extraction
results <- explore_patterns(
  seq_data,
  type = "ngrams",
  min_length = 2,
  max_length = 4,
  min_support = 0.05,
  correction = "fdr"
)

print(results)

# Get significant patterns
sig <- significant_patterns(results)
cat("\nSignificant patterns:\n")
print(head(sig, 10))

# ------------------------------------------------------------------------------
# 2. FILTERING BY STATE
# ------------------------------------------------------------------------------
cat("\n=== 2. FILTERING BY STATE ===\n\n")

# Patterns starting with "plan"
plan_patterns <- explore_patterns(
  seq_data,
  type = "ngrams",
  max_length = 3,
  start_state = "plan",
  verbose = FALSE
)
cat("Patterns starting with 'plan':", nrow(plan_patterns$patterns), "\n")
print(head(plan_patterns$patterns[, c("pattern", "count", "support", "significant")], 10))

# Patterns ending with "consensus"
consensus_patterns <- explore_patterns(
  seq_data,
  type = "ngrams",
  max_length = 3,
  end_state = "consensus",
  verbose = FALSE
)
cat("\nPatterns ending with 'consensus':", nrow(consensus_patterns$patterns), "\n")
print(head(consensus_patterns$patterns[, c("pattern", "count", "support")], 10))

# Patterns containing "emotion"
emotion_patterns <- explore_patterns(
  seq_data,
  type = "ngrams",
  max_length = 3,
  contains_state = "emotion",
  verbose = FALSE
)
cat("\nPatterns containing 'emotion':", nrow(emotion_patterns$patterns), "\n")
print(head(emotion_patterns$patterns[, c("pattern", "count", "support")], 10))

# ------------------------------------------------------------------------------
# 3. ABSTRACT PATTERN DETECTION
# ------------------------------------------------------------------------------
cat("\n=== 3. ABSTRACT PATTERN DETECTION ===\n\n")

abstract <- explore_patterns(
  seq_data,
  type = "abstract",
  min_gap = 1,
  max_gap = 3,
  min_support = 0.05
)

cat("Abstract patterns found:", nrow(abstract$patterns), "\n\n")
print(head(abstract$patterns, 15))

# ------------------------------------------------------------------------------
# 4. FULL SEQUENCE ANALYSIS
# ------------------------------------------------------------------------------
cat("\n=== 4. FULL SEQUENCE ANALYSIS ===\n\n")

full_seqs <- explore_patterns(
  seq_data,
  type = "full",
  min_support = 0.01,
  verbose = FALSE
)

cat("Unique sequences:", nrow(full_seqs$patterns), "\n\n")
cat("Top 10 most frequent sequences:\n")
print(head(full_seqs$patterns[, c("pattern", "count", "support")], 10))

# ------------------------------------------------------------------------------
# 5. VISUALIZATION
# ------------------------------------------------------------------------------
cat("\n=== 5. VISUALIZATION ===\n\n")

# Plot top patterns
# Note: Plotting requires a graphics device. In non-interactive mode, this might not show.
if (interactive()) {
  plot(results, type = "patterns", top_n = 15)
}

# ------------------------------------------------------------------------------
# 6. HELPER FUNCTIONS
# ------------------------------------------------------------------------------
cat("\n=== 6. HELPER FUNCTIONS ===\n\n")

# Filter patterns programmatically
filtered <- filter_patterns(
  results,
  min_support = 0.1,
  pattern_length = 2,
  significant_only = TRUE
)
cat("Filtered patterns (support > 0.1, length 2, significant):\n")
print(filtered)

# Get top sequences
top <- top_sequences(results, top_n = 5, significant_only = TRUE)
cat("\nTop 5 significant patterns:\n")
print(top)

cat("\n======================================================\n")
cat("DEMO COMPLETE\n")
cat("======================================================\n")

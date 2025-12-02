# ==============================================================================
# DEMO: Sequence Motif Analysis
# ==============================================================================
#
# This demo shows how to use the sequence motif analysis functions:
# - detect_abstract_patterns(): Find returns, oscillations, progressions
# - find_gapped_patterns(): Find patterns with wildcards (A-*-B)
# - find_meta_paths(): Discover meta-paths across node types
#
# Using regulation data from the tna package
# ==============================================================================

# Load packages
library(tnaExtras)

# Load regulation data from tna package
# If tna is not installed, use the included regulation_grouped dataset
if (requireNamespace("tna", quietly = TRUE)) {
  data <- tna::group_regulation
  cat("Using tna::group_regulation data\n")
} else {
  data(regulation_grouped)
  data <- regulation_grouped
  cat("Using tnaExtras::regulation_grouped data\n")
}

# Remove group column for analysis
seq_data <- data[, -1]

cat("\n=== Data Overview ===\n")
cat("Sequences:", nrow(seq_data), "\n")
cat("Time points:", ncol(seq_data), "\n")
cat("Sample sequence:", paste(head(as.character(seq_data[1, ]), 10), collapse = " -> "), "...\n")

# ==============================================================================
# 1. ABSTRACT PATTERN DETECTION
# ==============================================================================

cat("\n\n")
cat("=" ,rep("=", 60), "\n", sep = "")
cat("1. ABSTRACT PATTERN DETECTION\n")
cat("=" ,rep("=", 60), "\n", sep = "")

# Detect all abstract patterns
abstract <- detect_abstract_patterns(
  data = seq_data,
  patterns = "all",        # Detect all pattern types
  min_gap = 1,
  max_gap = 3,
  min_support = 0.05,
  verbose = TRUE
)

print(abstract)

# Access specific pattern types
cat("\n--- Return Patterns (A->*->A) ---\n")
print(head(abstract$returns, 10))

cat("\n--- Repetition Patterns (A->A) ---\n")
print(head(abstract$repetitions, 10))

cat("\n--- Oscillation Patterns (A->B->A->B) ---\n")
if (nrow(abstract$oscillations) > 0) {
  print(head(abstract$oscillations, 10))
} else {
  cat("No oscillation patterns found with current thresholds\n")
}

# Detect only returns with specific gap range
cat("\n--- Returns with Gap 1-2 ---\n")
returns_only <- detect_abstract_patterns(
  seq_data,
  patterns = "returns",
  min_gap = 1,
  max_gap = 2,
  min_support = 0.1,
  verbose = FALSE
)
print(head(returns_only$returns, 10))

# ==============================================================================
# 2. GAP-CONSTRAINED PATTERN FINDING
# ==============================================================================

cat("\n\n")
cat("=" ,rep("=", 60), "\n", sep = "")
cat("2. GAP-CONSTRAINED PATTERN FINDING\n")
cat("=" ,rep("=", 60), "\n", sep = "")

# Auto-discover gapped patterns
cat("\n--- Auto-discovered Gapped Patterns ---\n")
gapped <- find_gapped_patterns(
  data = seq_data,
  pattern = NULL,          # Auto-discover
  min_gap = 1,
  max_gap = 2,
  min_support = 0.05,
  verbose = TRUE
)

print(gapped)

# View top patterns
cat("\n--- Top 15 Gapped Patterns ---\n")
print(head(gapped$patterns, 15))

# Search for specific patterns
cat("\n--- Specific Pattern: plan-*-consensus ---\n")
specific <- find_gapped_patterns(
  seq_data,
  pattern = "plan-*-consensus",
  verbose = FALSE
)
print(specific$patterns)
cat("\nInstances found:\n")
print(head(specific$instances, 10))

# Multi-gap pattern
cat("\n--- Multi-gap Pattern: plan-**-plan (plan returns) ---\n")
multi_gap <- find_gapped_patterns(
  seq_data,
  pattern = "plan-**-plan",
  verbose = FALSE
)
print(multi_gap$patterns)

# ==============================================================================
# 3. META-PATH DISCOVERY
# ==============================================================================

cat("\n\n")
cat("=" ,rep("=", 60), "\n", sep = "")
cat("3. META-PATH DISCOVERY\n")
cat("=" ,rep("=", 60), "\n", sep = "")

# Define node types for regulation states
node_types <- list(
  cognitive = c("plan", "monitor", "adapt"),
  social = c("discuss", "consensus", "coregulate", "synthesis"),
  emotional = c("emotion", "cohesion")
)

cat("\nNode Types Defined:\n")
for (type in names(node_types)) {
  cat(sprintf("  %s: %s\n", type, paste(node_types[[type]], collapse = ", ")))
}

# Auto-discover all meta-paths (hybrid discovery)
cat("\n--- Auto-discovering Meta-Paths ---\n")
meta <- find_meta_paths(
  data = seq_data,
  node_types = node_types,
  schema = NULL,           # Auto-discover (hybrid)
  min_length = 2,
  max_length = 4,
  min_support = 0.05,
  verbose = TRUE
)

print(meta)

# Detailed summary
cat("\n--- Detailed Summary ---\n")
summary(meta)

# View meta-paths by length
cat("\n--- Meta-Paths by Length ---\n")
for (len in 2:4) {
  len_paths <- meta$meta_paths[meta$meta_paths$length == len, ]
  if (nrow(len_paths) > 0) {
    cat(sprintf("\nLength %d (%d patterns):\n", len, nrow(len_paths)))
    print(head(len_paths[order(len_paths$count, decreasing = TRUE), 
                         c("schema", "count", "support", "significant")], 5))
  }
}

# View type transitions
cat("\n--- Type-to-Type Transitions ---\n")
print(meta$type_transitions)

# Search for specific schema (using vector format - recommended)
cat("\n--- Specific Schema: cognitive->social->cognitive ---\n")
meta_specific <- find_meta_paths(
  seq_data,
  node_types = node_types,
  schema = c("cognitive", "social", "cognitive"),  # Vector format
  verbose = FALSE
)
print(meta_specific$meta_paths)
cat("\nTop instances:\n")
print(head(meta_specific$instances, 10))

# Schema with wildcards (using vector format)
cat("\n--- Schema with Wildcard: cognitive->**->cognitive (cognitive returns) ---\n")
meta_return <- find_meta_paths(
  seq_data,
  node_types = node_types,
  schema = c("cognitive", "**", "cognitive"),  # Vector with wildcard
  verbose = FALSE
)
print(meta_return$meta_paths)
cat("\nSample instances:\n")
print(head(meta_return$instances, 5))

# ==============================================================================
# 4. COMBINING ANALYSES
# ==============================================================================

cat("\n\n")
cat("=" ,rep("=", 60), "\n", sep = "")
cat("4. COMBINING ANALYSES\n")
cat("=" ,rep("=", 60), "\n", sep = "")

# Find significant patterns across all analyses
cat("\n--- Significant Abstract Patterns ---\n")
if (!is.null(abstract$returns) && "significant" %in% names(abstract$returns)) {
  sig_returns <- abstract$returns[abstract$returns$significant, ]
  cat(sprintf("Significant returns: %d\n", nrow(sig_returns)))
}

cat("\n--- Significant Gapped Patterns ---\n")
if ("significant" %in% names(gapped$patterns)) {
  sig_gapped <- gapped$patterns[gapped$patterns$significant, ]
  cat(sprintf("Significant gapped patterns: %d\n", nrow(sig_gapped)))
}

cat("\n--- Significant Meta-Paths ---\n")
if ("significant" %in% names(meta$meta_paths)) {
  sig_meta <- meta$meta_paths[meta$meta_paths$significant, ]
  cat(sprintf("Significant meta-paths: %d\n", nrow(sig_meta)))
  cat("\nTop 5 significant meta-paths:\n")
  print(head(sig_meta[order(sig_meta$lift, decreasing = TRUE), 
                       c("schema", "count", "support", "lift", "p_adjusted")], 5))
}

# ==============================================================================
# SUMMARY
# ==============================================================================

cat("\n\n")
cat("=" ,rep("=", 60), "\n", sep = "")
cat("ANALYSIS SUMMARY\n")
cat("=" ,rep("=", 60), "\n", sep = "")

cat("\nAbstract Patterns:\n")
cat(sprintf("  Returns: %d (significant: %d)\n", 
           abstract$summary$n_returns,
           if (!is.null(abstract$returns) && "significant" %in% names(abstract$returns))
             sum(abstract$returns$significant) else 0))
cat(sprintf("  Repetitions: %d\n", abstract$summary$n_repetitions))
cat(sprintf("  Oscillations: %d\n", abstract$summary$n_oscillations))
cat(sprintf("  Progressions: %d\n", abstract$summary$n_progressions))

cat("\nGapped Patterns:\n")
cat(sprintf("  Total: %d (significant: %d)\n",
           gapped$summary$n_patterns,
           if (!is.na(gapped$summary$n_significant)) gapped$summary$n_significant else 0))

cat("\nMeta-Paths:\n")
cat(sprintf("  Total: %d (significant: %d)\n",
           meta$summary$n_meta_paths,
           if (!is.na(meta$summary$n_significant)) meta$summary$n_significant else 0))

cat("\n=== Demo Complete ===\n")


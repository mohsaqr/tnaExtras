# ==============================================================================
# DEMO: Sequence Pattern Exploration
# ==============================================================================
# 
# This demo shows how to use explore_sequence_patterns() to discover frequent
# patterns, compute support frequencies, and test statistical significance
# in sequence data.
#
# ==============================================================================

# Load the package
library(tnaExtras)

# Load sample data
data(seqdata)

# View the data structure
cat("=== Sample Data Overview ===\n")
cat("Rows (sequences):", nrow(seqdata), "\n")
cat("Columns (time points + group):", ncol(seqdata), "\n")
cat("\nColumn names:", paste(names(seqdata), collapse = ", "), "\n")
cat("\nFirst 5 rows:\n")
print(head(seqdata, 5))

# ==============================================================================
# BASIC USAGE
# ==============================================================================

cat("\n\n=== Basic Pattern Exploration ===\n")

# Remove group column for single-group analysis
sequence_data <- seqdata[, -1]  # Remove Group column

# Run pattern exploration
results <- explore_sequence_patterns(
  data = sequence_data,
  min_length = 1,
  max_length = 4,
  min_support = 0.05,  # Patterns appearing in at least 5% of sequences
  min_count = 3,       # Patterns appearing at least 3 times
  verbose = TRUE
)

# View results
print(results)

# ==============================================================================
# DETAILED SUMMARY
# ==============================================================================

cat("\n\n=== Detailed Summary ===\n")
summary(results)

# ==============================================================================
# WORKING WITH PATTERNS
# ==============================================================================

cat("\n\n=== Working with Patterns ===\n")

# View all patterns
cat("All patterns (first 20):\n")
print(head(results$patterns, 20))

# Get significant patterns only
cat("\n\nSignificant patterns:\n")
sig_patterns <- significant_patterns(results)
print(sig_patterns)

# Filter patterns by criteria
cat("\n\nFiltering patterns (length=2, min_support=0.1):\n")
filtered <- filter_patterns(
  results, 
  min_support = 0.1, 
  pattern_length = 2
)
print(filtered)

# ==============================================================================
# STATE AND TRANSITION ANALYSIS
# ==============================================================================

cat("\n\n=== State Frequencies ===\n")
print(results$state_frequencies)

cat("\n\n=== Top Transitions ===\n")
print(head(results$transitions, 15))

# Get transition matrix
cat("\n\n=== Transition Probability Matrix ===\n")
trans_matrix <- get_transition_matrix(results, type = "probability")
print(round(trans_matrix, 3))

# ==============================================================================
# FULL SEQUENCE ANALYSIS
# ==============================================================================

cat("\n\n=== Most Common Complete Sequences ===\n")
print(head(results$full_sequences, 10))

# ==============================================================================
# VISUALIZATION
# ==============================================================================

cat("\n\n=== Visualizations ===\n")

# Plot top patterns
par(mfrow = c(2, 2))

# Top patterns by frequency
plot(results, type = "patterns", top_n = 15)

# Top states
plot(results, type = "states", top_n = 10)

# Top transitions
plot(results, type = "transitions", top_n = 10)

par(mfrow = c(1, 1))

# ==============================================================================
# ANALYSIS WITH DIFFERENT PARAMETERS
# ==============================================================================

cat("\n\n=== Analysis with Stricter Thresholds ===\n")

# Stricter analysis
results_strict <- explore_sequence_patterns(
  data = sequence_data,
  min_length = 2,
  max_length = 3,
  min_support = 0.10,    # 10% support
  min_count = 5,         # At least 5 occurrences
  correction = "bonferroni",  # Stricter correction
  alpha = 0.01,          # Stricter significance
  verbose = TRUE
)

print(results_strict)

# ==============================================================================
# GROUP-SPECIFIC ANALYSIS
# ==============================================================================

cat("\n\n=== Group-Specific Pattern Analysis ===\n")

# Analyze Group A only
group_A_data <- seqdata[seqdata$Group == "A", -1]
cat("Analyzing Group A (", nrow(group_A_data), " sequences):\n")

results_A <- explore_sequence_patterns(
  data = group_A_data,
  min_length = 2,
  max_length = 3,
  min_support = 0.05,
  verbose = FALSE
)

cat("Group A - Significant patterns:", results_A$summary$n_significant_patterns, "\n")
cat("Most common pattern:", results_A$summary$most_common_pattern, "\n")

# Analyze Group B only
group_B_data <- seqdata[seqdata$Group == "B", -1]
cat("\nAnalyzing Group B (", nrow(group_B_data), " sequences):\n")

results_B <- explore_sequence_patterns(
  data = group_B_data,
  min_length = 2,
  max_length = 3,
  min_support = 0.05,
  verbose = FALSE
)

cat("Group B - Significant patterns:", results_B$summary$n_significant_patterns, "\n")
cat("Most common pattern:", results_B$summary$most_common_pattern, "\n")

# ==============================================================================
# ACCESSING RAW DATA
# ==============================================================================

cat("\n\n=== Accessing Analysis Components ===\n")

# Summary statistics
cat("Summary statistics:\n")
cat("  - Mean sequence length:", round(results$summary$mean_sequence_length, 2), "\n")
cat("  - State entropy:", round(results$summary$state_entropy, 3), "\n")
cat("  - Number of unique patterns:", results$summary$n_unique_patterns, "\n")
cat("  - Number of significant patterns:", results$summary$n_significant_patterns, "\n")

# Get parameters used
cat("\nAnalysis parameters:\n")
print(results$parameters[c("min_length", "max_length", "min_support", "correction")])

cat("\n=== Demo Complete ===\n")


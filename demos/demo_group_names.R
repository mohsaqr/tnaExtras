# =============================================================================
# DEMONSTRATION: ACTUAL GROUP NAMES IN ANALYSIS
# =============================================================================
# 
# This demo showcases how tnaExtras v0.3.0+ displays actual group names
# throughout all analysis functions instead of generic "A", "B", "C" labels.
#
# =============================================================================

library(tnaExtras)

cat("=== ACTUAL GROUP NAMES IN tnaExtras v0.3.0+ ===\n")
cat("Demonstrating how actual group names are preserved in all analysis functions\n\n")

# =============================================================================
# SAMPLE DATA WITH MEANINGFUL GROUP NAMES
# =============================================================================

educational_data <- data.frame(
  T1 = c('planning', 'discussing', 'monitoring', 'planning', 'discussing', 'evaluating'),
  T2 = c('consensus', 'emotion', 'planning', 'consensus', 'emotion', 'planning'),
  T3 = c('discussing', 'planning', 'consensus', 'discussing', 'planning', 'monitoring'),
  Group = c('High_Achievers', 'Average_Students', 'Struggling_Learners', 
            'High_Achievers', 'Average_Students', 'Struggling_Learners')
)

cat("Sample data with meaningful group names:\n")
print(educational_data)
cat("\n")

# =============================================================================
# 1. MULTI-GROUP PATTERN ANALYSIS
# =============================================================================

cat("1. MULTI-GROUP PATTERN ANALYSIS\n")
cat("================================\n")

patterns_multi <- analyze_patterns_multi(educational_data, 
                                        group_col = 'Group',
                                        min_frequency = 1)

cat("\nNotice how support columns use actual group names:\n")
support_columns <- grep("^support_", names(patterns_multi$support), value = TRUE)
cat("Support columns:", paste(support_columns, collapse = ", "), "\n")

cat("\nFirst few patterns with actual group names in columns:\n")
print(head(patterns_multi$support[c("pattern", support_columns)], 3))

# =============================================================================
# 2. TWO-GROUP PATTERN ANALYSIS
# =============================================================================

cat("\n\n2. TWO-GROUP PATTERN ANALYSIS\n")
cat("==============================\n")

# Use two groups for detailed analysis
two_group_data <- educational_data[educational_data$Group %in% c('High_Achievers', 'Struggling_Learners'), ]

patterns_two <- analyze_patterns(two_group_data, 
                               group_col = 'Group',
                               min_frequency = 1)

cat("\nTwo-group analysis with actual names in all measures:\n")
cat("Support columns:", paste(grep("^support_", names(patterns_two$support), value = TRUE), collapse = ", "), "\n")
cat("Lift columns:", paste(grep("^lift_", names(patterns_two$lift), value = TRUE), collapse = ", "), "\n")
cat("Confidence columns:", paste(grep("^confidence_", names(patterns_two$confidence), value = TRUE), collapse = ", "), "\n")

cat("\nSample support data:\n")
print(head(patterns_two$support[c("pattern", "support_High_Achievers", "support_Struggling_Learners")], 2))

# =============================================================================
# 3. SEQUENCE COMPARISON
# =============================================================================

cat("\n\n3. SEQUENCE COMPARISON\n")
cat("======================\n")

sequence_comp <- compare_sequences(two_group_data, 
                                 'Group', 
                                 detailed = FALSE, 
                                 statistical = FALSE)

cat("\nSequence comparison with actual group names in frequency columns:\n")
freq_columns <- grep("^freq_", names(sequence_comp$summary), value = TRUE)
cat("Frequency columns:", paste(freq_columns, collapse = ", "), "\n")

cat("\nSample frequency data:\n")
print(head(sequence_comp$summary[c("pattern", freq_columns)], 3))

# =============================================================================
# 4. SEQUENCE INDICES
# =============================================================================

cat("\n\n4. SEQUENCE INDICES\n")
cat("===================\n")

indices_result <- compute_sequence_indices(educational_data,
                                         group_col = 'Group',
                                         favorable_states = c('consensus', 'planning', 'discussing'),
                                         return_group_summary = TRUE)

cat("\nActual group names preserved in individual indices:\n")
print(unique(indices_result$individual_indices$Group))

cat("\nActual group names in group summaries:\n")
print(indices_result$group_summaries$Group)

# =============================================================================
# SUMMARY
# =============================================================================

cat("\n\n=== SUMMARY ===\n")
cat("✓ Multi-group pattern analysis: Uses actual group names in column headers\n")
cat("✓ Two-group pattern analysis: Uses actual group names across all measures\n")
cat("✓ Sequence comparison: Uses actual group names in frequency and proportion columns\n")
cat("✓ Sequence indices: Preserves actual group names in summaries and outputs\n")
cat("✓ Visualizations: Display actual group names on axes and legends\n")
cat("\nNo more generic 'A', 'B', 'C' labels - your meaningful group names are preserved!\n") 
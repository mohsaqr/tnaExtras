# =============================================================================
# DEMO: Multi-Group Analysis with tnaExtras v0.2.0
# =============================================================================

# Load the functions (in a real package, you would use library(tnaExtras))
source('R/pattern_analysis.R')
source('R/sequence_comparison.R') 
source('R/sequence_indices.R')

# =============================================================================
# SAMPLE DATA WITH 4 GROUPS
# =============================================================================

# Create example data with multiple groups representing different learning approaches
demo_data <- data.frame(
  T1 = c("plan", "discuss", "monitor", "evaluate", "plan", "discuss", "monitor", "evaluate",
         "plan", "discuss", "monitor", "evaluate", "plan", "discuss", "monitor", "evaluate"),
  T2 = c("discuss", "emotion", "plan", "reflect", "consensus", "emotion", "plan", "reflect",
         "emotion", "plan", "consensus", "reflect", "discuss", "plan", "emotion", "consensus"),
  T3 = c("consensus", "plan", "consensus", "plan", "discuss", "plan", "consensus", "plan",
         "monitor", "evaluate", "discuss", "plan", "consensus", "emotion", "plan", "discuss"),
  T4 = c("emotion", "consensus", "discuss", "emotion", "plan", "consensus", "discuss", "emotion",
         "plan", "consensus", "emotion", "plan", "evaluate", "consensus", "discuss", "emotion"),
  Group = c("Collaborative", "Reflective", "Strategic", "Adaptive",
            "Collaborative", "Reflective", "Strategic", "Adaptive",
            "Collaborative", "Reflective", "Strategic", "Adaptive",
            "Collaborative", "Reflective", "Strategic", "Adaptive")
)

cat("=== MULTI-GROUP ANALYSIS DEMO ===\n\n")
cat("Data overview:\n")
print(demo_data)
cat("\nGroups:", paste(unique(demo_data$Group), collapse = ", "), "\n")
cat("Group sizes:", table(demo_data$Group), "\n\n")

# =============================================================================
# 1. MULTI-GROUP PATTERN ANALYSIS
# =============================================================================

cat("1. MULTI-GROUP PATTERN ANALYSIS\n")
cat("================================\n")

pattern_result <- analyze_patterns_multi(demo_data, 
                                        group_col = "Group",
                                        min_length = 2,
                                        max_length = 3)

cat("\nPattern Analysis Results:\n")
print(pattern_result)

cat("\nDetailed Summary:\n")
summary(pattern_result, measure = "support")

# =============================================================================
# 2. MULTI-GROUP SEQUENCE COMPARISON
# =============================================================================

cat("\n\n2. MULTI-GROUP SEQUENCE COMPARISON\n")
cat("===================================\n")

sequence_result <- compare_sequences_multi(demo_data, 
                                          "Group",
                                          min_length = 2,
                                          max_length = 3,
                                          top_n = 10)

cat("\nSequence Comparison Results:\n")
print(sequence_result)

cat("\nDetailed Summary:\n")
summary(sequence_result)

# =============================================================================
# 3. SEQUENCE INDICES (works with any number of groups)
# =============================================================================

cat("\n\n3. SEQUENCE INDICES ANALYSIS\n")
cat("=============================\n")

# Define favorable states for learning
favorable_states <- c("consensus", "plan", "discuss", "evaluate")

indices_result <- compute_sequence_indices(demo_data,
                                         group_col = "Group",
                                         favorable_states = favorable_states,
                                         return_group_summary = TRUE)

cat("\nSequence Indices Results:\n")
print_indices_summary(indices_result)

# =============================================================================
# 4. BACKWARD COMPATIBILITY DEMO
# =============================================================================

cat("\n\n4. BACKWARD COMPATIBILITY\n")
cat("==========================\n")

# Create 2-group subset
two_group_data <- demo_data[demo_data$Group %in% c("Collaborative", "Reflective"), ]

cat("Using original analyze_patterns() with 2 groups:\n")
original_result <- analyze_patterns(two_group_data, group_col = "Group")
print(original_result)

cat("\nUsing original analyze_patterns() with 4 groups (shows warning):\n")
multi_with_warning <- analyze_patterns(demo_data, group_col = "Group")
print(multi_with_warning)

cat("\n=== DEMO COMPLETE ===\n")
cat("The package now supports both multi-group and traditional two-group analysis!\n")
cat("Use *_multi() functions for 3+ groups, original functions for detailed 2-group analysis.\n") 
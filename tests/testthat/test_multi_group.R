# =============================================================================
# TESTS FOR MULTI-GROUP FUNCTIONALITY
# =============================================================================

library(testthat)

# Create test data with 3 groups
test_data_multi <- data.frame(
  T1 = c("Active", "Average", "Disengaged", "Active", "Average", "Disengaged"),
  T2 = c("Active", "Active", "Disengaged", "Average", "Average", "Active"),
  T3 = c("Average", "Active", "Average", "Average", "Active", "Disengaged"),
  Group = c("A", "B", "C", "A", "B", "C")
)

test_that("analyze_patterns_multi works with 3 groups", {
  result <- analyze_patterns_multi(test_data_multi, group_col = "Group")
  
  expect_s3_class(result, "pattern_analysis_multi")
  expect_true("support" %in% names(result))
  expect_true("metadata" %in% names(result))
  expect_equal(result$metadata$n_groups, 3)
  expect_equal(length(result$metadata$group_names), 3)
  expect_true(all(c("A", "B", "C") %in% result$metadata$group_names))
})

test_that("compare_sequences_multi works with 3 groups", {
  result <- compare_sequences_multi(test_data_multi, "Group", min_length = 2, max_length = 3)
  
  expect_s3_class(result, "compare_sequences_multi")
  expect_true("results" %in% names(result))
  expect_true("summary" %in% names(result))
  expect_true("stats" %in% names(result))
  expect_equal(result$stats$n_groups, 3)
  expect_true(all(c("A", "B", "C") %in% result$stats$group_names))
})

test_that("analyze_patterns detects multiple groups and issues warning", {
  expect_warning(
    result <- analyze_patterns(test_data_multi, group_col = "Group"),
    "More than 2 groups detected"
  )
  
  # Should still work but only use first 2 groups
  expect_s3_class(result, "pattern_analysis")
})

test_that("compare_sequences detects multiple groups and issues warning", {
  expect_warning(
    result <- compare_sequences(test_data_multi, "Group"),
    "More than 2 groups detected"
  )
  
  # Should still work but only use first 2 groups
  expect_s3_class(result, "compare_sequences")
})

test_that("print methods work for multi-group objects", {
  result_patterns <- analyze_patterns_multi(test_data_multi, group_col = "Group")
  result_sequences <- compare_sequences_multi(test_data_multi, "Group")
  
  expect_output(print(result_patterns), "Multi-Group Pattern Analysis")
  expect_output(print(result_sequences), "Multi-Group Sequence Comparison")
})

test_that("summary methods work for multi-group objects", {
  result_patterns <- analyze_patterns_multi(test_data_multi, group_col = "Group")
  result_sequences <- compare_sequences_multi(test_data_multi, "Group")
  
  expect_output(summary(result_patterns), "Multi-Group Pattern Analysis Summary")
  expect_output(summary(result_sequences), "Multi-Group Sequence Comparison - Detailed Summary")
})

test_that("compute_sequence_indices works with multiple groups", {
  result <- compute_sequence_indices(test_data_multi, 
                                   group_col = "Group",
                                   return_group_summary = TRUE)
  
  expect_true(is.list(result))
  expect_true("individual_indices" %in% names(result))
  expect_true("group_summaries" %in% names(result))
  expect_true("Group" %in% names(result$individual_indices))
  
  # Check that all groups are represented
  groups_in_data <- unique(result$individual_indices$Group)
  expect_true(all(c("A", "B", "C") %in% groups_in_data))
}) 
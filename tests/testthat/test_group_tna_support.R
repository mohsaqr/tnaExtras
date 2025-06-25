# =============================================================================
# TESTS FOR GROUP_TNA OBJECT SUPPORT
# =============================================================================

library(testthat)

# Source all required functions
source("../../R/group_tna_support.R")
source("../../R/pattern_analysis.R")
source("../../R/sequence_comparison.R")
source("../../R/sequence_indices.R")

test_that("is_group_tna correctly identifies group_tna objects", {
  # Create a mock group_tna object
  mock_tna <- create_mock_group_tna()
  
  # Create a regular data.frame
  regular_data <- data.frame(
    T1 = c("A", "B", "C"),
    T2 = c("B", "A", "A"),
    Group = c("X", "Y", "Z")
  )
  
  expect_true(is_group_tna(mock_tna))
  expect_false(is_group_tna(regular_data))
  expect_false(is_group_tna(list()))
  expect_false(is_group_tna("not an object"))
})

test_that("create_mock_group_tna creates valid group_tna objects", {
  # Test default creation
  mock_tna <- create_mock_group_tna()
  
  expect_s3_class(mock_tna, "group_tna")
  expect_equal(length(mock_tna), 3)  # Default 3 groups
  expect_equal(attr(mock_tna, "levels"), c("GroupA", "GroupB", "GroupC"))
  expect_equal(attr(mock_tna, "label"), "Treatment")
  
  # Test custom creation
  custom_tna <- create_mock_group_tna(
    groups = c("Control", "Treatment"), 
    n_sequences = 3,
    n_timepoints = 5,
    states = c("Low", "High")
  )
  
  expect_equal(length(custom_tna), 2)
  expect_equal(attr(custom_tna, "levels"), c("Control", "Treatment"))
  expect_equal(ncol(custom_tna[[1]]$data), 5)
  expect_equal(nrow(custom_tna[[1]]$data), 3)
})

test_that("extract_group_tna_info extracts metadata correctly", {
  mock_tna <- create_mock_group_tna(
    groups = c("A", "B", "C"),
    states = c("State1", "State2", "State3")
  )
  
  info <- extract_group_tna_info(mock_tna)
  
  expect_equal(info$label, "Treatment")
  expect_equal(info$levels, c("A", "B", "C"))
  expect_equal(info$na_rm, FALSE)
  expect_equal(info$n_groups, 3)
  expect_equal(info$state_labels, c("State1", "State2", "State3"))
})

test_that("convert_group_tna converts objects correctly", {
  mock_tna <- create_mock_group_tna(
    groups = c("GroupA", "GroupB"),
    n_sequences = 4,
    n_timepoints = 3
  )
  
  converted <- convert_group_tna(mock_tna)
  
  expect_true(is.list(converted))
  expect_true("data" %in% names(converted))
  expect_true("group_info" %in% names(converted))
  expect_true("group_col" %in% names(converted))
  
  # Check converted data
  expect_true(is.data.frame(converted$data))
  expect_true(".group" %in% names(converted$data))
  expect_equal(converted$group_col, ".group")
  
  # Check data structure
  expect_equal(nrow(converted$data), 8)  # 4 sequences * 2 groups
  expect_equal(ncol(converted$data), 4)  # 3 timepoints + 1 group column
  
  # Check group values
  groups_in_data <- unique(converted$data$.group)
  expect_true(all(c("GroupA", "GroupB") %in% groups_in_data))
})

test_that("convert_group_tna fails on non-group_tna objects", {
  regular_data <- data.frame(T1 = c("A", "B"), T2 = c("B", "A"))
  
  expect_error(
    convert_group_tna(regular_data),
    "Object is not a group_tna object"
  )
})

test_that("analyze_patterns_multi works with group_tna objects", {
  mock_tna <- create_mock_group_tna(
    groups = c("GroupA", "GroupB", "GroupC"),
    n_sequences = 5,
    n_timepoints = 4,
    states = c("Active", "Average", "Disengaged")
  )
  
  # Should work without errors
  result <- analyze_patterns_multi(mock_tna)
  
  expect_s3_class(result, "pattern_analysis_multi")
  expect_true("support" %in% names(result))
  expect_true("metadata" %in% names(result))
  expect_equal(result$metadata$n_groups, 3)
  expect_true(all(c("GroupA", "GroupB", "GroupC") %in% result$metadata$group_names))
  
  # Check that group_tna_info is preserved
  expect_true("group_tna_info" %in% names(result$metadata))
  expect_equal(result$metadata$group_tna_info$label, "Treatment")
})

test_that("compare_sequences_multi works with group_tna objects", {
  mock_tna <- create_mock_group_tna(
    groups = c("GroupA", "GroupB", "GroupC"),
    n_sequences = 6,
    n_timepoints = 4,
    states = c("Active", "Average", "Disengaged")
  )
  
  # Should work without errors
  result <- compare_sequences_multi(mock_tna, min_length = 2, max_length = 3)
  
  expect_s3_class(result, "compare_sequences_multi")
  expect_true("results" %in% names(result))
  expect_true("summary" %in% names(result))
  expect_true("stats" %in% names(result))
  expect_equal(result$stats$n_groups, 3)
  expect_true(all(c("GroupA", "GroupB", "GroupC") %in% result$stats$group_names))
  
  # Check that group_tna_info is preserved
  expect_true("group_tna_info" %in% names(result$stats))
  expect_equal(result$stats$group_tna_info$label, "Treatment")
})

test_that("compute_sequence_indices works with group_tna objects", {
  mock_tna <- create_mock_group_tna(
    groups = c("GroupA", "GroupB"),
    n_sequences = 4,
    n_timepoints = 5,
    states = c("Active", "Average", "Disengaged")
  )
  
  # Test with group summary
  result <- compute_sequence_indices(
    mock_tna, 
    favorable_states = c("Active"),
    return_group_summary = TRUE
  )
  
  expect_true(is.list(result))
  expect_true("individual_indices" %in% names(result))
  expect_true("group_summaries" %in% names(result))
  expect_true("parameters" %in% names(result))
  
  # Check that .group column exists in individual_indices
  expect_true(".group" %in% names(result$individual_indices))
  
  # Check that group_tna_info is preserved
  expect_true("group_tna_info" %in% names(result$parameters))
  expect_equal(result$parameters$group_tna_info$label, "Treatment")
  
  # Check groups in data
  groups_in_data <- unique(result$individual_indices$.group)
  expect_true(all(c("GroupA", "GroupB") %in% groups_in_data))
})

test_that("group_tna objects use state labels appropriately", {
  mock_tna <- create_mock_group_tna(
    groups = c("GroupA", "GroupB"),
    n_sequences = 3,
    n_timepoints = 3,
    states = c("Low", "Medium", "High")
  )
  
  # The state labels should be used when not specified - use return_group_summary to get parameters
  result <- compute_sequence_indices(mock_tna, return_group_summary = TRUE)
  
  # Should have used the state labels from group_tna
  expect_true("Low" %in% result$parameters$all_possible_states ||
              "Medium" %in% result$parameters$all_possible_states ||
              "High" %in% result$parameters$all_possible_states)
  
  # Also test that group_tna_info is preserved when returning simple results
  simple_result <- compute_sequence_indices(mock_tna)
  expect_true(!is.null(attr(simple_result, "group_tna_info")))
  expect_equal(attr(simple_result, "group_tna_info")$label, "Treatment")
})

test_that("print method works for group_tna objects", {
  mock_tna <- create_mock_group_tna()
  
  expect_output(print(mock_tna), "Mock group_tna object")
  expect_output(print(mock_tna), "Groups: 3")
  expect_output(print(mock_tna), "GroupA, GroupB, GroupC")
  expect_output(print(mock_tna), "Treatment")
})

test_that("null coalescing operator works correctly", {
  expect_equal(NULL %||% "default", "default")
  expect_equal("value" %||% "default", "value")
  expect_equal(NA %||% "default", NA)
  expect_equal(0 %||% "default", 0)
  expect_equal("" %||% "default", "")
})

test_that("backwards compatibility maintained", {
  # Test that regular data.frame input still works with all functions
  regular_data <- data.frame(
    T1 = c("Active", "Average", "Disengaged", "Active", "Average", "Disengaged"),
    T2 = c("Active", "Active", "Disengaged", "Average", "Average", "Active"),
    T3 = c("Average", "Active", "Average", "Average", "Active", "Disengaged"),
    Group = c("A", "B", "C", "A", "B", "C")
  )
  
  # All functions should still work normally
  result1 <- analyze_patterns_multi(regular_data, group_col = "Group")
  expect_s3_class(result1, "pattern_analysis_multi")
  
  result2 <- compare_sequences_multi(regular_data, "Group", min_length = 2, max_length = 3)
  expect_s3_class(result2, "compare_sequences_multi")
  
  result3 <- compute_sequence_indices(regular_data, group_col = "Group")
  expect_true(is.data.frame(result3))
})

test_that("error handling works correctly", {
  # Test with invalid input
  expect_error(is_group_tna(), "argument \"x\" is missing")
  
  # Test convert_group_tna with wrong input
  expect_error(
    convert_group_tna(data.frame(x = 1:3)),
    "Object is not a group_tna object"
  )
  
  # Test with empty group_tna-like object
  empty_obj <- list()
  class(empty_obj) <- "group_tna"
  
  expect_error(
    convert_group_tna(empty_obj),
    NA  # Should not error, but may produce warnings
  )
}) 
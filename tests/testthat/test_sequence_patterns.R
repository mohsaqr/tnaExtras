# ==============================================================================
# TESTS FOR SEQUENCE PATTERN EXPLORATION
# ==============================================================================

test_that("explore_sequence_patterns works with basic input", {
  # Create simple test data
  test_data <- data.frame(
    T1 = c("A", "B", "A", "A", "B", "A", "B", "A", "B", "A"),
    T2 = c("B", "A", "B", "B", "A", "B", "A", "B", "A", "B"),
    T3 = c("A", "B", "A", "A", "B", "A", "B", "A", "B", "A"),
    T4 = c("B", "A", "B", "B", "A", "B", "A", "B", "A", "B"),
    stringsAsFactors = FALSE
  )
  
  result <- explore_sequence_patterns(
    test_data, 
    min_length = 1, 
    max_length = 2,
    min_support = 0.1,
    min_count = 1,
    verbose = FALSE
  )
  
  expect_s3_class(result, "sequence_pattern_analysis")
  expect_true("patterns" %in% names(result))
  expect_true("full_sequences" %in% names(result))
  expect_true("state_frequencies" %in% names(result))
  expect_true("transitions" %in% names(result))
  expect_true("summary" %in% names(result))
  expect_true("parameters" %in% names(result))
})

test_that("explore_sequence_patterns computes correct state frequencies", {
  test_data <- data.frame(
    T1 = c("A", "A", "A", "A", "A"),
    T2 = c("B", "B", "B", "B", "B"),
    T3 = c("A", "A", "A", "A", "A"),
    stringsAsFactors = FALSE
  )
  
  result <- explore_sequence_patterns(
    test_data,
    min_length = 1,
    max_length = 1,
    min_support = 0,
    min_count = 1,
    test_significance = FALSE,
    verbose = FALSE
  )
  
  # Should have 10 A's and 5 B's
  state_freq <- result$state_frequencies
  a_freq <- state_freq$count[state_freq$state == "A"]
  b_freq <- state_freq$count[state_freq$state == "B"]
  
  expect_equal(a_freq, 10)
  expect_equal(b_freq, 5)
})

test_that("explore_sequence_patterns respects min_support threshold", {
  test_data <- data.frame(
    T1 = c("A", "A", "A", "A", "A", "X", "Y", "Z", "W", "V"),
    T2 = c("B", "B", "B", "B", "B", "A", "A", "A", "A", "A"),
    stringsAsFactors = FALSE
  )
  
  result <- explore_sequence_patterns(
    test_data,
    min_length = 1,
    max_length = 2,
    min_support = 0.3,  # 30% = at least 3 sequences
    min_count = 1,
    test_significance = FALSE,
    verbose = FALSE
  )
  
  # All patterns should have support >= 0.3
  expect_true(all(result$patterns$support >= 0.3))
})

test_that("explore_sequence_patterns respects min_count threshold", {
  test_data <- data.frame(
    T1 = c("A", "A", "B", "C", "D"),
    T2 = c("B", "B", "C", "D", "E"),
    stringsAsFactors = FALSE
  )
  
  result <- explore_sequence_patterns(
    test_data,
    min_length = 2,
    max_length = 2,
    min_support = 0,
    min_count = 2,  # At least 2 occurrences
    test_significance = FALSE,
    verbose = FALSE
  )
  
  # All patterns should have count >= 2
  if (nrow(result$patterns) > 0) {
    expect_true(all(result$patterns$count >= 2))
  }
})

test_that("explore_sequence_patterns performs significance testing", {
  test_data <- data.frame(
    T1 = c("A", "A", "A", "A", "A", "A", "A", "A", "A", "A"),
    T2 = c("B", "B", "B", "B", "B", "B", "B", "B", "B", "B"),
    T3 = c("A", "A", "A", "A", "A", "A", "A", "A", "A", "A"),
    stringsAsFactors = FALSE
  )
  
  result <- explore_sequence_patterns(
    test_data,
    min_length = 1,
    max_length = 2,
    min_support = 0.1,
    min_count = 1,
    test_significance = TRUE,
    correction = "bonferroni",
    alpha = 0.05,
    verbose = FALSE
  )
  
  expect_true("p_value" %in% names(result$patterns))
  expect_true("p_adjusted" %in% names(result$patterns))
  expect_true("significant" %in% names(result$patterns))
})

test_that("significant_patterns returns only significant patterns", {
  test_data <- data.frame(
    T1 = c("A", "A", "A", "A", "A", "A", "A", "A", "A", "A"),
    T2 = c("B", "B", "B", "B", "B", "B", "B", "B", "B", "B"),
    stringsAsFactors = FALSE
  )
  
  result <- explore_sequence_patterns(
    test_data,
    min_length = 1,
    max_length = 2,
    min_support = 0.1,
    min_count = 1,
    test_significance = TRUE,
    verbose = FALSE
  )
  
  sig <- significant_patterns(result)
  
  if (nrow(sig) > 0) {
    expect_true(all(sig$significant == TRUE))
  }
})

test_that("filter_patterns works correctly", {
  test_data <- data.frame(
    T1 = c("A", "A", "A", "A", "A", "A", "A", "A", "A", "A"),
    T2 = c("B", "B", "B", "B", "B", "B", "B", "B", "B", "B"),
    T3 = c("A", "A", "A", "A", "A", "A", "A", "A", "A", "A"),
    stringsAsFactors = FALSE
  )
  
  result <- explore_sequence_patterns(
    test_data,
    min_length = 1,
    max_length = 3,
    min_support = 0.1,
    min_count = 1,
    verbose = FALSE
  )
  
  # Filter by length
  filtered <- filter_patterns(result, pattern_length = 2)
  if (nrow(filtered) > 0) {
    expect_true(all(filtered$length == 2))
  }
  
  # Filter by support
  filtered <- filter_patterns(result, min_support = 0.5)
  if (nrow(filtered) > 0) {
    expect_true(all(filtered$support >= 0.5))
  }
})

test_that("get_transition_matrix returns correct format", {
  test_data <- data.frame(
    T1 = c("A", "A", "B"),
    T2 = c("B", "B", "A"),
    T3 = c("A", "A", "B"),
    stringsAsFactors = FALSE
  )
  
  result <- explore_sequence_patterns(
    test_data,
    min_length = 1,
    max_length = 2,
    min_support = 0,
    min_count = 1,
    test_significance = FALSE,
    verbose = FALSE
  )
  
  # Get count matrix
  count_mat <- get_transition_matrix(result, type = "count")
  expect_true(is.matrix(count_mat))
  expect_equal(nrow(count_mat), ncol(count_mat))
  
  # Get probability matrix
  prob_mat <- get_transition_matrix(result, type = "probability")
  expect_true(is.matrix(prob_mat))
  
  # Probabilities should sum to 1 for each row (or 0 if no transitions)
  row_sums <- rowSums(prob_mat)
  for (rs in row_sums) {
    expect_true(rs == 0 || abs(rs - 1) < 0.01)
  }
})

test_that("explore_sequence_patterns handles missing values", {
  test_data <- data.frame(
    T1 = c("A", "A", "A", NA, "A"),
    T2 = c("B", NA, "B", "B", "B"),
    T3 = c("A", "A", NA, "A", "A"),
    stringsAsFactors = FALSE
  )
  
  result <- explore_sequence_patterns(
    test_data,
    min_length = 1,
    max_length = 2,
    min_support = 0,
    min_count = 1,
    test_significance = FALSE,
    verbose = FALSE
  )
  
  expect_s3_class(result, "sequence_pattern_analysis")
  expect_true(result$summary$n_sequences_valid > 0)
})

test_that("explore_sequence_patterns validates parameters", {
  test_data <- data.frame(
    T1 = c("A", "B", "A"),
    T2 = c("B", "A", "B"),
    stringsAsFactors = FALSE
  )
  
  # Invalid min_length
  expect_error(
    explore_sequence_patterns(test_data, min_length = 0, verbose = FALSE)
  )
  
  # Invalid max_length
  expect_error(
    explore_sequence_patterns(test_data, min_length = 3, max_length = 2, verbose = FALSE)
  )
  
  # Invalid min_support
  expect_error(
    explore_sequence_patterns(test_data, min_support = -0.1, verbose = FALSE)
  )
  
  # Invalid correction method
  expect_error(
    explore_sequence_patterns(test_data, correction = "invalid", verbose = FALSE)
  )
})

test_that("print and summary methods work", {
  test_data <- data.frame(
    T1 = c("A", "A", "A", "A", "A"),
    T2 = c("B", "B", "B", "B", "B"),
    stringsAsFactors = FALSE
  )
  
  result <- explore_sequence_patterns(
    test_data,
    min_length = 1,
    max_length = 2,
    min_support = 0.1,
    min_count = 1,
    verbose = FALSE
  )
  
  # Print should work without error
  expect_output(print(result))
  
  # Summary should work without error
  expect_output(summary(result))
})

test_that("explore_sequence_patterns works with real data", {
  skip_if_not_installed("tnaExtras")
  
  # Load sample data
  data(seqdata, package = "tnaExtras")
  
  # Remove group column
  seq_data <- seqdata[, -1]
  
  result <- explore_sequence_patterns(
    seq_data,
    min_length = 2,
    max_length = 3,
    min_support = 0.05,
    min_count = 3,
    verbose = FALSE
  )
  
  expect_s3_class(result, "sequence_pattern_analysis")
  expect_true(nrow(result$patterns) > 0)
  expect_true(result$summary$n_significant_patterns > 0)
})


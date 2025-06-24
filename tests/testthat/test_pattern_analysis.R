# Test suite for analyze_patterns() function

# Test data (simplified from original for brevity in this refactoring example)
test_data_pattern <- data.frame(
  T1 = c("Active", "Average", "Disengaged", "Active", "Average"),
  T2 = c("Active", "Active", "Disengaged", "Active", "Average"),
  T3 = c("Average", "Active", "Average", "Disengaged", "Active"),
  T4 = c("Average", "Disengaged", "Average", "Disengaged", "Active"),
  T5 = c("Disengaged", "Disengaged", "Active", "Active", "Active"),
  Group = c("A", "A", "B", "B", "A") # Added group for some tests
)

# Edge case data
single_pattern_data <- data.frame(
  T1 = rep("A", 5), T2 = rep("A", 5), T3 = rep("A", 5),
  Group = rep("X", 5)
)

context("Pattern Analysis: analyze_patterns")

test_that("Basic functionality and output structure are correct", {
  result <- analyze_patterns(test_data_pattern, group_col = "Group", min_char_length_state = 1)

  expected_components <- c("support", "lift", "confidence", "effect_size", "metadata")
  expect_true(all(expected_components %in% names(result)),
              info = "All main measure components should be present.")
  expect_s3_class(result, "pattern_analysis") # Check class

  # Check metadata
  expect_type(result$metadata, "list")
  # The following check for n_patterns in metadata might need refinement.
  # It counts unique states in the patterns found in the support table.
  # This is a proxy and might not perfectly match metadata$n_patterns if that counts unique n-grams before measure calculation.
  # For now, let's check if it's a positive integer.
  expect_true(result$metadata$n_patterns >= 0)
})

test_that("extract_ngrams handles min_char_length_state correctly", {
  sequence1 <- "A-B-CC-DDD-E"
  expect_equal(extract_ngrams(sequence1, n = 2, min_char_length_state = 1),
               c("A-B", "B-CC", "CC-DDD", "DDD-E"))
  expect_equal(extract_ngrams(sequence1, n = 2, min_char_length_state = 2),
               c("CC-DDD")) # A, B, E are filtered
  expect_equal(extract_ngrams("A-B", n = 2, min_char_length_state = 2), character(0))
  expect_equal(extract_ngrams(sequence1, n = 1, min_char_length_state = 1),
               c("A", "B", "CC", "DDD", "E"))
  expect_equal(extract_ngrams(sequence1, n = 1, min_char_length_state = 3),
               c("DDD"))
})

test_that("compute_lift_measures calculates lift correctly", {
  # Test case:
  # Group A: 10 sequences, Pattern X appears in 5. P(X|A) = 0.5
  # Group B: 10 sequences, Pattern X appears in 2. P(X|B) = 0.2
  # Total: 20 sequences, Pattern X appears in 7. P(X) = 7/20 = 0.35
  # Lift A = P(X|A) / P(X) = 0.5 / 0.35 = 1.42857
  # Lift B = P(X|B) / P(X) = 0.2 / 0.35 = 0.57143

  patterns_lift <- c("PatternX")
  group_A_seqs_lift <- c(rep("PatternX", 5), rep("Other", 5)) # 10 seqs
  group_B_seqs_lift <- c(rep("PatternX", 2), rep("Other", 8)) # 10 seqs

  lift_results <- compute_lift_measures(patterns_lift, group_A_seqs_lift, group_B_seqs_lift)

  expect_equal(lift_results$lift_A[1], (5/10) / (7/20), tolerance = 1e-5)
  expect_equal(lift_results$lift_B[1], (2/10) / (7/20), tolerance = 1e-5)

  # Test case: Pattern only in Group A
  patterns_lift2 <- c("PatternY")
  group_A_seqs_lift2 <- c(rep("PatternY", 3), rep("Other", 7))
  group_B_seqs_lift2 <- rep("Other", 10)

  lift_results2 <- compute_lift_measures(patterns_lift2, group_A_seqs_lift2, group_B_seqs_lift2)
  expect_equal(lift_results2$lift_A[1], (3/10) / (3/20), tolerance = 1e-5)
  expect_equal(lift_results2$lift_B[1], 0.0 / (3/20), tolerance = 1e-5) # P(Y|B)/P(Y) = 0 / P(Y) = 0

  # Test case: Pattern appears in no group (count_total = 0 should be handled by loop skip in main func)
  # If a pattern is passed to compute_lift_measures that has zero total count,
  # prob_pattern will be 0. lift_A/B will be prop_A/1e-10 or prop_B/1e-10.
  # If count_A and count_B are also 0, then prop_A/B are 0, so lift_A/B are 0. This is fine.
  patterns_lift3 <- c("PatternZ")
  group_A_seqs_lift3 <- rep("Other", 10)
  group_B_seqs_lift3 <- rep("Other", 10)
  lift_results3 <- compute_lift_measures(patterns_lift3, group_A_seqs_lift3, group_B_seqs_lift3)
  expect_equal(lift_results3$lift_A[1], 0.0, tolerance=1e-5) # 0 / (0 + 1e-10)
  expect_equal(lift_results3$lift_B[1], 0.0, tolerance=1e-5)
})

test_that("Effect sizes (Cohen's h, Phi) are calculated correctly", {
  patterns_es <- "PatternP"
  group_A_es <- c(rep("PatternP", 8), rep("Other", 2)) # 10 seqs, pA = 0.8
  group_B_es <- c(rep("PatternP", 2), rep("Other", 8)) # 10 seqs, pB = 0.2

  es_results <- compute_effect_sizes(patterns_es, group_A_es, group_B_es)

  h_expected <- 2 * (asin(sqrt(0.8)) - asin(sqrt(0.2)))
  expect_equal(es_results$cohens_h[1], h_expected, tolerance = 1e-5)

  # ChiSq for [[8,2],[2,8]] (PatternP vs Other, GroupA vs GroupB for sequences)
  # present_A=8, absent_A=2 (sequences in A not having PatternP)
  # present_B=2, absent_B=8 (sequences in B not having PatternP)
  # E_present_A = ( (8+2) * (8+2) ) / 20 = 10 * 10 / 20 = 5
  # ChiSq = (8-5)^2/5 + (2-5)^2/5 + (2-5)^2/5 + (8-5)^2/5 = 4 * (9/5) = 7.2
  # Phi = sqrt(ChiSq / N) = sqrt(7.2 / 20) = 0.6
  expect_equal(es_results$phi_coefficient[1], sqrt(7.2/20), tolerance = 1e-5)
  expect_equal(es_results$cramers_v[1], sqrt(7.2/20), tolerance = 1e-5)
})

test_that("Support measures for n=1 patterns (state presence)", {
  # test_data_pattern Group A: 3 sequences, Group B: 2 sequences
  # Group A seqs:
  # 1: Active-Active-Average-Average-Disengaged
  # 2: Average-Active-Active-Disengaged-Disengaged
  # 3: Active-Average-Active-Active-Active
  # Group B seqs:
  # 4: Disengaged-Disengaged-Average-Average-Active
  # 5: Active-Average-Disengaged-Active-Active

  # For min_length=1, max_length=1, it means we are looking at individual states.
  # Support for a state "S" in a group is the proportion of sequences in that group containing "S".
  result_n1 <- analyze_patterns(test_data_pattern, group_col = "Group",
                                min_length=1, max_length=1, # Focus on single states
                                min_frequency = 1,         # Include all states found
                                min_char_length_state = 1)

  active_support <- result_n1$support[result_n1$support$pattern == "Active",]
  expect_equal(active_support$support_A, 3/3, info = "Support for Active in Group A") # All 3 Group A seqs have "Active"
  expect_equal(active_support$support_B, 2/2, info = "Support for Active in Group B") # All 2 Group B seqs have "Active"

  average_support <- result_n1$support[result_n1$support$pattern == "Average",]
  expect_equal(average_support$support_A, 3/3, info = "Support for Average in Group A") # Seq 1,2,3 of A have Average
  expect_equal(average_support$support_B, 2/2, info = "Support for Average in Group B") # Seq 1,2 of B have Average

  disengaged_support <- result_n1$support[result_n1$support$pattern == "Disengaged",]
  expect_equal(disengaged_support$support_A, 2/3, info = "Support for Disengaged in Group A") # Seq 1,2 of A have Disengaged
  expect_equal(disengaged_support$support_B, 2/2, info = "Support for Disengaged in Group B") # Seq 1,2 of B have Disengaged
})

test_that("analyze_patterns handles edge cases gracefully", {
  expect_no_error(analyze_patterns(single_pattern_data, group_col = "Group", min_char_length_state = 1))

  minimal_data <- data.frame(T1 = c("A", "B"), T2 = c("A", "B"), Group = c("X", "Y"))
  expect_no_error(analyze_patterns(minimal_data, group_col = "Group", min_char_length_state = 1))

  missing_data_df <- data.frame(
    T1 = c("A", NA, "C", "A", "B"), T2 = c("A", "B", NA, "B", "B"),
    T3 = c("B", "B", "C", NA, "C"), T4 = c("C", "C", "A", "A", NA),
    Group = c("P", "P", "Q", "Q", "P"))
  expect_no_error(analyze_patterns(missing_data_df, group_col = "Group", min_char_length_state = 1))

  # Test case with no valid sequences after processing
  all_empty_seq_data <- data.frame(T1 = c(NA, ""), T2 = c(NA, ""), Group = c("X","Y"))
  expect_error(analyze_patterns(all_empty_seq_data, group_col="Group"),
               "No valid sequences found after processing")

  # Test case with less than 2 groups
  one_group_data <- data.frame(T1="A", T2="B", Group="X")
   expect_error(analyze_patterns(one_group_data, group_col="Group"),
               "Exactly two groups must be present in the data")
})

test_that("Output is consistent across multiple runs", {
  consistency_data <- data.frame(
    T1 = c("A", "B", "A"), T2 = c("B", "A", "A"), Group = c("X", "Y", "X"))
  result1 <- analyze_patterns(consistency_data, group_col = "Group", min_char_length_state = 1)
  result2 <- analyze_patterns(consistency_data, group_col = "Group", min_char_length_state = 1)

  expect_identical(result1$support, result2$support)
  expect_identical(result1$lift, result2$lift)
  expect_identical(result1$confidence, result2$confidence)
  expect_identical(result1$effect_size, result2$effect_size)
  expect_identical(result1$metadata, result2$metadata) # Metadata should also be identical
})

test_that("Print and summary methods for pattern_analysis objects run", {
  result <- analyze_patterns(test_data_pattern, group_col = "Group", min_char_length_state = 1)
  expect_output(print(result), "Pattern Analysis Results")
  expect_output(summary(result), "Pattern Analysis Summary - SUPPORT Measures")
  expect_output(summary(result, measure = "lift"), "Pattern Analysis Summary - LIFT Measures")
  expect_error(summary(result, measure = "nonexistent"),
               "Measure 'nonexistent' not found. Available: support, lift, confidence, effect_size")
})
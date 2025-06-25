library(testthat)
library(tnaExtras)

# Helper function to create sample sequence data (can be shared if in helper-file)
create_sample_pattern_data <- function(n_rows, n_cols, groups_vec) {
  set.seed(456) # for reproducibility
  data_matrix <- matrix(sample(c("s1", "s2", "s3", "s4", "s5"), n_rows * n_cols, replace = TRUE),
                        nrow = n_rows, ncol = n_cols)
  data_df <- as.data.frame(data_matrix)
  names(data_df) <- paste0("Time", 1:n_cols)
  data_df$GroupCol <- groups_vec
  return(data_df)
}

test_that("analyze_patterns_multiple works with 3 groups", {
  groups_3 <- sample(c("GroupA", "GroupB", "GroupC"), 30, replace = TRUE)
  data_3_groups <- create_sample_pattern_data(n_rows = 30, n_cols = 4, groups_vec = groups_3)

  result <- analyze_patterns_multiple(data_3_groups, group_col = "GroupCol",
                                    measures = c("support", "lift"))

  expect_s3_class(result, "analyze_patterns_multiple_result")
  expect_length(result, 3) # Should be 3 pairs

  # Check names (order might vary based on factor levels)
  expected_pairs <- c("GroupA_vs_GroupB", "GroupA_vs_GroupC", "GroupB_vs_GroupC")
  actual_pairs <- names(result)
  expect_true(all(sort(expected_pairs) == sort(actual_pairs)))


  for (pair_name in names(result)) {
    expect_s3_class(result[[pair_name]], "pattern_analysis",
                    info = paste("Checking pair:", pair_name))
    expect_true("support" %in% names(result[[pair_name]]),
                info = paste("Support check for pair:", pair_name))
    expect_true("lift" %in% names(result[[pair_name]]),
                info = paste("Lift check for pair:", pair_name))
    expect_false("confidence" %in% names(result[[pair_name]]),
                 info = paste("Confidence check for pair:", pair_name)) # Since not requested
  }
})

test_that("analyze_patterns_multiple works with 2 groups (defers to analyze_patterns)", {
  groups_2 <- sample(c("X1", "X2"), 20, replace = TRUE)
  data_2_groups <- create_sample_pattern_data(n_rows = 20, n_cols = 3, groups_vec = groups_2)

  result_multi <- analyze_patterns_multiple(data_2_groups, group_col = "GroupCol",
                                          min_frequency = 1) # Lower freq for small sample

  expect_s3_class(result_multi, "analyze_patterns_multiple_result")
  expect_length(result_multi, 1)

  group_levels <- levels(as.factor(data_2_groups$GroupCol))
  expected_pair_name <- paste(group_levels[1], "vs", group_levels[2], sep = "_")
  expect_named(result_multi, expected_pair_name)

  expect_s3_class(result_multi[[expected_pair_name]], "pattern_analysis")
  expect_true(all(c("support", "lift", "confidence", "effect_size") %in% names(result_multi[[expected_pair_name]])))
})

test_that("analyze_patterns_multiple handles errors correctly", {
  groups_1 <- rep("Solo", 10)
  data_1_group <- create_sample_pattern_data(n_rows = 10, n_cols = 3, groups_vec = groups_1)

  expect_error(analyze_patterns_multiple(data_1_group, group_col = "GroupCol"),
               "At least two groups are required for pattern analysis. Found 1 groups.")

  data_valid <- create_sample_pattern_data(10, 3, sample(c("Y","Z"),10,replace=T))
  expect_error(analyze_patterns_multiple(data_valid, group_col = "NonExistent"),
               "'group_col' must be a single valid column name in 'data'.")

  expect_error(analyze_patterns_multiple(as.matrix(data_valid), group_col = "GroupCol"),
               "'data' must be a data frame.")
})

test_that("analyze_patterns_multiple passes parameters down to analyze_patterns", {
  groups_3 <- sample(c("Cat1", "Cat2", "Cat3"), 45, replace = TRUE)
  # Ensure enough data for min_frequency to be met for some patterns
  data_3_groups <- create_sample_pattern_data(n_rows = 45, n_cols = 5, groups_vec = groups_3)
  data_3_groups[1:15, 1:2] <- "pA-pB" # Plant a common pattern

  result <- analyze_patterns_multiple(data_3_groups, group_col = "GroupCol",
                                    min_length = 2, max_length = 3,
                                    min_frequency = 3, measures = "support")

  expect_length(result, 3)
  for (pair_name in names(result)) {
    pair_res <- result[[pair_name]]
    expect_s3_class(pair_res, "pattern_analysis")
    if (!is.null(pair_res)) { # If a pair didn't error out
        expect_equal(pair_res$metadata$pattern_length_range, c(2, 3))
        expect_equal(pair_res$metadata$min_pattern_frequency, 3)
        expect_true("support" %in% names(pair_res))
        expect_false("lift" %in% names(pair_res))
        # Check if dynamic column names for support are present
        group_names_in_pair <- pair_res$metadata$group_names
        expect_true(paste0("support_", group_names_in_pair[1]) %in% names(pair_res$support))
        expect_true(paste0("support_", group_names_in_pair[2]) %in% names(pair_res$support))
    }
  }
})


# Test with a pair that might fail due to no common patterns or insufficient data for a pair
test_that("analyze_patterns_multiple handles pairs with no common patterns or errors gracefully", {
  # Group G4 has very different patterns and fewer sequences
  data_problematic <- data.frame(
    T1 = c(rep("s1", 5), rep("s2", 5), rep("s3", 5), "x1", "x1"),
    T2 = c(rep("s1", 5), rep("s2", 5), rep("s3", 5), "x2", "x2"),
    GroupCol = c(rep("G1", 5), rep("G2", 5), rep("G3", 5), rep("G4", 2))
  )

  # Using min_frequency = 1 to ensure some patterns are found within G1, G2, G3
  # G4 sequences are "x1-x2", others are "s1-s1", "s2-s2", "s3-s3"
  # Comparisons involving G4 might lead to "No patterns meet the minimum frequency threshold"
  # if min_frequency is higher for the pairwise analysis, or no common patterns.

  # The error is caught by analyze_patterns, and analyze_patterns_multiple stores NULL.
  expect_warning({
    result <- analyze_patterns_multiple(data_problematic, group_col = "GroupCol", min_frequency = 1, min_length=2)
  }, regexp="Error during pattern analysis of G1_vs_G4.*Skipping this pair.|Error during pattern analysis of G2_vs_G4.*Skipping this pair.|Error during pattern analysis of G3_vs_G4.*Skipping this pair.")

  expect_s3_class(result, "analyze_patterns_multiple_result")
  expect_length(result, 6) # 4C2 = 6 pairs

  # Valid pairs should have results
  expect_s3_class(result[["G1_vs_G2"]], "pattern_analysis")
  expect_s3_class(result[["G1_vs_G3"]], "pattern_analysis")
  expect_s3_class(result[["G2_vs_G3"]], "pattern_analysis")

  # Pairs involving G4 should be NULL due to errors in analyze_patterns
  expect_null(result[["G1_vs_G4"]])
  expect_null(result[["G2_vs_G4"]])
  expect_null(result[["G3_vs_G4"]])
})

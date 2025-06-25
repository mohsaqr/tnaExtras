library(testthat)
library(tnaExtras) # Assuming your package name is tnaExtras

# Helper function to create sample sequence data
create_sample_seq_data <- function(n_rows, n_cols, groups_vec) {
  set.seed(123) # for reproducibility
  data_matrix <- matrix(sample(c("A", "B", "C", "D", NA), n_rows * n_cols, replace = TRUE, prob = c(0.25, 0.25, 0.25, 0.15, 0.1)),
                        nrow = n_rows, ncol = n_cols)
  data_df <- as.data.frame(data_matrix)
  names(data_df) <- paste0("T", 1:n_cols)
  data_df$Group <- groups_vec
  return(data_df)
}

test_that("compare_sequences_multiple works with 3 groups", {
  groups_3 <- sample(c("G1", "G2", "G3"), 30, replace = TRUE)
  data_3_groups <- create_sample_seq_data(n_rows = 30, n_cols = 5, groups_vec = groups_3)

  result <- compare_sequences_multiple(data_3_groups, group = "Group", statistical = FALSE)

  expect_s3_class(result, "compare_sequences_multiple_result")
  expect_length(result, 3) # Should be 3 pairs: G1_vs_G2, G1_vs_G3, G2_vs_G3
  expect_named(result, c("G1_vs_G2", "G1_vs_G3", "G2_vs_G3")) # Order might vary based on levels

  for (pair_name in names(result)) {
    expect_s3_class(result[[pair_name]], "compare_sequences",
                    info = paste("Checking pair:", pair_name))
    expect_true(is.list(result[[pair_name]]$summary) || is.data.frame(result[[pair_name]]$summary),
                info = paste("Summary check for pair:", pair_name))
  }
})

test_that("compare_sequences_multiple works with 2 groups (defers to compare_sequences)", {
  groups_2 <- sample(c("Alpha", "Beta"), 20, replace = TRUE)
  data_2_groups <- create_sample_seq_data(n_rows = 20, n_cols = 4, groups_vec = groups_2)

  result_multi <- compare_sequences_multiple(data_2_groups, group = "Group", statistical = TRUE, top_n = 3)

  expect_s3_class(result_multi, "compare_sequences_multiple_result")
  expect_length(result_multi, 1)

  # Determine expected pair name based on factor levels
  group_levels <- levels(as.factor(data_2_groups$Group))
  expected_pair_name <- paste(group_levels[1], "vs", group_levels[2], sep = "_")
  expect_named(result_multi, expected_pair_name)

  expect_s3_class(result_multi[[expected_pair_name]], "compare_sequences")
  expect_equal(result_multi[[expected_pair_name]]$parameters$top_n, 3)
})

test_that("compare_sequences_multiple handles errors correctly", {
  groups_1 <- rep("G1", 10)
  data_1_group <- create_sample_seq_data(n_rows = 10, n_cols = 3, groups_vec = groups_1)

  expect_error(compare_sequences_multiple(data_1_group, group = "Group"),
               "Number of unique groups (1) is less than the minimum required (2).")

  data_valid <- create_sample_seq_data(10, 3, sample(c("A","B"),10,replace=T))
  expect_error(compare_sequences_multiple(data_valid, group = "NonExistentGroup"),
               "Group column 'NonExistentGroup' not found in 'data'.")

  expect_error(compare_sequences_multiple(as.matrix(data_valid), group = "Group"),
               "'data' must be a data frame.")
})

test_that("compare_sequences_multiple passes parameters down to compare_sequences", {
  groups_3 <- sample(c("X", "Y", "Z"), 30, replace = TRUE)
  data_3_groups <- create_sample_seq_data(n_rows = 30, n_cols = 5, groups_vec = groups_3)

  result <- compare_sequences_multiple(data_3_groups, group = "Group",
                                     min_length = 3, max_length = 4,
                                     statistical = TRUE, correction = "holm", top_n = 5)

  expect_length(result, 3)
  for (pair_name in names(result)) {
    pair_res <- result[[pair_name]]
    expect_s3_class(pair_res, "compare_sequences")
    expect_equal(pair_res$parameters$min_length, 3)
    expect_equal(pair_res$parameters$max_length, 4)
    expect_true(pair_res$parameters$statistical)
    expect_equal(pair_res$parameters$correction, "holm")
    expect_equal(pair_res$parameters$top_n, 5)
  }
})

test_that("compare_sequences_multiple works with group vector input", {
  n_r <- 25
  data_df_no_group <- create_sample_seq_data(n_rows = n_r, n_cols = 3, groups_vec = rep("dummy", n_r))
  data_df_no_group$Group <- NULL # Remove dummy group column

  group_v <- sample(c("P1", "P2", "P3"), n_r, replace = TRUE)

  result <- compare_sequences_multiple(data_df_no_group, group = group_v, statistical = FALSE)

  expect_s3_class(result, "compare_sequences_multiple_result")
  expect_length(result, 3) # P1_vs_P2, P1_vs_P3, P2_vs_P3

  # Check one of the pairs
  pair_levels <- levels(as.factor(group_v))
  expected_name_example <- paste(pair_levels[1], "vs", pair_levels[2], sep = "_")
  if(expected_name_example %in% names(result)){
    expect_s3_class(result[[expected_name_example]], "compare_sequences")
  } else {
      # Fallback if level order is different than expected by simple indexing
      expect_true(length(names(result)[grepl(paste0(pair_levels[1], "_vs_", pair_levels[2]), names(result)) |
                                      grepl(paste0(pair_levels[2], "_vs_", pair_levels[1]), names(result))]) > 0)
  }
})

# Test with a pair that might fail to ensure overall function completes
test_that("compare_sequences_multiple handles failing pairs gracefully", {
  data_mixed <- data.frame(
    T1 = c("A", "A", "B", "B", "C", "C", "D"),
    T2 = c("A", "A", "B", "B", "C", "C", "D"),
    Group = c("G1", "G1", "G2", "G2", "G3", "G3", "G4") # G4 has only 1 sequence
  )

  # G4 will likely cause issues in pairwise comparisons due to insufficient data
  # Expect warnings, but the function should complete for other pairs.
  expect_warning({
    result <- compare_sequences_multiple(data_mixed, group = "Group", statistical = FALSE)
  }, regexp = "Error during comparison of G1_vs_G4.*Skipping this pair.|Error during comparison of G2_vs_G4.*Skipping this pair.|Error during comparison of G3_vs_G4.*Skipping this pair.") # Regex for any failing pair with G4

  expect_s3_class(result, "compare_sequences_multiple_result")
  expect_length(result, 6) # G1vG2, G1vG3, G1vG4, G2vG3, G2vG4, G3vG4

  # Check that valid pairs have results
  expect_s3_class(result[["G1_vs_G2"]], "compare_sequences")
  expect_s3_class(result[["G1_vs_G3"]], "compare_sequences")
  expect_s3_class(result[["G2_vs_G3"]], "compare_sequences")

  # Check that failed pairs might be NULL (or contain error objects if we change error handling)
  # Current implementation stores NULL for failed pairs.
  expect_null(result[["G1_vs_G4"]]) # Or any pair involving G4
  expect_null(result[["G2_vs_G4"]])
  expect_null(result[["G3_vs_G4"]])
})

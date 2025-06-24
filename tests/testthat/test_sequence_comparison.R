# Test suite for compare_sequences() in R/sequence_comparison.R

context("Sequence Comparison: compare_sequences")

# Test data
test_data_comp <- data.frame(
  T1 = c("Active", "Average", "Disengaged", "Active", "Average", "Disengaged"),
  T2 = c("Active", "Active", "Disengaged", "Active", "Average", "Disengaged"),
  T3 = c("Average", "Active", "Average", "Disengaged", "Active", "Average"),
  Group = c("A", "A", "B", "B", "A", "B")
)

test_that("compare_sequences basic functionality (discrimination mode)", {
  result <- compare_sequences(
    data = test_data_comp,
    group = "Group",
    min_length = 2,
    max_length = 2,
    top_n = 5,
    detailed = FALSE,
    statistical = FALSE,
    min_char_length_state = 1 # Explicitly set
  )
  
  expect_s3_class(result, "compare_sequences")
  expect_named(result, c("results", "summary", "parameters", "stats"))
  expect_true(is.data.frame(result$summary) || length(result$summary) == 0) # summary can be empty df
  if (nrow(result$summary) > 0) {
    expect_true("discrimination" %in% names(result$summary))
  }
  expect_equal(result$parameters$groups, c("A", "B"))
  expect_equal(result$parameters$min_char_length_state, 1)
})

test_that("compare_sequences statistical mode", {
  result_stat <- compare_sequences(
    data = test_data_comp,
    group = "Group",
    min_length = 2,
    max_length = 2,
    top_n = 5,
    detailed = FALSE,
    statistical = TRUE,
    correction = "none", # Use "none" for predictable p-values
    min_char_length_state = 1
  )
  
  expect_s3_class(result_stat, "compare_sequences")
  if (nrow(result_stat$summary) > 0) {
    expect_true("p_value" %in% names(result_stat$summary))
    expect_true("significant" %in% names(result_stat$summary))
  }
  expect_true(result_stat$parameters$statistical)
})

test_that("compare_sequences handles min_char_length_state", {
  data_single_char <- data.frame(
    T1 = c("A", "B", "C", "A"),
    T2 = c("B", "C", "A", "B"),
    Group = c("G1", "G1", "G2", "G2")
  )
  
  # Default min_char_length_state = 1 (as per function default) should find patterns
  result_mcls1 <- compare_sequences(data_single_char, "Group", min_length = 2, max_length = 2, statistical = FALSE, min_char_length_state = 1)
  # Check if any patterns are found. Subsequences like "A-B", "B-C" should be formed.
  # The exact number of patterns depends on internal logic of table creation and filtering.
  # A simple check is that the summary is not empty if patterns are expected.
  # seqs G1: A-B, B-C. seqs G2: C-A, A-B.
  # Patterns: A-B (G1:1, G2:1), B-C (G1:1, G2:0), C-A (G1:0, G2:1)
  expect_true(nrow(result_mcls1$summary) > 0, info = "Patterns should be found with min_char_length_state = 1")

  # min_char_length_state = 2 should find no patterns from single chars
  result_mcls2 <- compare_sequences(data_single_char, "Group", min_length = 2, max_length = 2, statistical = FALSE, min_char_length_state = 2)
  expect_true(nrow(result_mcls2$summary) == 0, info = "No patterns should be found with min_char_length_state = 2")
})

test_that("compare_sequences handles group as vector", {
  data_only <- test_data_comp[, -which(names(test_data_comp) == "Group")]
  group_vector <- test_data_comp$Group

  result <- compare_sequences(
    data = data_only,
    group = group_vector,
    min_length = 2,
    max_length = 2,
    statistical = FALSE,
    min_char_length_state = 1
  )
  expect_s3_class(result, "compare_sequences")
  expect_equal(result$parameters$groups, c("A", "B"))
})

test_that("compare_sequences error handling for invalid inputs", {
  # Not a dataframe
  expect_error(compare_sequences(data = 1:10, group = "Group", min_char_length_state = 1),
               "'data' must be a data frame")
  # Empty dataframe
  expect_error(compare_sequences(data = data.frame(), group = "Group", min_char_length_state = 1),
               "'data' cannot be empty")
  # Invalid group column
  expect_error(compare_sequences(data = test_data_comp, group = "InvalidGroup", min_char_length_state = 1),
               "Group column 'InvalidGroup' not found in data")
  # Group vector wrong length
  expect_error(compare_sequences(data = test_data_comp, group = c("A", "B"), min_char_length_state = 1),
               "'group' must have the same length as number of rows in 'data'")
  # Not exactly 2 groups
  data_one_group <- test_data_comp[test_data_comp$Group == "A",]
  expect_error(compare_sequences(data = data_one_group, group = "Group", min_char_length_state = 1),
               "Group variable must have exactly 2 levels, found: 1")
})

test_that("Plotting methods for compare_sequences run", {
  # Discrimination mode
  result_disc <- compare_sequences(test_data_comp, "Group", statistical = FALSE, min_char_length_state = 1)
  if (nrow(result_disc$summary) > 0) {
    expect_silent(plot(result_disc))
  } else {
    expect_silent(plot(result_disc))
  }

  # Statistical mode
  result_stat <- compare_sequences(test_data_comp, "Group", statistical = TRUE, min_char_length_state = 1)
   if (nrow(result_stat$summary) > 0) {
    expect_silent(plot(result_stat))
  } else {
    expect_silent(plot(result_stat))
  }
})

test_that("compare_sequences detailed output", {
  result_detailed <- compare_sequences(
    data = test_data_comp,
    group = "Group",
    min_length = 2,
    max_length = 3,
    detailed = TRUE,
    statistical = FALSE,
    min_char_length_state = 1
  )
  expect_type(result_detailed$summary, "list")
  if(length(result_detailed$summary) > 0) {
    expect_true(all(sapply(result_detailed$summary, is.data.frame)))
    # Check if names of the list items follow the pattern "length_X"
    # This depends on whether min_length can equal max_length for detailed view.
    # If min_length=2, max_length=3, names should be "length_2", "length_3"
    if(result_detailed$parameters$min_length < result_detailed$parameters$max_length) {
        expect_true(all(grepl("^length_[0-9]+$", names(result_detailed$summary))))
    } else if (result_detailed$parameters$min_length == result_detailed$parameters$max_length && length(result_detailed$summary) > 0) {
        # If min_length == max_length, still could be named length_X or just be a data frame
        # The current implementation of compare_sequences returns a list of dataframes for detailed=TRUE
        # even if only one length is analyzed.
        expect_true(all(grepl("^length_[0-9]+$", names(result_detailed$summary))))
    }
  }
})

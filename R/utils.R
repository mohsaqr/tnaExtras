# ==============================================================================
# UTILITY FUNCTIONS
# ==============================================================================
# This file contains helper functions used across different modules of the package.
# ==============================================================================

#' Validate Group Input
#'
#' Checks if the group input is valid for analysis functions.
#' Ensures the group vector is a factor, has the correct length, and meets
#' the minimum number of groups required.
#'
#' @param group_vector A vector representing group assignments for each sequence.
#' @param n_sequences The total number of sequences in the dataset.
#' @param min_groups The minimum number of unique groups required (default is 2).
#' @return A factor representing the validated group levels.
#' @examples
#' \dontrun{
#'   # Assuming 'sequences_df' is your data frame of sequences
#'   # and 'group_labels' is a vector of group assignments
#'   n_seqs <- nrow(sequences_df)
#'   validated_groups <- validate_group_input(group_labels, n_seqs)
#' }
#' @noRd
validate_group_input <- function(group_vector, n_sequences, min_groups = 2) {
  if (is.null(group_vector)) {
    stop("Group vector cannot be NULL.", call. = FALSE)
  }
  if (length(group_vector) != n_sequences) {
    stop(paste0("Length of group vector (", length(group_vector),
                ") does not match the number of sequences (", n_sequences, ")."),
         call. = FALSE)
  }

  group_factor <- as.factor(group_vector)
  n_unique_groups <- length(levels(group_factor))

  if (n_unique_groups < min_groups) {
    stop(paste0("Number of unique groups (", n_unique_groups,
                ") is less than the minimum required (", min_groups, ")."),
         call. = FALSE)
  }

  # Check for NA in group assignments after conversion to factor
  if (any(is.na(group_factor))) {
    stop("Group vector contains NA values. Please remove or impute them.", call. = FALSE)
  }

  # Check for empty group levels (e.g., "" if original data had empty strings)
  if (any(levels(group_factor) == "")) {
      stop("Group vector contains empty string ('') as a group level. Please provide meaningful group names.", call. = FALSE)
  }

  return(group_factor)
}

#' Split Data by Group
#'
#' Splits the main data frame into a list of data frames, one for each group.
#'
#' @param data A data frame containing the sequence data (or any data to be split).
#' @param group_factor A factor vector representing group assignments for each row in 'data'.
#' @return A list of data frames, where each element corresponds to a group.
#'         The list is named with the group levels.
#' @examples
#' \dontrun{
#'   # Assuming 'sequences_df' is your data and 'validated_group_factor' is from validate_group_input
#'   grouped_data_list <- split_data_by_group(sequences_df, validated_group_factor)
#' }
#' @noRd
split_data_by_group <- function(data, group_factor) {
  if (!is.data.frame(data)) {
    stop("Input 'data' must be a data frame.", call. = FALSE)
  }
  if (!is.factor(group_factor)) {
    stop("Input 'group_factor' must be a factor.", call. = FALSE)
  }
  if (nrow(data) != length(group_factor)) {
    stop("Number of rows in 'data' does not match the length of 'group_factor'.", call. = FALSE)
  }

  group_levels <- levels(group_factor)

  # Use split() which is efficient for this
  # Ensure the factor levels are preserved in the list names
  # drop = FALSE is not directly applicable to split, but split preserves empty groups if they exist
  # However, validate_group_input should prevent empty or NA groups by this point.
  data_list <- split(data, group_factor)

  # Ensure the list is named correctly according to group_levels
  # split() already names the list elements by the levels of the factor

  # Check if all expected groups are present as list elements
  if (!all(group_levels %in% names(data_list))) {
      stop("Mismatch between group factor levels and split data list names. This should not happen.", call. = FALSE)
  }

  return(data_list)
}

#' Generate Unique Group Pairs
#'
#' Generates all unique pairs of group levels for pairwise comparisons.
#'
#' @param group_levels A character vector of unique group levels.
#' @return A list of character vectors, where each vector contains a pair of group levels.
#'         Returns an empty list if less than two group levels are provided.
#' @examples
#' \dontrun{
#'   gr_levels <- c("A", "B", "C")
#'   group_pairs <- generate_group_pairs(gr_levels)
#'   # group_pairs would be list(c("A", "B"), c("A", "C"), c("B", "C"))
#' }
#' @noRd
generate_group_pairs <- function(group_levels) {
  if (!is.character(group_levels) && !is.factor(group_levels)) {
    stop("group_levels must be a character vector or a factor.", call. = FALSE)
  }

  unique_levels <- unique(as.character(group_levels))

  if (length(unique_levels) < 2) {
    # stop("At least two unique group levels are required to generate pairs.", call. = FALSE)
    # Return empty list if not enough groups to form pairs, consistent with no comparisons
    return(list())
  }

  # Use combn to get all combinations of 2
  pairs_matrix <- utils::combn(unique_levels, 2)

  # Convert matrix columns to a list of vectors
  pairs_list <- lapply(seq_len(ncol(pairs_matrix)), function(i) pairs_matrix[, i])

  return(pairs_list)
}

#' Preprocess Sequence Data (from pattern_analysis.R, generalized)
#'
#' Converts wide format sequence data with a group column into a list of processed
#' sequences for each group. This is an internal helper.
#'
#' @param data Data frame with sequences in wide format.
#' @param group_col_name Character string, the name of the group column.
#' @param data_cols_indices Numeric vector, indices of columns containing sequence data.
#' @param min_length_seq Minimum sequence length to include.
#' @return A list where names are group levels and values are character vectors of sequences for that group.
#'         Also includes an element `all_valid_sequences` (all sequences passing min_length)
#'         and `group_assignments` (factor of group assignments for all_valid_sequences).
#'         Returns NULL if no valid sequences are found or groups are problematic.
#'
#' @noRd
preprocess_sequences_for_analysis <- function(data_df, group_col_name, data_cols_indices, min_length_seq = 1) {
  if (!is.data.frame(data_df)) stop("data_df must be a data frame.", call. = FALSE)
  if (!is.character(group_col_name) || length(group_col_name) != 1) {
    stop("group_col_name must be a single string.", call. = FALSE)
  }
  if (!group_col_name %in% names(data_df)) {
    stop(paste0("Group column '", group_col_name, "' not found in data_df."), call. = FALSE)
  }
  if (!is.numeric(data_cols_indices) || any(data_cols_indices < 1) || any(data_cols_indices > ncol(data_df))) {
    stop("data_cols_indices are invalid.", call. = FALSE)
  }
  if (group_col_name %in% names(data_df)[data_cols_indices]) {
      stop("Group column cannot also be a data column for sequences.", call. = FALSE)
  }

  groups_factor <- as.factor(data_df[[group_col_name]])
  sequence_data_df <- data_df[, data_cols_indices, drop = FALSE]

  # Convert to sequences (hyphen-separated strings)
  all_sequences_char <- apply(sequence_data_df, 1, function(row) {
    valid_states <- row[!is.na(row) & row != "" & nchar(as.character(row)) > 0] # Ensure states have substance
    if (length(valid_states) < min_length_seq) return("") # Return empty if too short
    paste(valid_states, collapse = "-")
  })

  # Filter out sequences that became empty (due to min_length_seq or all NA/empty states)
  valid_sequence_indices <- nchar(all_sequences_char) > 0

  if(sum(valid_sequence_indices) == 0) {
    stop(paste0("No valid sequences found after processing (min_length_seq = ", min_length_seq, ")."), call. = FALSE)
  }

  processed_sequences <- all_sequences_char[valid_sequence_indices]
  processed_groups_factor <- groups_factor[valid_sequence_indices]

  # Drop unused levels from the processed_groups_factor
  processed_groups_factor <- droplevels(processed_groups_factor)

  if (length(levels(processed_groups_factor)) < 1) { # Should be at least 1 if sum(valid_sequence_indices) > 0
      stop("No groups remain after filtering sequences.", call. = FALSE)
  }

  # Split sequences by the processed group factor
  sequences_by_group_list <- split(processed_sequences, processed_groups_factor)

  return(list(
    sequences_by_group = sequences_by_group_list,
    all_valid_sequences = processed_sequences, # All sequences that passed length filter
    group_assignments = processed_groups_factor # Group assignment for each of all_valid_sequences
  ))
}


#' Safe Statistical Test
#'
#' A wrapper for Fisher's Exact Test and Chi-squared Test that handles
#' common errors and selects tests based on expected frequencies.
#'
#' @param contingency_table A 2x2 contingency table.
#' @param test_method Character, statistical test selection ("auto", "fisher", "chi.squared").
#' @param min_expected Numeric, minimum expected count for automatic test selection (default: 5).
#' @return A list containing the p-value, test used, and statistic.
#'         Returns NA for p-value and statistic if the test fails or table is unsuitable.
#'
#' @noRd
safe_statistical_test <- function(contingency_table, test_method = "auto", min_expected = 5) {
  if (!is.matrix(contingency_table) || !all(dim(contingency_table) == c(2, 2))) {
    # stop("Contingency table must be a 2x2 matrix.", call. = FALSE)
    # This can happen if a pattern is not present in one group or not present at all.
    # The calling function should handle this (e.g., by not testing such patterns).
    # For robustness here, return NAs.
    return(list(p_value = NA_real_, test_used = "Invalid Table", statistic = NA_real_))
  }
  if (any(is.na(contingency_table)) || any(contingency_table < 0)) {
    return(list(p_value = NA_real_, test_used = "Invalid Table Values", statistic = NA_real_))
  }

  # If all values are zero, no test is meaningful
  if (sum(contingency_table) == 0) {
    return(list(p_value = 1, test_used = "No Data", statistic = NA_real_)) # p=1 implies no difference
  }

  use_fisher <- FALSE
  test_result_obj <- NULL # Placeholder for the test result object

  if (test_method == "fisher") {
    use_fisher <- TRUE
  } else if (test_method == "auto") {
    # Try chi-squared to get expected values, but catch errors if table is unsuitable
    chisq_attempt <- try(stats::chisq.test(contingency_table), silent = TRUE)
    if (inherits(chisq_attempt, "try-error") || any(chisq_attempt$expected < 1e-9)) {
        # If chi-sq fails (e.g. sum of row/col is 0), or expected values are near zero, default to Fisher
        use_fisher <- TRUE
    } else {
        use_fisher <- any(chisq_attempt$expected < min_expected)
    }
  } # else, test_method == "chi.squared", so use_fisher remains FALSE

  p_val <- NA_real_
  stat_val <- NA_real_
  test_name <- NA_character_

  tryCatch({
    if (use_fisher) {
      test_result_obj <- stats::fisher.test(contingency_table)
      p_val <- test_result_obj$p.value
      test_name <- "Fisher"
      # Fisher's test doesn't have a 'statistic' in the same way chi-sq does; often odds ratio is reported.
      # For consistency, we might leave statistic NA or calculate odds ratio. Here, NA.
    } else {
      # suppressWarnings for Yates' continuity correction if not desired, but default is fine
      test_result_obj <- stats::chisq.test(contingency_table)
      p_val <- test_result_obj$p.value
      stat_val <- test_result_obj$statistic
      test_name <- "Chi-squared"
    }
  }, error = function(e) {
    # This error block will catch issues not handled by preliminary checks
    # e.g., if a row/column sum is zero leading to issues in chi-squared that weren't caught.
    # In such cases, try Fisher as a fallback if not already tried.
    if (!use_fisher && test_method == "auto") {
        warning(paste("Chi-squared test failed, attempting Fisher's Exact Test. Error:", e$message), call.=FALSE)
        tryCatch({
            test_result_obj <- stats::fisher.test(contingency_table)
            p_val <<- test_result_obj$p.value
            test_name <<- "Fisher (fallback)"
        }, error = function(e2) {
            warning(paste("Fallback Fisher's test also failed. Error:", e2$message), call.=FALSE)
            # p_val, stat_val remain NA
            test_name <<- "Test Failed"
        })
    } else {
        warning(paste("Statistical test failed. Error:", e$message), call.=FALSE)
        # p_val, stat_val remain NA
        test_name <<- "Test Failed"
    }
  })

  return(list(p_value = p_val, test_used = test_name, statistic = stat_val))
}

# (More helper functions will be added here in subsequent steps)
# For example:
# - format_pairwise_results()
# - extract_ngrams_from_sequences_list()
# - etc.

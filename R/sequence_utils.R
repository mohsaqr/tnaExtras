# ==============================================================================
# SEQUENCE ANALYSIS UTILITIES
# ==============================================================================
#
# Shared utility functions for sequence pattern analysis.
# These functions eliminate code duplication across sequence_patterns.R
# and sequence_motifs.R
#
# ==============================================================================

# ==============================================================================
# DATA PREPARATION
# ==============================================================================

#' Prepare Sequence Data from Wide Format
#'
#' Converts wide-format data to a list of sequences, handling group_tna objects
#' and filtering invalid entries.
#'
#' @param data Data frame in wide format or group_tna object
#' @param verbose Whether to print progress messages
#' @return List with: sequences, n_sequences, all_states, n_states
#' @keywords internal
prepare_sequence_data <- function(data, verbose = TRUE) {
  

  # Handle group_tna objects
  if (is_group_tna(data)) {
    if (verbose) cat("Converting group_tna object...\n")
    converted <- convert_group_tna(data)
    data <- converted$data
    group_col <- converted$group_col
    if (!is.null(group_col) && group_col %in% names(data)) {
      data <- data[, names(data) != group_col, drop = FALSE]
    }
  }
  
  # Ensure data frame

  if (!is.data.frame(data) && !is.matrix(data)) {
    stop("data must be a data frame or matrix")
  }
  data <- as.data.frame(data)
  
  if (nrow(data) < 2) {
    stop("data must contain at least 2 sequences")
  }
  
  # Convert to sequences
  if (verbose) cat("Processing sequences...\n")
  sequences <- apply(data, 1, function(row) {
    valid <- row[!is.na(row) & row != "" & row != "NA"]
    if (length(valid) == 0) return(NULL)
    as.character(valid)
  })
  sequences <- sequences[!sapply(sequences, is.null)]
  n_sequences <- length(sequences)
  
  if (n_sequences == 0) {
    stop("No valid sequences found in data")
  }
  
  all_states <- unique(unlist(sequences))
  n_states <- length(all_states)
  
  if (verbose) {
    cat("  Sequences:", n_sequences, "\n")
    cat("  Unique states:", n_states, "\n")
  }
  
  list(
    sequences = sequences,
    n_sequences = n_sequences,
    all_states = all_states,
    n_states = n_states,
    original_data = data
  )
}

# ==============================================================================
# FILTERING UTILITIES
# ==============================================================================

#' Filter Patterns by State
#'
#' Filter patterns based on starting state, ending state, or containing state.
#'
#' @param df Data frame with pattern column
#' @param start_state Filter patterns starting with this state
#' @param end_state Filter patterns ending with this state
#' @param contains_state Filter patterns containing this state
#' @param pattern_col Name of the pattern column (default: "pattern")
#' @param separator Pattern separator (default: "->")
#' @return Filtered data frame
#' @keywords internal
filter_by_state <- function(df, start_state = NULL, end_state = NULL, 
                            contains_state = NULL, pattern_col = "pattern",
                            separator = "->") {
  
  if (nrow(df) == 0) return(df)
  
  patterns <- df[[pattern_col]]
  keep <- rep(TRUE, nrow(df))
  
  if (!is.null(start_state)) {
    # Escape special regex characters and match start
    start_escaped <- gsub("([.+*?^${}()\\[\\]|\\\\])", "\\\\\\1", start_state)
    start_pattern <- paste0("^", start_escaped, "(", separator, "|$)")
    keep <- keep & grepl(start_pattern, patterns)
  }
  
  if (!is.null(end_state)) {
    # Escape special regex characters and match end
    end_escaped <- gsub("([.+*?^${}()\\[\\]|\\\\])", "\\\\\\1", end_state)
    end_pattern <- paste0("(", separator, "|^)", end_escaped, "$")
    keep <- keep & grepl(end_pattern, patterns)
  }
  
  if (!is.null(contains_state)) {
    # Escape special regex characters
    contains_escaped <- gsub("([.+*?^${}()\\[\\]|\\\\])", "\\\\\\1", contains_state)
    contains_pattern <- paste0("(^|", separator, ")", contains_escaped, "(", separator, "|$)")
    keep <- keep & grepl(contains_pattern, patterns)
  }
  
  df[keep, , drop = FALSE]
}

#' Filter Patterns by Schema (Type-Level)
#'
#' Filter patterns based on starting schema, ending schema, or containing schema.
#'
#' @param df Data frame with schema column
#' @param start_schema Filter patterns with schema starting with this type
#' @param end_schema Filter patterns with schema ending with this type
#' @param contains_schema Filter patterns with schema containing this type
#' @return Filtered data frame
#' @keywords internal
filter_by_schema <- function(df, start_schema = NULL, end_schema = NULL,
                             contains_schema = NULL) {
  
  if (nrow(df) == 0 || !"schema" %in% names(df)) return(df)
  
  schemas <- df$schema
  keep <- rep(TRUE, nrow(df))
  
  if (!is.null(start_schema)) {
    start_escaped <- gsub("([.+*?^${}()\\[\\]|\\\\])", "\\\\\\1", start_schema)
    start_pattern <- paste0("^", start_escaped, "(->|$)")
    keep <- keep & grepl(start_pattern, schemas)
  }
  
  if (!is.null(end_schema)) {
    end_escaped <- gsub("([.+*?^${}()\\[\\]|\\\\])", "\\\\\\1", end_schema)
    end_pattern <- paste0("(->|^)", end_escaped, "$")
    keep <- keep & grepl(end_pattern, schemas)
  }
  
  if (!is.null(contains_schema)) {
    contains_escaped <- gsub("([.+*?^${}()\\[\\]|\\\\])", "\\\\\\1", contains_schema)
    contains_pattern <- paste0("(^|->)", contains_escaped, "(->|$)")
    keep <- keep & grepl(contains_pattern, schemas)
  }
  
  df[keep, , drop = FALSE]
}

#' Filter by Support and Count Thresholds
#'
#' @param df Data frame with support and count columns
#' @param min_support Minimum support threshold
#' @param min_count Minimum count threshold
#' @return Filtered data frame
#' @keywords internal
filter_by_thresholds <- function(df, min_support = 0.01, min_count = 2) {
  
  if (nrow(df) == 0) return(df)
  
  keep <- rep(TRUE, nrow(df))
  
  if ("support" %in% names(df) && !is.null(min_support)) {
    keep <- keep & (df$support >= min_support)
  }
  
  if ("count" %in% names(df) && !is.null(min_count)) {
    keep <- keep & (df$count >= min_count)
  }
  
  df[keep, , drop = FALSE]
}

# ==============================================================================
# SIGNIFICANCE TESTING
# ==============================================================================

#' Compute Significance Statistics for Patterns
#'
#' Adds statistical measures to a patterns data frame including expected
#' probability, lift, chi-square, z-score, p-values, and significance flags.
#'
#' @param df Data frame with count and sequences_containing columns
#' @param n_sequences Total number of sequences
#' @param n_states Number of unique states (for expected probability)
#' @param correction Multiple testing correction method
#' @param alpha Significance level
#' @param pattern_length_col Column name for pattern length (default: "length")
#' @return Data frame with added significance columns
#' @keywords internal
compute_significance_stats <- function(df, n_sequences, n_states,
                                       correction = "fdr", alpha = 0.05,
                                       pattern_length_col = "length") {
  
  if (nrow(df) == 0) {
    # Add empty columns
    df$expected_prob <- numeric(0)
    df$lift <- numeric(0)
    df$chi_square <- numeric(0)
    df$z_score <- numeric(0)
    df$p_value <- numeric(0)
    df$p_adjusted <- numeric(0)
    df$significant <- logical(0)
    return(df)
  }
  
  # Calculate expected probability based on pattern length
  if (pattern_length_col %in% names(df)) {
    df$expected_prob <- round((1 / n_states) ^ df[[pattern_length_col]], 6)
  } else {
    df$expected_prob <- round(1 / n_states, 6)
  }
  
  # Ensure support is calculated
  if (!"support" %in% names(df) && "sequences_containing" %in% names(df)) {
    df$support <- round(df$sequences_containing / n_sequences, 4)
  }
  
  # Lift
  df$lift <- round(df$support / df$expected_prob, 2)
  
  # Chi-square test
  df$chi_square <- sapply(1:nrow(df), function(i) {
    if (!"sequences_containing" %in% names(df)) return(NA)
    observed <- c(df$sequences_containing[i], 
                  n_sequences - df$sequences_containing[i])
    expected <- c(n_sequences * df$expected_prob[i], 
                  n_sequences * (1 - df$expected_prob[i]))
    if (all(expected > 0)) {
      round(sum((observed - expected)^2 / expected), 2)
    } else {
      NA
    }
  })
  
  # Z-score
  df$z_score <- sapply(1:nrow(df), function(i) {
    prop <- df$support[i]
    exp_prob <- df$expected_prob[i]
    se <- sqrt(exp_prob * (1 - exp_prob) / n_sequences)
    if (se > 0) round((prop - exp_prob) / se, 2) else NA
  })
  
  # Binomial test p-values
  df$p_value <- sapply(1:nrow(df), function(i) {
    if (!"sequences_containing" %in% names(df)) return(NA)
    tryCatch({
      stats::binom.test(df$sequences_containing[i], n_sequences,
                        df$expected_prob[i], alternative = "greater")$p.value
    }, error = function(e) NA)
  })
  
  # Multiple testing correction
  df$p_adjusted <- stats::p.adjust(df$p_value, method = correction)
  
  # Significance flag
  df$significant <- !is.na(df$p_adjusted) & df$p_adjusted < alpha
  
  df
}

# ==============================================================================
# OUTPUT FORMATTING
# ==============================================================================
#' Create Empty Patterns Data Frame
#'
#' Creates a data frame with the standard column structure for patterns output.
#'
#' @param include_schema Whether to include schema column
#' @return Empty data frame with proper columns
#' @keywords internal
create_empty_patterns_df <- function(include_schema = FALSE) {
  
  df <- data.frame(
    pattern = character(0),
    length = integer(0),
    count = integer(0),
    sequences_containing = integer(0),
    support = numeric(0),
    proportion = numeric(0),
    expected_prob = numeric(0),
    lift = numeric(0),
    chi_square = numeric(0),
    z_score = numeric(0),
    p_value = numeric(0),
    p_adjusted = numeric(0),
    significant = logical(0),
    stringsAsFactors = FALSE
  )
  
  if (include_schema) {
    df$schema <- character(0)
    # Reorder to put schema after pattern
    cols <- c("pattern", "schema", setdiff(names(df), c("pattern", "schema")))
    df <- df[, cols, drop = FALSE]
  }
  
  df
}

#' Standardize Output Column Order
#'
#' Ensures consistent column ordering in output data frames.
#'
#' @param df Data frame to reorder
#' @param include_schema Whether schema column should be included
#' @return Data frame with standardized column order
#' @keywords internal
standardize_columns <- function(df, include_schema = FALSE) {
  
  if (nrow(df) == 0) return(df)
  
  # Standard column order
  if (include_schema) {
    preferred_order <- c("pattern", "schema", "length", "count", 
                         "sequences_containing", "support", "proportion",
                         "expected_prob", "lift", "chi_square", "z_score",
                         "p_value", "p_adjusted", "significant")
  } else {
    preferred_order <- c("pattern", "length", "count", 
                         "sequences_containing", "support", "proportion",
                         "expected_prob", "lift", "chi_square", "z_score",
                         "p_value", "p_adjusted", "significant")
  }
  
  # Get columns that exist in df
  existing_cols <- intersect(preferred_order, names(df))
  # Add any extra columns at the end
  extra_cols <- setdiff(names(df), preferred_order)
  
  df[, c(existing_cols, extra_cols), drop = FALSE]
}

# ==============================================================================
# PRINT/SUMMARY UTILITIES
# ==============================================================================

#' Print Patterns Summary
#'
#' @param patterns_df Data frame of patterns
#' @param title Title for the section
#' @param n_show Number of patterns to show
#' @keywords internal
print_patterns_summary <- function(patterns_df, title = "Patterns", n_show = 15) {
  
  cat(title, "\n")
  cat(paste(rep("-", 70), collapse = ""), "\n")
  
  if (is.null(patterns_df) || nrow(patterns_df) == 0) {
    cat("  No patterns found.\n\n")
    return(invisible(NULL))
  }
  
  # Select display columns
  display_cols <- c("pattern", "schema", "count", "sequences_containing", 
                    "support", "lift", "significant")
  display_cols <- display_cols[display_cols %in% names(patterns_df)]
  
  top_patterns <- head(patterns_df[order(patterns_df$count, decreasing = TRUE), ], n_show)
  print(top_patterns[, display_cols, drop = FALSE], row.names = FALSE)
  
  if (nrow(patterns_df) > n_show) {
    cat(sprintf("\n... and %d more patterns\n", nrow(patterns_df) - n_show))
  }
  
  cat("\n  support = sequences_containing / total_sequences\n")
  
  if ("significant" %in% names(patterns_df)) {
    n_sig <- sum(patterns_df$significant, na.rm = TRUE)
    cat(sprintf("  Significant patterns: %d / %d (%.1f%%)\n",
                n_sig, nrow(patterns_df), 100 * n_sig / nrow(patterns_df)))
  }
  cat("\n")
}

cat("Sequence utilities loaded.\n")


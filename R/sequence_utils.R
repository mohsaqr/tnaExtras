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
# DATA PREPARATION & VALIDATION
# ==============================================================================

#' Prepare Sequence Input
#'
#' Converts wide-format data to a list of sequences, handling group_tna objects
#' and filtering invalid entries.
#'
#' @param data Data frame in wide format or group_tna object
#' @param verbose Whether to print progress messages
#' @return List with: sequences, n_sequences, all_states, n_states
#' @keywords internal
prepare_sequence_input <- function(data, verbose = TRUE) {
  
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

#' Validate Sequence Analysis Parameters
#' 
#' Validates common parameters for sequence analysis functions.
#' 
#' @param params Named list of parameters to validate
#' @keywords internal
validate_sequence_params <- function(params) {
  
  if (!is.null(params$min_length) && params$min_length < 1) {
    stop("min_length must be at least 1")
  }
  
  if (!is.null(params$max_length) && !is.null(params$min_length)) {
    if (params$max_length < params$min_length) {
      stop("max_length must be >= min_length")
    }
  }
  
  if (!is.null(params$min_support)) {
    if (params$min_support < 0 || params$min_support > 1) {
      stop("min_support must be between 0 and 1")
    }
  }
  
  if (!is.null(params$min_count) && params$min_count < 1) {
    stop("min_count must be at least 1")
  }
  
  if (!is.null(params$correction)) {
    valid_corrections <- c("bonferroni", "holm", "hochberg", "hommel", "BH", "BY", "fdr", "none")
    if (!params$correction %in% valid_corrections) {
      stop("correction must be one of: ", paste(valid_corrections, collapse = ", "))
    }
  }
}

# ==============================================================================
# VECTORIZED EXTRACTION HELPERS (Optimization)
# ==============================================================================

#' Extract N-grams Vectorized
#' 
#' Uses vector shifting to extract n-grams without inner loops.
#' 
#' @param sequences List of character vectors
#' @param n N-gram length
#' @return Character vector of all n-grams (with repetition)
#' @keywords internal
extract_ngrams_fast <- function(sequences, n) {
  if (n < 1) return(character(0))
  
  # Vectorized extraction using lapply and paste
  # For n=2: paste(x[1:L-1], x[2:L])
  ngrams_list <- lapply(sequences, function(seq) {
    len <- length(seq)
    if (len < n) return(NULL)
    
    # Create shifted versions
    args <- lapply(0:(n-1), function(i) seq[(1 + i):(len - n + 1 + i)])
    do.call(paste, c(args, sep = "->"))
  })
  
  unlist(ngrams_list)
}

#' Extract Gapped Patterns Vectorized
#' 
#' Uses vector shifting to extract patterns with a fixed gap.
#' Pattern: A -> *gap* -> B
#' 
#' @param sequences List of character vectors
#' @param gap Gap size (number of wildcards)
#' @return Character vector of all gapped patterns
#' @keywords internal
extract_gapped_fast <- function(sequences, gap) {
  # Pattern length is gap + 2 (A, gap, B)
  # Lag is gap + 1
  lag <- gap + 1
  sep_str <- paste0("->", paste(rep("*", gap), collapse = "->"), "->")
  
  patterns_list <- lapply(sequences, function(seq) {
    len <- length(seq)
    if (len <= lag) return(NULL)
    
    # A is at 1:(len-lag), B is at (1+lag):len
    paste(seq[1:(len - lag)], seq[(1 + lag):len], sep = sep_str)
  })
  
  unlist(patterns_list)
}

# ==============================================================================
# FILTERING UTILITIES
# ==============================================================================

#' Filter Patterns by Text Matches
#'
#' Generalized filter for patterns based on start/end/contains criteria.
#' Replaces specific filter_by_state and filter_by_schema functions.
#'
#' @param df Data frame with pattern column
#' @param start Filter patterns starting with this text
#' @param end Filter patterns ending with this text
#' @param contains Filter patterns containing this text
#' @param text_col Name of the column to filter (default: "pattern")
#' @param separator Pattern separator (default: "->")
#' @return Filtered data frame
#' @keywords internal
filter_patterns_by_text <- function(df, start = NULL, end = NULL, 
                                    contains = NULL, text_col = "pattern",
                                    separator = "->") {
  
  if (nrow(df) == 0 || !text_col %in% names(df)) return(df)
  
  patterns <- df[[text_col]]
  keep <- rep(TRUE, nrow(df))
  
  # Escape special regex characters helper
  esc <- function(x) gsub("([.+*?^${}()\\[\\]|\\\\])", "\\\\\\1", x)
  
  if (!is.null(start)) {
    start_pattern <- paste0("^", esc(start), "(", separator, "|$)")
    keep <- keep & grepl(start_pattern, patterns)
  }
  
  if (!is.null(end)) {
    end_pattern <- paste0("(", separator, "|^)", esc(end), "$")
    keep <- keep & grepl(end_pattern, patterns)
  }
  
  if (!is.null(contains)) {
    contains_pattern <- paste0("(^|", separator, ")", esc(contains), "(", separator, "|$)")
    keep <- keep & grepl(contains_pattern, patterns)
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
# AGGREGATION & STATISTICS
# ==============================================================================

#' Aggregate Instances into Statistics Table
#' 
#' Converts a list of raw instance strings into a statistical data frame.
#' Handles counting, support calculation, and significance testing.
#' 
#' @param instance_list Vector of instance strings (e.g., "A->B->C")
#' @param instance_seqs List tracking which sequences contain each instance (optional)
#' @param n_sequences Total number of sequences
#' @param n_states Total number of states/types (for expected probability)
#' @param schema Associated schema string (optional)
#' @param test_significance Whether to compute stats
#' @param correction Correction method
#' @param alpha Significance level
#' @return Data frame with pattern statistics
#' @keywords internal
aggregate_instances <- function(instance_list, instance_seqs = NULL, n_sequences, n_states,
                                schema = NULL, test_significance = TRUE,
                                correction = "fdr", alpha = 0.05) {
  
  if (length(instance_list) == 0) {
    return(create_empty_patterns_df(include_schema = !is.null(schema)))
  }
  
  # If instance_seqs is NULL, we assume we need to calculate it or infer support differently
  # For optimized vectorization, we might not pass instance_seqs explicitly to save memory
  # but we need it for accurate support (sequences containing).
  
  # Optimized aggregation using data.table logic with base R
  inst_df <- data.frame(instance = instance_list, stringsAsFactors = FALSE)
  
  # Filter out empty or invalid instances
  instance_list <- instance_list[nchar(instance_list) > 0 & !is.na(instance_list)]

  if (length(instance_list) == 0) {
    return(create_empty_patterns_df(include_schema = !is.null(schema)))
  }

  # Count frequency
  counts <- table(instance_list)

  # Handle empty table (shouldn't happen after filtering, but safety check)
  if (length(counts) == 0) {
    return(create_empty_patterns_df(include_schema = !is.null(schema)))
  }

  # Support (sequences containing)
  if (is.null(instance_seqs)) {
    # Fallback if we don't have sequence IDs tracked: assume count ~ support (bad)
    # Or, optimized functions should pass a mapping of instance -> sequence_id
    stop("instance_seqs required for accurate support calculation") 
  } else {
    if (is.list(instance_seqs)) {
        # Map from list
        seqs_per_pattern <- sapply(names(counts), function(p) length(unique(instance_seqs[[p]])))
    } else {
        # instance_seqs is a vector aligned with instance_list
        df_map <- data.frame(pattern = instance_list, seq_id = instance_seqs, stringsAsFactors = FALSE)
        seqs_per_pattern <- tapply(df_map$seq_id, df_map$pattern, function(x) length(unique(x)))
        # Reorder to match counts
        seqs_per_pattern <- seqs_per_pattern[names(counts)]
    }
  }
  
  # Final safety check
  n_patterns <- length(counts)
  if (n_patterns == 0) {
    return(create_empty_patterns_df(include_schema = !is.null(schema)))
  }

  # Ensure seqs_per_pattern has correct length
  if (length(seqs_per_pattern) != n_patterns) {
    seqs_per_pattern <- rep(1L, n_patterns)
  }

  # Create data frame with explicit length checks
  pattern_vec <- names(counts)
  if (is.null(pattern_vec)) pattern_vec <- character(length(counts))
  length_vec <- sapply(strsplit(as.character(pattern_vec), "->"), length)
  count_vec <- as.integer(counts)
  seqs_vec <- as.integer(seqs_per_pattern)

  # Ensure all vectors have the same length
  vec_lengths <- c(length(pattern_vec), length(length_vec), length(count_vec), length(seqs_vec))
  if (length(unique(vec_lengths)) > 1) {
    warning(sprintf("Vector length mismatch: %s", paste(vec_lengths, collapse = ", ")))
    # Use the minimum length
    min_len <- min(vec_lengths)
    pattern_vec <- pattern_vec[1:min_len]
    length_vec <- length_vec[1:min_len]
    count_vec <- count_vec[1:min_len]
    seqs_vec <- seqs_vec[1:min_len]
  }

  df <- data.frame(
    pattern = pattern_vec,
    length = length_vec,
    count = count_vec,
    sequences_containing = seqs_vec,
    stringsAsFactors = FALSE
  )
  
  if (!is.null(schema)) {
    df$schema <- rep(schema, nrow(df))
  }
  
  # Basic metrics
  df$support <- round(df$sequences_containing / n_sequences, 4)
  df$proportion <- round(df$count / sum(df$count), 4)
  
  # Significance stats
  if (test_significance) {
    df <- compute_significance_stats(df, n_sequences, n_states, correction, alpha)
  } else {
    # Add empty placeholders
    cols <- c("expected_prob", "lift", "chi_square", "z_score", "p_value", "p_adjusted", "significant")
    df[cols] <- NA
  }
  
  # Standardize column order
  standardize_columns(df, include_schema = !is.null(schema))
}

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
  
  if (nrow(df) == 0) return(df)
  
  # Calculate expected probability
  # For state patterns: (1/n_states)^length
  # For mixed patterns: 1 / n_unique_patterns (heuristic)
  if (pattern_length_col %in% names(df) && !is.null(n_states)) {
    df$expected_prob <- round((1 / n_states) ^ df[[pattern_length_col]], 6)
  } else {
    # Fallback for when length isn't available or n_states is undefined
    n_unique <- nrow(df)
    df$expected_prob <- round(1 / max(n_unique, 1), 6)
  }
  
  # Ensure support is calculated
  if (!"support" %in% names(df) && "sequences_containing" %in% names(df)) {
    df$support <- round(df$sequences_containing / n_sequences, 4)
  }
  
  # Lift
  df$lift <- round(df$support / df$expected_prob, 2)
  
  # Chi-square test (only use when expected counts >= 5)
  df$chi_square <- sapply(seq_len(nrow(df)), function(i) {
    if (!"sequences_containing" %in% names(df)) return(NA)
    exp_count <- n_sequences * df$expected_prob[i]
    if (exp_count >= 5 && (n_sequences - exp_count) >= 5) {
      obs <- c(df$sequences_containing[i], n_sequences - df$sequences_containing[i])
      exp <- c(exp_count, n_sequences - exp_count)
      round(sum((obs - exp)^2 / exp), 2)
    } else {
      NA  # Chi-square not appropriate for low expected counts
    }
  })
  
  # Z-score
  df$z_score <- sapply(seq_len(nrow(df)), function(i) {
    prop <- df$support[i]
    exp_prob <- df$expected_prob[i]
    se <- sqrt(exp_prob * (1 - exp_prob) / n_sequences)
    if (se > 0) round((prop - exp_prob) / se, 2) else NA
  })
  
  # Binomial test p-values (use when appropriate)
  df$p_value <- sapply(seq_len(nrow(df)), function(i) {
    if (!"sequences_containing" %in% names(df)) return(NA)
    tryCatch({
      # Only use binomial test if expected count >= 5 and <= n-5
      exp_count <- n_sequences * df$expected_prob[i]
      if (exp_count >= 5 && exp_count <= n_sequences - 5) {
        stats::binom.test(df$sequences_containing[i], n_sequences,
                          df$expected_prob[i], alternative = "greater")$p.value
      } else {
        # For extreme probabilities, use normal approximation or Fisher exact test
        # For now, use z-test p-value as approximation
        prop <- df$support[i]
        exp_prob <- df$expected_prob[i]
        se <- sqrt(exp_prob * (1 - exp_prob) / n_sequences)
        if (se > 0) {
          z <- (prop - exp_prob) / se
          1 - stats::pnorm(z)  # One-tailed p-value
        } else {
          NA
        }
      }
    }, error = function(e) NA)
  })
  
  # Multiple testing correction
  df$p_adjusted <- stats::p.adjust(df$p_value, method = correction)
  df$significant <- !is.na(df$p_adjusted) & df$p_adjusted < alpha
  
  df
}

# ==============================================================================
# OUTPUT FORMATTING & DISPLAY
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
  }
  
  standardize_columns(df, include_schema)
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
  
  preferred_order <- c("pattern", "schema", "length", "count", 
                       "sequences_containing", "support", "proportion",
                       "expected_prob", "lift", "chi_square", "z_score",
                       "p_value", "p_adjusted", "significant")
  
  if (!include_schema) preferred_order <- setdiff(preferred_order, "schema")
  
  # Get columns that exist in df
  existing_cols <- intersect(preferred_order, names(df))
  extra_cols <- setdiff(names(df), preferred_order)
  
  df[, c(existing_cols, extra_cols), drop = FALSE]
}

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

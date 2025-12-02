# ==============================================================================
# SEQUENCE PATTERN EXPLORATION TOOLKIT
# ==============================================================================
# 
# A comprehensive toolkit for exploring sequence frequency, pattern discovery,
# and statistical significance testing in wide-format sequence data.
#
# Main function: explore_sequence_patterns()
#
# ==============================================================================

#' Explore Sequence Patterns with Statistical Significance
#'
#' Analyzes sequences in wide format to discover frequent patterns, compute support
#' frequencies, and test for statistical significance. This function is designed
#' for exploratory analysis of sequence data without group comparisons.
#'
#' @param data Data frame with sequences in wide format (rows are sequences, 
#'   columns are time points). Can also be a group_tna object.
#' @param min_length Minimum pattern length (n-gram size) to analyze (default: 1)
#' @param max_length Maximum pattern length (n-gram size) to analyze (default: 5)
#' @param min_support Minimum support threshold to include a pattern (default: 0.01, i.e., 1%)
#' @param min_count Minimum count threshold to include a pattern (default: 2)
#' @param test_significance Whether to perform statistical significance tests (default: TRUE)
#' @param correction Method for multiple testing correction (default: "BH" for Benjamini-Hochberg).
#'   Options: "bonferroni", "holm", "hochberg", "hommel", "BH", "BY", "fdr", "none"
#' @param alpha Significance level for statistical tests (default: 0.05)
#' @param expected_prob Expected probability under null hypothesis for binomial test.
#'   If NULL (default), uses uniform probability based on state space.
#' @param verbose Whether to print progress messages (default: TRUE)
#'
#' @return An object of class "sequence_pattern_analysis" containing:
#'   \item{patterns}{Data frame with pattern frequencies, support, and significance}
#'   \item{full_sequences}{Data frame with complete sequence frequencies}
#'   \item{state_frequencies}{Data frame with individual state frequencies}
#'   \item{summary}{List with summary statistics}
#'   \item{parameters}{List of analysis parameters}
#'
#' @examples
#' \dontrun{
#' # Load sequence data in wide format
#' data <- read.csv("sequences.csv")
#' 
#' # Basic pattern exploration
#' results <- explore_sequence_patterns(data)
#' print(results)
#' 
#' # Custom analysis with stricter thresholds
#' results <- explore_sequence_patterns(
#'   data, 
#'   min_length = 2, 
#'   max_length = 4,
#'   min_support = 0.05,
#'   correction = "bonferroni"
#' )
#' 
#' # View significant patterns only
#' significant_patterns(results)
#' 
#' # Plot most frequent sequences
#' plot(results, type = "sequences", top_n = 15)
#' }
#'
#' @export
explore_sequence_patterns <- function(data,
                                     min_length = 1,
                                     max_length = 5,
                                     min_support = 0.01,
                                     min_count = 2,
                                     test_significance = TRUE,
                                     correction = "BH",
                                     alpha = 0.05,
                                     expected_prob = NULL,
                                     verbose = TRUE) {
  
  # =====================================================================
  # INPUT VALIDATION
  # =====================================================================
  
  # Check if input is a group_tna object
  if (is_group_tna(data)) {
    if (verbose) cat("Detected group_tna object, converting to tnaExtras format...\n")
    converted <- convert_group_tna(data)
    data <- converted$data
    # Remove group column for single-group analysis
    group_col <- converted$group_col
    if (!is.null(group_col) && group_col %in% names(data)) {
      data <- data[, names(data) != group_col, drop = FALSE]
    }
  }
  
  if (!is.data.frame(data) && !is.matrix(data)) {
    stop("data must be a data frame or matrix")
  }
  
  data <- as.data.frame(data)
  
  if (nrow(data) < 2) {
    stop("data must contain at least 2 sequences")
  }
  
  if (ncol(data) < 1) {
    stop("data must contain at least 1 time point")
  }
  
  # Validate parameters
  if (min_length < 1) {
    stop("min_length must be at least 1")
  }
  
  if (max_length < min_length) {
    stop("max_length must be greater than or equal to min_length")
  }
  
  if (min_support < 0 || min_support > 1) {
    stop("min_support must be between 0 and 1")
  }
  
  if (min_count < 1) {
    stop("min_count must be at least 1")
  }
  
  valid_corrections <- c("bonferroni", "holm", "hochberg", "hommel", "BH", "BY", "fdr", "none")
  if (!correction %in% valid_corrections) {
    stop("correction must be one of: ", paste(valid_corrections, collapse = ", "))
  }
  
  # =====================================================================
  # SEQUENCE EXTRACTION
  # =====================================================================
  
  if (verbose) cat("Extracting sequences...\n")
  
  n_sequences <- nrow(data)
  
  # Convert wide format to sequence strings
  sequences <- apply(data, 1, function(row) {
    valid_states <- row[!is.na(row) & row != "" & row != "NA"]
    if (length(valid_states) == 0) return(NA_character_)
    paste(valid_states, collapse = "-")
  })
  
  # Remove empty sequences
  valid_indices <- !is.na(sequences) & nchar(sequences) > 0
  sequences <- sequences[valid_indices]
  n_valid <- length(sequences)
  
  if (n_valid == 0) {
    stop("No valid sequences found in the data")
  }
  
  if (verbose) {
    cat("  Total sequences:", n_sequences, "\n")
    cat("  Valid sequences:", n_valid, "\n")
  }
  
  # Extract all unique states
  all_states <- unique(unlist(strsplit(sequences, "-")))
  all_states <- all_states[!is.na(all_states) & all_states != ""]
  n_states <- length(all_states)
  
  if (verbose) {
    cat("  Unique states:", n_states, "\n")
    cat("  States:", paste(head(all_states, 10), collapse = ", "), 
        if (n_states > 10) paste("...", "(", n_states - 10, "more)") else "", "\n")
  }
  
  # =====================================================================
  # PATTERN EXTRACTION
  # =====================================================================
  
  if (verbose) cat("\nExtracting patterns...\n")
  
  # Extract n-grams of all specified lengths
  all_patterns <- list()
  
  for (n in min_length:max_length) {
    if (verbose) cat("  Extracting", n, "-grams...\n")
    
    ngrams <- character(0)
    for (seq in sequences) {
      seq_ngrams <- extract_ngrams_for_pattern(seq, n)
      ngrams <- c(ngrams, seq_ngrams)
    }
    
    if (length(ngrams) > 0) {
      all_patterns[[as.character(n)]] <- ngrams
    }
  }
  
  # =====================================================================
  # FREQUENCY COMPUTATION
  # =====================================================================
  
  if (verbose) cat("\nComputing frequencies...\n")
  
  # 1. State frequencies (unigrams)
  state_freq_table <- table(unlist(strsplit(sequences, "-")))
  total_states <- sum(state_freq_table)
  
  state_frequencies <- data.frame(
    state = names(state_freq_table),
    count = as.numeric(state_freq_table),
    frequency = as.numeric(state_freq_table) / total_states,
    stringsAsFactors = FALSE
  )
  state_frequencies <- state_frequencies[order(state_frequencies$count, decreasing = TRUE), ]
  rownames(state_frequencies) <- NULL
  
  # 2. Full sequence frequencies
  full_seq_table <- table(sequences)
  
  full_sequences <- data.frame(
    sequence = names(full_seq_table),
    count = as.numeric(full_seq_table),
    support = as.numeric(full_seq_table) / n_valid,
    stringsAsFactors = FALSE
  )
  full_sequences$length <- sapply(strsplit(full_sequences$sequence, "-"), length)
  full_sequences <- full_sequences[order(full_sequences$count, decreasing = TRUE), ]
  rownames(full_sequences) <- NULL
  
  # 3. Pattern (n-gram) frequencies
  patterns_df <- compute_pattern_frequencies(all_patterns, sequences, n_valid, 
                                            min_support, min_count, verbose)
  
  # =====================================================================
  # STATISTICAL SIGNIFICANCE TESTING
  # =====================================================================
  
  if (test_significance && nrow(patterns_df) > 0) {
    if (verbose) cat("\nPerforming statistical significance tests...\n")
    
    patterns_df <- compute_pattern_significance(
      patterns_df, 
      n_valid, 
      n_states,
      expected_prob,
      correction,
      alpha,
      verbose
    )
    
    # Also test full sequences
    if (nrow(full_sequences) > 0) {
      full_sequences <- compute_sequence_significance(
        full_sequences,
        n_valid,
        n_states,
        correction,
        alpha
      )
    }
  }
  
  # =====================================================================
  # SUMMARY STATISTICS
  # =====================================================================
  
  summary_stats <- list(
    n_sequences_total = n_sequences,
    n_sequences_valid = n_valid,
    n_unique_states = n_states,
    n_unique_patterns = nrow(patterns_df),
    n_unique_full_sequences = nrow(full_sequences),
    n_significant_patterns = if (test_significance && "significant" %in% names(patterns_df)) 
      sum(patterns_df$significant) else NA,
    mean_sequence_length = mean(sapply(strsplit(sequences, "-"), length)),
    median_sequence_length = median(sapply(strsplit(sequences, "-"), length)),
    min_sequence_length = min(sapply(strsplit(sequences, "-"), length)),
    max_sequence_length = max(sapply(strsplit(sequences, "-"), length)),
    most_common_state = state_frequencies$state[1],
    most_common_pattern = if (nrow(patterns_df) > 0) patterns_df$pattern[1] else NA_character_,
    state_entropy = compute_entropy(state_frequencies$frequency)
  )
  
  # =====================================================================
  # CREATE RESULT OBJECT
  # =====================================================================
  
  result <- list(
    patterns = patterns_df,
    full_sequences = full_sequences,
    state_frequencies = state_frequencies,
    summary = summary_stats,
    parameters = list(
      min_length = min_length,
      max_length = max_length,
      min_support = min_support,
      min_count = min_count,
      test_significance = test_significance,
      correction = correction,
      alpha = alpha,
      all_states = all_states
    )
  )
  
  class(result) <- "sequence_pattern_analysis"
  
  if (verbose) {
    cat("\n=== Analysis Complete ===\n")
    cat("Unique patterns found:", nrow(patterns_df), "\n")
    cat("Unique full sequences:", nrow(full_sequences), "\n")
    if (test_significance && "significant" %in% names(patterns_df)) {
      cat("Significant patterns:", sum(patterns_df$significant), "\n")
    }
  }
  
  return(result)
}

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

#' Extract n-grams from a sequence string
#' @param sequence Sequence string with states separated by "-"
#' @param n N-gram length
#' @return Character vector of n-grams
extract_ngrams_for_pattern <- function(sequence, n) {
  if (is.na(sequence) || nchar(sequence) == 0) return(character(0))
  
  parts <- unlist(strsplit(sequence, "-"))
  parts <- parts[!is.na(parts) & parts != ""]
  
  if (length(parts) < n) return(character(0))
  
  ngrams <- character(length(parts) - n + 1)
  for (i in 1:(length(parts) - n + 1)) {
    ngrams[i] <- paste(parts[i:(i + n - 1)], collapse = "-")
  }
  
  return(ngrams)
}

#' Compute pattern frequencies
#' @param all_patterns List of patterns by n-gram length
#' @param sequences Vector of sequence strings
#' @param n_valid Number of valid sequences
#' @param min_support Minimum support threshold
#' @param min_count Minimum count threshold
#' @param verbose Whether to print progress
#' @return Data frame with pattern frequencies
compute_pattern_frequencies <- function(all_patterns, sequences, n_valid, 
                                       min_support, min_count, verbose) {
  
  # Combine all patterns
  all_ngrams <- unlist(all_patterns, use.names = FALSE)
  
  if (length(all_ngrams) == 0) {
    return(data.frame(
      pattern = character(0),
      length = integer(0),
      count = integer(0),
      support = numeric(0),
      sequences_containing = integer(0),
      stringsAsFactors = FALSE
    ))
  }
  
  # Count pattern occurrences
  pattern_counts <- table(all_ngrams)
  
  # Create data frame
  patterns_df <- data.frame(
    pattern = names(pattern_counts),
    count = as.numeric(pattern_counts),
    stringsAsFactors = FALSE
  )
  
  # Calculate pattern lengths
  patterns_df$length <- sapply(strsplit(patterns_df$pattern, "-"), length)
  
  # Calculate support (proportion of sequences containing the pattern)
  patterns_df$sequences_containing <- sapply(patterns_df$pattern, function(p) {
    sum(grepl(p, sequences, fixed = TRUE))
  })
  
  patterns_df$support <- patterns_df$sequences_containing / n_valid
  
  # Filter by thresholds
  min_count_threshold <- max(min_count, ceiling(min_support * n_valid))
  patterns_df <- patterns_df[patterns_df$sequences_containing >= min_count & 
                            patterns_df$support >= min_support, ]
  
  # Sort by count (descending)
  patterns_df <- patterns_df[order(patterns_df$count, decreasing = TRUE), ]
  rownames(patterns_df) <- NULL
  
  if (verbose) {
    cat("  Patterns meeting thresholds:", nrow(patterns_df), "\n")
  }
  
  return(patterns_df)
}

#' Compute pattern significance
#' @param patterns_df Data frame with pattern frequencies
#' @param n_sequences Number of sequences
#' @param n_states Number of unique states
#' @param expected_prob Expected probability under null
#' @param correction Multiple testing correction method
#' @param alpha Significance level
#' @param verbose Whether to print progress
#' @return Data frame with significance results
compute_pattern_significance <- function(patterns_df, n_sequences, n_states,
                                        expected_prob, correction, alpha, verbose) {
  
  if (nrow(patterns_df) == 0) return(patterns_df)
  
  # Calculate p-values using binomial test
  # H0: Pattern occurs with expected probability
  # H1: Pattern occurs more frequently than expected by chance
  
  p_values <- numeric(nrow(patterns_df))
  expected_probs <- numeric(nrow(patterns_df))
  observed_expected <- numeric(nrow(patterns_df))
  
  for (i in 1:nrow(patterns_df)) {
    pattern_length <- patterns_df$length[i]
    sequences_with_pattern <- patterns_df$sequences_containing[i]
    
    # Expected probability under null hypothesis
    # If not specified, use uniform distribution assumption
    if (is.null(expected_prob)) {
      # Expected probability that a sequence contains this specific pattern
      # Under uniform distribution: (1/n_states)^pattern_length * possible_positions
      # Simplified: we use a conservative estimate based on state space
      exp_prob <- (1 / n_states) ^ pattern_length
      # Adjust for multiple possible positions in a sequence
      avg_seq_length <- mean(sapply(strsplit(patterns_df$pattern, "-"), length))
      exp_prob <- min(exp_prob * max(1, avg_seq_length - pattern_length + 1), 1)
    } else {
      exp_prob <- expected_prob
    }
    
    expected_probs[i] <- exp_prob
    observed_expected[i] <- patterns_df$support[i] / exp_prob
    
    # Binomial test: is the observed count significantly higher than expected?
    if (sequences_with_pattern > 0) {
      # Use one-sided binomial test for "greater than expected"
      test_result <- stats::binom.test(
        x = sequences_with_pattern,
        n = n_sequences,
        p = exp_prob,
        alternative = "greater"
      )
      p_values[i] <- test_result$p.value
    } else {
      p_values[i] <- 1
    }
  }
  
  # Apply multiple testing correction
  if (correction != "none") {
    adjusted_p <- stats::p.adjust(p_values, method = correction)
  } else {
    adjusted_p <- p_values
  }
  
  # Add results to data frame
  patterns_df$expected_prob <- expected_probs
  patterns_df$observed_expected_ratio <- round(observed_expected, 2)
  patterns_df$p_value <- p_values
  patterns_df$p_adjusted <- adjusted_p
  patterns_df$significant <- adjusted_p < alpha
  
  # Sort by significance then by count
  patterns_df <- patterns_df[order(patterns_df$p_adjusted, -patterns_df$count), ]
  rownames(patterns_df) <- NULL
  
  if (verbose) {
    cat("  Significant patterns (p <", alpha, "):", sum(patterns_df$significant), "\n")
  }
  
  return(patterns_df)
}

#' Compute sequence significance
#' @param sequences_df Data frame with sequence frequencies
#' @param n_sequences Number of sequences
#' @param n_states Number of unique states
#' @param correction Multiple testing correction method
#' @param alpha Significance level
#' @return Data frame with significance results
compute_sequence_significance <- function(sequences_df, n_sequences, n_states,
                                         correction, alpha) {
  
  if (nrow(sequences_df) == 0) return(sequences_df)
  
  p_values <- numeric(nrow(sequences_df))
  
  for (i in 1:nrow(sequences_df)) {
    seq_length <- sequences_df$length[i]
    count <- sequences_df$count[i]
    
    # Expected probability under uniform distribution
    exp_prob <- (1 / n_states) ^ seq_length
    
    # Binomial test
    if (count > 0) {
      test_result <- stats::binom.test(
        x = count,
        n = n_sequences,
        p = exp_prob,
        alternative = "greater"
      )
      p_values[i] <- test_result$p.value
    } else {
      p_values[i] <- 1
    }
  }
  
  # Apply multiple testing correction
  if (correction != "none") {
    adjusted_p <- stats::p.adjust(p_values, method = correction)
  } else {
    adjusted_p <- p_values
  }
  
  sequences_df$p_value <- p_values
  sequences_df$p_adjusted <- adjusted_p
  sequences_df$significant <- adjusted_p < alpha
  
  return(sequences_df)
}

#' Compute Shannon entropy
#' @param probs Vector of probabilities
#' @return Entropy value
compute_entropy <- function(probs) {
  probs <- probs[probs > 0]
  if (length(probs) == 0) return(0)
  -sum(probs * log2(probs))
}

# ==============================================================================
# S3 METHODS
# ==============================================================================

#' Print method for sequence_pattern_analysis objects
#' @param x sequence_pattern_analysis object
#' @param ... additional arguments
#' @export
print.sequence_pattern_analysis <- function(x, ...) {
  cat("Sequence Pattern Analysis Results\n")
  cat("=================================\n\n")
  
  cat("Data Summary:\n")
  cat("  Total sequences:", x$summary$n_sequences_total, "\n")
  cat("  Valid sequences:", x$summary$n_sequences_valid, "\n")
  cat("  Unique states:", x$summary$n_unique_states, "\n")
  cat("  Sequence length: mean =", round(x$summary$mean_sequence_length, 1),
      ", range = [", x$summary$min_sequence_length, ",", x$summary$max_sequence_length, "]\n")
  cat("  State entropy:", round(x$summary$state_entropy, 3), "\n\n")
  
  cat("Pattern Analysis:\n")
  cat("  Unique patterns found:", x$summary$n_unique_patterns, "\n")
  cat("  Unique full sequences:", x$summary$n_unique_full_sequences, "\n")
  if (!is.na(x$summary$n_significant_patterns)) {
    cat("  Significant patterns:", x$summary$n_significant_patterns, "\n")
  }
  cat("\n")
  
  cat("Analysis Parameters:\n")
  cat("  Pattern length range:", x$parameters$min_length, "-", x$parameters$max_length, "\n")
  cat("  Min support:", x$parameters$min_support, "\n")
  cat("  Min count:", x$parameters$min_count, "\n")
  if (x$parameters$test_significance) {
    cat("  Correction method:", x$parameters$correction, "\n")
    cat("  Significance level:", x$parameters$alpha, "\n")
  }
  cat("\n")
  
  # Show top patterns
  if (nrow(x$patterns) > 0) {
    cat("Top 10 Most Frequent Patterns:\n")
    top_patterns <- head(x$patterns[order(x$patterns$count, decreasing = TRUE), ], 10)
    
    for (i in 1:nrow(top_patterns)) {
      sig_marker <- if ("significant" %in% names(top_patterns) && top_patterns$significant[i]) " ***" else ""
      cat(sprintf("  %2d. %s (count=%d, support=%.3f)%s\n", 
                 i, top_patterns$pattern[i], top_patterns$count[i], 
                 top_patterns$support[i], sig_marker))
    }
    
    if ("significant" %in% names(x$patterns)) {
      cat("\n*** = statistically significant (p <", x$parameters$alpha, "after", 
          x$parameters$correction, "correction)\n")
    }
  }
  
  cat("\nUse summary() for detailed statistics, or access components directly:\n")
  cat("  $patterns, $full_sequences, $state_frequencies\n")
  
  invisible(x)
}

#' Summary method for sequence_pattern_analysis objects
#' @param object sequence_pattern_analysis object
#' @param ... additional arguments
#' @export
summary.sequence_pattern_analysis <- function(object, ...) {
  cat("Sequence Pattern Analysis - Detailed Summary\n")
  cat("=============================================\n\n")
  
  # Data summary
  cat("DATA OVERVIEW\n")
  cat(paste(rep("-", 40), collapse = ""), "\n")
  cat("Sequences analyzed:", object$summary$n_sequences_valid, "of", 
      object$summary$n_sequences_total, "total\n")
  cat("Unique states:", object$summary$n_unique_states, "\n")
  cat("Sequence length:\n")
  cat("  Mean:", round(object$summary$mean_sequence_length, 2), "\n")
  cat("  Median:", object$summary$median_sequence_length, "\n")
  cat("  Range:", object$summary$min_sequence_length, "-", 
      object$summary$max_sequence_length, "\n\n")
  
  # State frequencies
  cat("STATE FREQUENCIES (Top 10)\n")
  cat(paste(rep("-", 40), collapse = ""), "\n")
  top_states <- head(object$state_frequencies, 10)
  for (i in 1:nrow(top_states)) {
    cat(sprintf("  %-20s count=%5d  freq=%.3f\n",
               top_states$state[i], top_states$count[i], top_states$frequency[i]))
  }
  cat("\n")
  
  # Pattern summary by length
  if (nrow(object$patterns) > 0) {
    cat("PATTERNS BY LENGTH\n")
    cat(paste(rep("-", 40), collapse = ""), "\n")
    
    for (len in unique(object$patterns$length)) {
      len_patterns <- object$patterns[object$patterns$length == len, ]
      n_patterns <- nrow(len_patterns)
      n_sig <- if ("significant" %in% names(len_patterns)) sum(len_patterns$significant) else NA
      
      cat(sprintf("  Length %d: %d patterns", len, n_patterns))
      if (!is.na(n_sig)) cat(sprintf(" (%d significant)", n_sig))
      cat("\n")
      
      # Top pattern of this length
      top <- len_patterns[which.max(len_patterns$count), ]
      cat(sprintf("    Top: %s (count=%d, support=%.3f)\n", 
                 top$pattern[1], top$count[1], top$support[1]))
    }
    cat("\n")
  }
  
  # Significant patterns
  if ("significant" %in% names(object$patterns) && any(object$patterns$significant)) {
    cat("SIGNIFICANT PATTERNS (Top 10)\n")
    cat(paste(rep("-", 40), collapse = ""), "\n")
    
    sig_patterns <- object$patterns[object$patterns$significant, ]
    sig_patterns <- sig_patterns[order(sig_patterns$p_adjusted), ]
    top_sig <- head(sig_patterns, 10)
    
    for (i in 1:nrow(top_sig)) {
      cat(sprintf("  %s\n", top_sig$pattern[i]))
      cat(sprintf("    count=%d, support=%.3f, O/E=%.1f, p_adj=%.2e\n",
                 top_sig$count[i], top_sig$support[i], 
                 top_sig$observed_expected_ratio[i], top_sig$p_adjusted[i]))
    }
    cat("\n")
  }
  
  # Most frequent full sequences
  if (nrow(object$full_sequences) > 0) {
    cat("MOST FREQUENT COMPLETE SEQUENCES (Top 5)\n")
    cat(paste(rep("-", 40), collapse = ""), "\n")
    top_seqs <- head(object$full_sequences, 5)
    for (i in 1:nrow(top_seqs)) {
      # Truncate long sequences for display
      seq_display <- top_seqs$sequence[i]
      if (nchar(seq_display) > 50) {
        seq_display <- paste0(substr(seq_display, 1, 47), "...")
      }
      cat(sprintf("  %d. %s\n", i, seq_display))
      cat(sprintf("     count=%d, support=%.3f, length=%d\n",
                 top_seqs$count[i], top_seqs$support[i], top_seqs$length[i]))
    }
  }
  
  invisible(object)
}

#' Extract significant patterns
#' @param x sequence_pattern_analysis object
#' @param sort_by Column to sort by (default: "p_adjusted")
#' @return Data frame of significant patterns
#' @export
significant_patterns <- function(x, sort_by = "p_adjusted") {
  if (!inherits(x, "sequence_pattern_analysis")) {
    stop("x must be a sequence_pattern_analysis object")
  }
  
  if (!"significant" %in% names(x$patterns)) {
    stop("No significance testing was performed. Run explore_sequence_patterns with test_significance = TRUE")
  }
  
  sig <- x$patterns[x$patterns$significant, ]
  
  if (nrow(sig) == 0) {
    message("No significant patterns found")
    return(sig)
  }
  
  if (sort_by %in% names(sig)) {
    sig <- sig[order(sig[[sort_by]]), ]
  }
  
  return(sig)
}

#' Plot sequence pattern analysis results
#' 
#' Visualize patterns, states, or full sequences from the analysis.
#' 
#' @param x sequence_pattern_analysis object
#' @param type Type of plot: "patterns", "states", or "sequences"
#' @param top_n Number of top items to plot (default: 20)
#' @param show_support For sequences plot, show support values (default: TRUE)
#' @param col Bar color(s). For patterns, can be a single color or will auto-color
#'   by significance if available.
#' @param ... Additional arguments passed to barplot
#' @export
plot.sequence_pattern_analysis <- function(x, type = "patterns", top_n = 20, 
                                          show_support = TRUE, col = NULL, ...) {
  
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par))
  
  if (type == "patterns") {
    if (nrow(x$patterns) == 0) {
      message("No patterns to plot")
      return(invisible(NULL))
    }
    
    data <- head(x$patterns[order(x$patterns$count, decreasing = TRUE), ], top_n)
    
    # Set colors
    if (is.null(col)) {
      if ("significant" %in% names(data)) {
        col <- ifelse(data$significant, "steelblue", "gray70")
      } else {
        col <- "steelblue"
      }
    }
    
    graphics::par(mar = c(10, 4, 4, 2))
    bp <- graphics::barplot(
      data$count,
      names.arg = data$pattern,
      las = 2,
      main = paste("Top", min(top_n, nrow(data)), "Patterns by Frequency"),
      ylab = "Count",
      col = col,
      ...
    )
    
    if ("significant" %in% names(data) && any(data$significant) && is.null(list(...)$col)) {
      graphics::legend("topright", 
                      legend = c("Significant", "Not significant"),
                      fill = c("steelblue", "gray70"),
                      bty = "n")
    }
    
  } else if (type == "states") {
    data <- head(x$state_frequencies, top_n)
    
    if (is.null(col)) col <- "coral"
    
    graphics::par(mar = c(8, 4, 4, 2))
    graphics::barplot(
      data$count,
      names.arg = data$state,
      las = 2,
      main = paste("Top", min(top_n, nrow(data)), "States by Frequency"),
      ylab = "Count",
      col = col,
      ...
    )
    
  } else if (type == "sequences") {
    if (nrow(x$full_sequences) == 0) {
      message("No sequences to plot")
      return(invisible(NULL))
    }
    
    data <- head(x$full_sequences[order(x$full_sequences$count, decreasing = TRUE), ], top_n)
    
    # Truncate long sequence labels
    seq_labels <- sapply(data$sequence, function(s) {
      if (nchar(s) > 40) {
        paste0(substr(s, 1, 37), "...")
      } else {
        s
      }
    })
    
    # Set colors based on significance if available
    if (is.null(col)) {
      if ("significant" %in% names(data)) {
        col <- ifelse(data$significant, "darkgreen", "gray70")
      } else {
        col <- "darkgreen"
      }
    }
    
    graphics::par(mar = c(12, 4, 4, 2))
    bp <- graphics::barplot(
      data$count,
      names.arg = seq_labels,
      las = 2,
      main = paste("Top", min(top_n, nrow(data)), "Most Frequent Sequences"),
      ylab = "Count",
      col = col,
      ...
    )
    
    # Add support labels on bars if requested
    if (show_support) {
      graphics::text(
        x = bp,
        y = data$count,
        labels = sprintf("%.1f%%", data$support * 100),
        pos = 3,
        cex = 0.7,
        col = "black"
      )
    }
    
    if ("significant" %in% names(data) && any(data$significant) && is.null(list(...)$col)) {
      graphics::legend("topright", 
                      legend = c("Significant", "Not significant"),
                      fill = c("darkgreen", "gray70"),
                      bty = "n")
    }
    
  } else {
    stop("type must be 'patterns', 'states', or 'sequences'")
  }
  
  invisible(NULL)
}

# ==============================================================================
# ADDITIONAL UTILITY FUNCTIONS
# ==============================================================================

#' Filter patterns by criteria
#' @param x sequence_pattern_analysis object
#' @param min_support Minimum support threshold
#' @param min_count Minimum count threshold
#' @param pattern_length Specific pattern length (or vector of lengths)
#' @param significant_only Return only significant patterns
#' @return Filtered data frame of patterns
#' @export
filter_patterns <- function(x, min_support = NULL, min_count = NULL,
                           pattern_length = NULL, significant_only = FALSE) {
  
  if (!inherits(x, "sequence_pattern_analysis")) {
    stop("x must be a sequence_pattern_analysis object")
  }
  
  patterns <- x$patterns
  
  if (!is.null(min_support)) {
    patterns <- patterns[patterns$support >= min_support, ]
  }
  
  if (!is.null(min_count)) {
    patterns <- patterns[patterns$count >= min_count, ]
  }
  
  if (!is.null(pattern_length)) {
    patterns <- patterns[patterns$length %in% pattern_length, ]
  }
  
  if (significant_only) {
    if (!"significant" %in% names(patterns)) {
      stop("No significance testing was performed")
    }
    patterns <- patterns[patterns$significant, ]
  }
  
  return(patterns)
}

#' Get most frequent sequences
#' 
#' Extract the most frequent complete sequences from the analysis.
#' 
#' @param x sequence_pattern_analysis object
#' @param top_n Number of top sequences to return (default: 10)
#' @param min_support Minimum support threshold (optional)
#' @param significant_only Return only significant sequences (default: FALSE)
#' @return Data frame of most frequent sequences
#' @export
top_sequences <- function(x, top_n = 10, min_support = NULL, significant_only = FALSE) {
  
  if (!inherits(x, "sequence_pattern_analysis")) {
    stop("x must be a sequence_pattern_analysis object")
  }
  
  seqs <- x$full_sequences
  
  if (!is.null(min_support)) {
    seqs <- seqs[seqs$support >= min_support, ]
  }
  
  if (significant_only) {
    if (!"significant" %in% names(seqs)) {
      stop("No significance testing was performed")
    }
    seqs <- seqs[seqs$significant, ]
  }
  
  seqs <- seqs[order(seqs$count, decreasing = TRUE), ]
  
  return(head(seqs, top_n))
}

cat("Sequence pattern exploration toolkit loaded.\n")
cat("Main function: explore_sequence_patterns()\n")
cat("Helper functions: significant_patterns(), filter_patterns(), top_sequences()\n")

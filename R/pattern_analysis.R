# ==============================================================================
# PATTERN ANALYSIS TOOLKIT
# ==============================================================================
# 
# A comprehensive toolkit for analyzing pattern differences between groups
# in sequential data. Implements multiple measures including support-based, 
# lift-based, confidence-based, and effect size measures.
#
# Functions:
# - analyze_patterns(): Main comprehensive analysis function
# - compute_support_measures(): Support-based difference measures
# - compute_lift_measures(): Lift-based difference measures  
# - compute_confidence_measures(): Confidence-based difference measures
# - compute_effect_sizes(): Effect size measures
# - Helper functions for data processing and n-gram extraction
#
# ==============================================================================

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

#' Extract n-grams from a sequence
#' 
#' @param sequence Character string with elements separated by "-"
#' @param n Integer specifying the n-gram length
#' @param min_char_length_state Minimum character length for a state to be considered valid (default: 2)
#' @return Character vector of n-grams
extract_ngrams <- function(sequence, n, min_char_length_state = 2) {
  if (nchar(sequence) == 0) return(character(0))
  
  seq_parts <- unlist(strsplit(sequence, "-"))
  # Filter based on min_char_length_state
  if (min_char_length_state > 0) {
    seq_parts <- seq_parts[nchar(seq_parts) >= min_char_length_state]
  }
  
  if (length(seq_parts) < n) return(character(0))
  
  ngrams <- character(0)
  for (i in 1:(length(seq_parts) - n + 1)) {
    ngram <- paste(seq_parts[i:(i + n - 1)], collapse = "-")
    ngrams <- c(ngrams, ngram)
  }
  
  return(ngrams)
}

#' Prepare data for analysis
#' 
#' @param data Data frame with sequences in wide format
#' @param group_col Column name or index containing group information
#' @param min_length Minimum sequence length to include
#' @return List with processed sequences and groups
prepare_sequence_data <- function(data, group_col = "Group", min_length = 2) {
  # Handle group column selection programmatically
  if (is.character(group_col)) {
    if (!group_col %in% names(data)) {
      stop("Group column '", group_col, "' not found in data")
    }
    group_idx <- which(names(data) == group_col)
  } else if (is.numeric(group_col)) {
    group_idx <- group_col
    group_col <- names(data)[group_idx]
  } else {
    stop("group_col must be either a column name or column index")
  }
  
  # Extract group information
  groups <- data[[group_col]]
  
  # Extract sequence data (all columns except group)
  sequence_data <- data[, -group_idx, drop = FALSE]
  
  # Convert to sequences
  all_sequences <- apply(sequence_data, 1, function(row) {
    valid_states <- row[!is.na(row) & row != ""]
    if (length(valid_states) < min_length) return("")
    paste(valid_states, collapse = "-")
  })
  
  # Remove empty sequences
  valid_indices <- nchar(all_sequences) > 0
  all_sequences <- all_sequences[valid_indices]
  groups <- groups[valid_indices]
  
  if (length(all_sequences) == 0) {
    stop("No valid sequences found after processing")
  }
  
  if (length(unique(groups)) != 2) {
    stop("Exactly two groups must be present in the data")
  }
  
  return(list(
    sequences = all_sequences,
    groups = groups,
    group_names = unique(groups)
  ))
}

# ==============================================================================
# SUPPORT MEASURES
# ==============================================================================

#' Compute Support-based Differences
#'
#' Support measures how frequently a pattern appears in each group.
#' Support is defined as the proportion of sequences containing the pattern.
#'
#' @param patterns Character vector of patterns to analyze
#' @param group_A_seqs Character vector of sequences for group A
#' @param group_B_seqs Character vector of sequences for group B
#' @return Data frame with support measures
compute_support_measures <- function(patterns, group_A_seqs, group_B_seqs) {
  # Input validation
  if (!is.character(patterns) || length(patterns) == 0) {
    stop("patterns must be a non-empty character vector")
  }
  if (!is.character(group_A_seqs) || !is.character(group_B_seqs)) {
    stop("group sequences must be character vectors")
  }
  if (length(group_A_seqs) == 0 || length(group_B_seqs) == 0) {
    stop("group sequences cannot be empty")
  }
  
  results <- data.frame(
    pattern = patterns,
    support_A = numeric(length(patterns)),
    support_B = numeric(length(patterns)),
    support_diff = numeric(length(patterns)),
    support_ratio = numeric(length(patterns)),
    relative_support = numeric(length(patterns)),
    stringsAsFactors = FALSE
  )
  
  n_A <- length(group_A_seqs)
  n_B <- length(group_B_seqs)
  
  for (i in seq_along(patterns)) {
    pattern <- patterns[i]
    
    # Count occurrences in each group
    count_A <- sum(grepl(pattern, group_A_seqs, fixed = TRUE))
    count_B <- sum(grepl(pattern, group_B_seqs, fixed = TRUE))
    
    # Support = frequency / total sequences
    support_A <- count_A / n_A
    support_B <- count_B / n_B
    
    results$support_A[i] <- support_A
    results$support_B[i] <- support_B
    results$support_diff[i] <- support_A - support_B
    
    # Handle support ratio with small constant to avoid infinity
    results$support_ratio[i] <- (support_A + 0.001) / (support_B + 0.001)
    
    # Relative support (normalized difference)
    total_support <- support_A + support_B
    results$relative_support[i] <- ifelse(total_support == 0, 0, 
                                         (support_A - support_B) / total_support)
  }
  
  # Sort by absolute support difference
  results <- results[order(abs(results$support_diff), decreasing = TRUE), ]
  rownames(results) <- NULL
  
  return(results)
}

# ==============================================================================
# LIFT MEASURES
# ==============================================================================

#' Compute Lift-based Differences
#'
#' Lift measures how much more likely a pattern is in one group versus the other
#' compared to what would be expected by chance.
#'
#' @param patterns Character vector of patterns to analyze
#' @param group_A_seqs Character vector of sequences for group A
#' @param group_B_seqs Character vector of sequences for group B
#' @return Data frame with lift measures
compute_lift_measures <- function(patterns, group_A_seqs, group_B_seqs) {
  # Input validation
  if (!is.character(patterns) || length(patterns) == 0) {
    stop("patterns must be a non-empty character vector")
  }
  if (!is.character(group_A_seqs) || !is.character(group_B_seqs)) {
    stop("group sequences must be character vectors")
  }
  if (length(group_A_seqs) == 0 || length(group_B_seqs) == 0) {
    stop("group sequences cannot be empty")
  }
  
  results <- data.frame(
    pattern = patterns,
    observed_A = numeric(length(patterns)),
    observed_B = numeric(length(patterns)),
    expected_A = numeric(length(patterns)),
    expected_B = numeric(length(patterns)),
    lift_A = numeric(length(patterns)),
    lift_B = numeric(length(patterns)),
    lift_ratio = numeric(length(patterns)),
    lift_difference = numeric(length(patterns)),
    stringsAsFactors = FALSE
  )
  
  n_A <- length(group_A_seqs)
  n_B <- length(group_B_seqs)
  n_total <- n_A + n_B
  
  for (i in seq_along(patterns)) {
    pattern <- patterns[i]
    
    # Count occurrences
    count_A <- sum(grepl(pattern, group_A_seqs, fixed = TRUE))
    count_B <- sum(grepl(pattern, group_B_seqs, fixed = TRUE))
    count_total <- count_A + count_B
    
    # Skip patterns with no occurrences
    if (count_total == 0) {
      results[i, 2:9] <- 0
      next
    }
    
    # Expected counts under independence (null hypothesis)
    expected_A <- count_total * (n_A / n_total)
    expected_B <- count_total * (n_B / n_total)
    
    # Observed proportions P(Pattern | Group)
    prop_A <- count_A / n_A # P(Pattern | Group A)
    prop_B <- count_B / n_B # P(Pattern | Group B)
    
    # Overall probability of the pattern P(Pattern)
    prob_pattern <- count_total / n_total
    
    # Lift = P(Pattern | Group) / P(Pattern)
    # Add small constant to prob_pattern to avoid division by zero if pattern never occurs
    # (though loop skips if count_total == 0, this is for safety if n_total is huge and count_total is small)
    lift_A <- prop_A / (prob_pattern + 1e-10)
    lift_B <- prop_B / (prob_pattern + 1e-10)
    
    results$observed_A[i] <- count_A
    results$observed_B[i] <- count_B
    results$expected_A[i] <- n_A * prob_pattern # Expected count in Group A if pattern is independent of group
    results$expected_B[i] <- n_B * prob_pattern # Expected count in Group B
    results$lift_A[i] <- lift_A
    results$lift_B[i] <- lift_B
    
    # Lift ratio (capped to avoid extreme values)
    # Ensure lift_B is not zero for ratio calculation; add small constant
    lift_ratio_raw <- lift_A / (lift_B + 1e-10)
    results$lift_ratio[i] <- pmin(pmax(lift_ratio_raw, 0.01), 100)
    
    # Lift difference
    results$lift_difference[i] <- lift_A - lift_B
  }
  
  # Sort by lift difference
  results <- results[order(abs(results$lift_difference), decreasing = TRUE), ]
  rownames(results) <- NULL
  
  return(results)
}

# ==============================================================================
# CONFIDENCE MEASURES
# ==============================================================================

#' Compute Confidence-based Differences
#'
#' Confidence measures the conditional probability of patterns within groups.
#'
#' @param patterns Character vector of patterns to analyze
#' @param group_A_seqs Character vector of sequences for group A
#' @param group_B_seqs Character vector of sequences for group B
#' @return Data frame with confidence measures
compute_confidence_measures <- function(patterns, group_A_seqs, group_B_seqs) {
  # Input validation
  if (!is.character(patterns) || length(patterns) == 0) {
    stop("patterns must be a non-empty character vector")
  }
  if (!is.character(group_A_seqs) || !is.character(group_B_seqs)) {
    stop("group sequences must be character vectors")
  }
  if (length(group_A_seqs) == 0 || length(group_B_seqs) == 0) {
    stop("group sequences cannot be empty")
  }
  
  results <- data.frame(
    pattern = patterns,
    confidence_A = numeric(length(patterns)),
    confidence_B = numeric(length(patterns)),
    confidence_diff = numeric(length(patterns)),
    confidence_ratio = numeric(length(patterns)),
    balanced_confidence = numeric(length(patterns)),
    stringsAsFactors = FALSE
  )
  
  n_A <- length(group_A_seqs)
  n_B <- length(group_B_seqs)
  
  for (i in seq_along(patterns)) {
    pattern <- patterns[i]
    
    # Count occurrences
    count_A <- sum(grepl(pattern, group_A_seqs, fixed = TRUE))
    count_B <- sum(grepl(pattern, group_B_seqs, fixed = TRUE))
    count_total <- count_A + count_B
    
    # Skip patterns with no occurrences
    if (count_total == 0) {
      results[i, 2:6] <- 0
      next
    }
    
    # Confidence = P(pattern | group)
    conf_A <- count_A / count_total
    conf_B <- count_B / count_total
    
    results$confidence_A[i] <- conf_A
    results$confidence_B[i] <- conf_B
    results$confidence_diff[i] <- conf_A - conf_B
    
    # Confidence ratio with small constant
    results$confidence_ratio[i] <- (conf_A + 0.001) / (conf_B + 0.001)
    
    # Balanced confidence:
    # This measure compares the "gain" in confidence for Group A (P(GroupA|Pattern) - P(GroupA))
    # with the "gain" for Group B (P(GroupB|Pattern) - P(GroupB)).
    # A positive value suggests the pattern is more discriminative for Group A compared to its baseline,
    # relative to how discriminative it is for Group B compared to its baseline.
    # P(GroupA) is the expected confidence if pattern and group are independent.
    expected_A_prior <- n_A / (n_A + n_B) # P(GroupA)
    expected_B_prior <- n_B / (n_A + n_B) # P(GroupB)
    
    results$balanced_confidence[i] <- (conf_A - expected_A_prior) - (conf_B - expected_B_prior)
  }
  
  # Sort by absolute confidence difference
  results <- results[order(abs(results$confidence_diff), decreasing = TRUE), ]
  rownames(results) <- NULL
  
  return(results)
}

# ==============================================================================
# EFFECT SIZE MEASURES
# ==============================================================================

#' Compute Effect Size Measures
#'
#' Various effect size measures for pattern differences.
#'
#' @param patterns Character vector of patterns to analyze
#' @param group_A_seqs Character vector of sequences for group A
#' @param group_B_seqs Character vector of sequences for group B
#' @return Data frame with effect size measures
compute_effect_sizes <- function(patterns, group_A_seqs, group_B_seqs) {
  # Input validation
  if (!is.character(patterns) || length(patterns) == 0) {
    stop("patterns must be a non-empty character vector")
  }
  if (!is.character(group_A_seqs) || !is.character(group_B_seqs)) {
    stop("group sequences must be character vectors")
  }
  if (length(group_A_seqs) == 0 || length(group_B_seqs) == 0) {
    stop("group sequences cannot be empty")
  }
  
  results <- data.frame(
    pattern = patterns,
    cohens_h = numeric(length(patterns)),
    cohens_d = numeric(length(patterns)),
    cramers_v = numeric(length(patterns)),
    phi_coefficient = numeric(length(patterns)),
    standardized_diff = numeric(length(patterns)),
    stringsAsFactors = FALSE
  )
  
  n_A <- length(group_A_seqs)
  n_B <- length(group_B_seqs)
  
  for (i in seq_along(patterns)) {
    pattern <- patterns[i]
    
    # Count occurrences
    count_A <- sum(grepl(pattern, group_A_seqs, fixed = TRUE))
    count_B <- sum(grepl(pattern, group_B_seqs, fixed = TRUE))
    
    # Proportions
    p_A <- count_A / n_A
    p_B <- count_B / n_B
    
    # Cohen's h (effect size for proportions)
    h <- 2 * (asin(sqrt(p_A)) - asin(sqrt(p_B)))
    results$cohens_h[i] <- h
    
    # Cohen's d (standardized mean difference)
    pooled_p <- (count_A + count_B) / (n_A + n_B)
    pooled_var <- pooled_p * (1 - pooled_p)
    d <- ifelse(pooled_var == 0, 0, (p_A - p_B) / sqrt(pooled_var))
    results$cohens_d[i] <- d
    
    # Contingency table for other measures
    present_A <- count_A
    absent_A <- n_A - count_A
    present_B <- count_B
    absent_B <- n_B - count_B
    
    # Chi-square components
    n_total <- n_A + n_B
    
    chi_sq <- tryCatch({
      if (n_total == 0) {
        0
      } else {
        expected_11 <- (present_A + present_B) * (present_A + absent_A) / n_total
        expected_12 <- (present_A + present_B) * (present_B + absent_B) / n_total
        expected_21 <- (absent_A + absent_B) * (present_A + absent_A) / n_total
        expected_22 <- (absent_A + absent_B) * (present_B + absent_B) / n_total
        
        # Avoid division by zero
        if (any(c(expected_11, expected_12, expected_21, expected_22) < 1e-10)) {
          0
        } else {
          ((present_A - expected_11)^2 / expected_11) +
          ((present_B - expected_12)^2 / expected_12) +
          ((absent_A - expected_21)^2 / expected_21) +
          ((absent_B - expected_22)^2 / expected_22)
        }
      }
    }, error = function(e) 0)
    
    # Cramer's V and Phi coefficient
    results$cramers_v[i] <- sqrt(chi_sq / n_total)
    results$phi_coefficient[i] <- sqrt(chi_sq / n_total)
    
    # Standardized difference (simple but robust)
    pooled_se <- sqrt(pooled_p * (1 - pooled_p) * (1/n_A + 1/n_B))
    results$standardized_diff[i] <- ifelse(pooled_se == 0, 0, (p_A - p_B) / pooled_se)
  }
  
  # Sort by Cohen's h
  results <- results[order(abs(results$cohens_h), decreasing = TRUE), ]
  rownames(results) <- NULL
  
  return(results)
}

# ==============================================================================
# MAIN ANALYSIS FUNCTION
# ==============================================================================

#' Comprehensive Pattern Analysis
#'
#' Compute multiple difference measures for patterns in sequential data.
#'
#' @param data Data frame with sequences in wide format
#' @param group_col Column name or index containing group information (default: "Group")
#' @param min_length Minimum pattern length to analyze (default: 2)
#' @param max_length Maximum pattern length to analyze (default: 5)
#' @param min_frequency Minimum frequency required to include a pattern (default: 2)
#' @param measures Character vector of measures to compute (default: all)
#' @param min_char_length_state Minimum character length for a state to be considered valid during n-gram extraction (default: 1, meaning all non-empty states are considered).
#' @return List containing analysis results for each measure
analyze_patterns <- function(data, group_col = "Group", min_length = 2, max_length = 5,
                             min_frequency = 2, 
                             measures = c("support", "lift", "confidence", "effect_size"),
                             min_char_length_state = 1) { # New parameter
  
  # Input validation
  if (!is.data.frame(data)) {
    stop("data must be a data frame")
  }
  if (nrow(data) == 0) {
    stop("data cannot be empty")
  }
  
  # Prepare data
  cat("Preparing sequence data...\n")
  prepared_data <- prepare_sequence_data(data, group_col, min_length)
  
  all_sequences <- prepared_data$sequences
  groups <- prepared_data$groups
  group_names <- prepared_data$group_names
  
  # Split by group
  group_A_seqs <- all_sequences[groups == group_names[1]]
  group_B_seqs <- all_sequences[groups == group_names[2]]
  
  cat(sprintf("Group %s: %d sequences\n", group_names[1], length(group_A_seqs)))
  cat(sprintf("Group %s: %d sequences\n", group_names[2], length(group_B_seqs)))
  
  # Extract all unique patterns
  cat("Extracting patterns...\n")
  all_patterns <- character(0)
  for (length in min_length:max_length) {
    for (seq in all_sequences) {
      if (nchar(seq) > 0) {
        patterns <- extract_ngrams(seq, length, min_char_length_state = min_char_length_state) # Pass it
        all_patterns <- c(all_patterns, patterns)
      }
    }
  }
  
  # Filter by frequency
  pattern_counts <- table(all_patterns)
  frequent_patterns <- names(pattern_counts[pattern_counts >= min_frequency])
  
  if (length(frequent_patterns) == 0) {
    stop("No patterns meet the minimum frequency threshold")
  }
  
  cat(sprintf("Found %d patterns meeting frequency threshold\n", length(frequent_patterns)))
  
  # Compute measures
  results <- list()
  
  valid_measures <- c("support", "lift", "confidence", "effect_size")
  if (!all(measures %in% valid_measures)) {
    stop("measures must be one or more of: ", paste(valid_measures, collapse = ", "))
  }
  
  if ("support" %in% measures) {
    cat("Computing support measures...\n")
    results$support <- compute_support_measures(frequent_patterns, group_A_seqs, group_B_seqs)
  }
  
  if ("lift" %in% measures) {
    cat("Computing lift measures...\n")
    results$lift <- compute_lift_measures(frequent_patterns, group_A_seqs, group_B_seqs)
  }
  
  if ("confidence" %in% measures) {
    cat("Computing confidence measures...\n")
    results$confidence <- compute_confidence_measures(frequent_patterns, group_A_seqs, group_B_seqs)
  }
  
  if ("effect_size" %in% measures) {
    cat("Computing effect size measures...\n")
    results$effect_size <- compute_effect_sizes(frequent_patterns, group_A_seqs, group_B_seqs)
  }
  
  # Add metadata
  results$metadata <- list(
    group_names = group_names,
    n_sequences = c(length(group_A_seqs), length(group_B_seqs)),
    n_patterns = length(frequent_patterns),
    min_length = min_length,
    max_length = max_length,
    min_frequency = min_frequency
  )
  
  cat("Analysis complete!\n")
  
  class(results) <- "pattern_analysis"
  return(results)
}

# ==============================================================================
# PRINT AND SUMMARY METHODS
# ==============================================================================

#' Print method for pattern_analysis objects
#' @param x pattern_analysis object
#' @param ... additional arguments
print.pattern_analysis <- function(x, ...) {
  cat("Pattern Analysis Results\n")
  cat("========================\n\n")
  
  if (!is.null(x$metadata)) {
    cat("Groups:", paste(x$metadata$group_names, collapse = " vs "), "\n")
    cat("Sequences:", paste(x$metadata$n_sequences, collapse = " vs "), "\n")
    cat("Patterns analyzed:", x$metadata$n_patterns, "\n")
    cat("Pattern length range:", x$metadata$min_length, "to", x$metadata$max_length, "\n\n")
  }
  
  cat("Available measures:\n")
  for (measure in names(x)) {
    if (measure != "metadata") {
      n_patterns <- nrow(x[[measure]])
      cat(sprintf("  %s: %d patterns\n", measure, n_patterns))
    }
  }
  
  cat("\nUse summary(result, 'measure_name') to view detailed results.\n")
}

#' Summary method for pattern_analysis objects
#' @param object pattern_analysis object
#' @param measure which measure to summarize
#' @param top_n number of top patterns to show
#' @param ... additional arguments
summary.pattern_analysis <- function(object, measure = NULL, top_n = 10, ...) {
  available_measures <- names(object)[names(object) != "metadata"]
  
  if (is.null(measure)) {
    measure <- available_measures[1]
  }
  
  if (!measure %in% available_measures) {
    stop("Measure '", measure, "' not found. Available: ", 
         paste(available_measures, collapse = ", "))
  }
  
  cat("Pattern Analysis Summary -", toupper(measure), "Measures\n")
  cat(paste(rep("=", 50), collapse = ""), "\n\n")
  
  data <- object[[measure]]
  n_show <- min(top_n, nrow(data))
  
  cat("Top", n_show, "patterns:\n")
  print(head(data, n_show))
  
  if (nrow(data) > n_show) {
    cat(sprintf("\n... and %d more patterns\n", nrow(data) - n_show))
  }
} 
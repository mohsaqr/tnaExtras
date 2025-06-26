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
#' @return Character vector of n-grams
extract_ngrams <- function(sequence, n) {
  if (nchar(sequence) == 0) return(character(0))
  
  seq_parts <- unlist(strsplit(sequence, "-"))
  seq_parts <- seq_parts[nchar(seq_parts) > 1]  # Filter single characters
  
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
  
  if (length(unique(groups)) < 2) {
    stop("At least two groups must be present in the data")
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

#' Compute Support-based Differences for Multiple Groups
#'
#' Support measures how frequently a pattern appears in each group.
#' Support is defined as the proportion of sequences containing the pattern.
#'
#' @param patterns Character vector of patterns to analyze
#' @param group_sequences List of character vectors, one for each group
#' @param group_names Character vector of group names
#' @return Data frame with support measures
compute_support_measures_multi <- function(patterns, group_sequences, group_names) {
  # Input validation
  if (!is.character(patterns) || length(patterns) == 0) {
    stop("patterns must be a non-empty character vector")
  }
  if (!is.list(group_sequences) || length(group_sequences) == 0) {
    stop("group_sequences must be a non-empty list")
  }
  if (length(group_names) != length(group_sequences)) {
    stop("group_names must have the same length as group_sequences")
  }
  
  n_groups <- length(group_sequences)
  n_patterns <- length(patterns)
  
  # Initialize results dataframe
  results <- data.frame(
    pattern = patterns,
    stringsAsFactors = FALSE
  )
  
  # Add support columns for each group
  for (i in 1:n_groups) {
    col_name <- paste0("support_", group_names[i])
    results[[col_name]] <- numeric(n_patterns)
  }
  
  # Add overall statistics columns
  results$max_support <- numeric(n_patterns)
  results$min_support <- numeric(n_patterns)
  results$support_range <- numeric(n_patterns)
  results$support_variance <- numeric(n_patterns)
  results$dominant_group <- character(n_patterns)
  
  # Get group sizes
  group_sizes <- sapply(group_sequences, length)
  
  # Compute support for each pattern
  for (i in seq_along(patterns)) {
    pattern <- patterns[i]
    supports <- numeric(n_groups)
    
    # Calculate support for each group
    for (j in 1:n_groups) {
      if (group_sizes[j] > 0) {
        count <- sum(grepl(pattern, group_sequences[[j]], fixed = TRUE))
        support <- count / group_sizes[j]
        supports[j] <- support
        col_name <- paste0("support_", group_names[j])
        results[[col_name]][i] <- support
      }
    }
    
    # Calculate summary statistics
    results$max_support[i] <- max(supports)
    results$min_support[i] <- min(supports)
    results$support_range[i] <- max(supports) - min(supports)
    results$support_variance[i] <- var(supports)
    results$dominant_group[i] <- group_names[which.max(supports)]
  }
  
  # Sort by support range (difference between max and min)
  results <- results[order(results$support_range, decreasing = TRUE), ]
  rownames(results) <- NULL
  
  return(results)
}

# Keep the original function for backward compatibility
compute_support_measures <- function(patterns, group_A_seqs, group_B_seqs, group_names = c("A", "B")) {
  group_sequences <- list(group_A_seqs, group_B_seqs)
  
  multi_result <- compute_support_measures_multi(patterns, group_sequences, group_names)
  
  # Create dynamic column names using actual group names
  support_col_A <- paste0("support_", group_names[1])
  support_col_B <- paste0("support_", group_names[2])
  
  # Convert to original format with actual group names
  results <- data.frame(
    pattern = multi_result$pattern,
    stringsAsFactors = FALSE
  )
  
  # Add support columns with actual group names
  results[[support_col_A]] <- multi_result[[support_col_A]]
  results[[support_col_B]] <- multi_result[[support_col_B]]
  
  # Calculate differences and ratios
  results$support_diff <- multi_result[[support_col_A]] - multi_result[[support_col_B]]
  results$support_ratio <- (multi_result[[support_col_A]] + 0.001) / (multi_result[[support_col_B]] + 0.001)
  results$relative_support <- ifelse((multi_result[[support_col_A]] + multi_result[[support_col_B]]) == 0, 0, 
                           (multi_result[[support_col_A]] - multi_result[[support_col_B]]) / 
                           (multi_result[[support_col_A]] + multi_result[[support_col_B]]))
  
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
#' @param group_names Character vector of actual group names (default: c("A", "B"))
#' @return Data frame with lift measures
compute_lift_measures <- function(patterns, group_A_seqs, group_B_seqs, group_names = c("A", "B")) {
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
  
  # Create dynamic column names using actual group names
  observed_col_A <- paste0("observed_", group_names[1])
  observed_col_B <- paste0("observed_", group_names[2])
  expected_col_A <- paste0("expected_", group_names[1])
  expected_col_B <- paste0("expected_", group_names[2])
  lift_col_A <- paste0("lift_", group_names[1])
  lift_col_B <- paste0("lift_", group_names[2])
  
  results <- data.frame(
    pattern = patterns,
    stringsAsFactors = FALSE
  )
  
  # Add columns with actual group names
  results[[observed_col_A]] <- numeric(length(patterns))
  results[[observed_col_B]] <- numeric(length(patterns))
  results[[expected_col_A]] <- numeric(length(patterns))
  results[[expected_col_B]] <- numeric(length(patterns))
  results[[lift_col_A]] <- numeric(length(patterns))
  results[[lift_col_B]] <- numeric(length(patterns))
  results$lift_ratio <- numeric(length(patterns))
  results$lift_difference <- numeric(length(patterns))
  
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
      next
    }
    
    # Expected counts under independence (null hypothesis)
    expected_A <- count_total * (n_A / n_total)
    expected_B <- count_total * (n_B / n_total)
    
    # Observed proportions
    prop_A <- count_A / n_A
    prop_B <- count_B / n_B
    
    # Expected proportions under independence
    expected_prop_A <- count_total / n_total * (n_A / n_total)
    expected_prop_B <- count_total / n_total * (n_B / n_total)
    
    # Lift = observed / expected (with small constant to avoid division by zero)
    lift_A <- prop_A / (expected_prop_A + 1e-10)
    lift_B <- prop_B / (expected_prop_B + 1e-10)
    
    results[[observed_col_A]][i] <- count_A
    results[[observed_col_B]][i] <- count_B
    results[[expected_col_A]][i] <- expected_A
    results[[expected_col_B]][i] <- expected_B
    results[[lift_col_A]][i] <- lift_A
    results[[lift_col_B]][i] <- lift_B
    
    # Lift ratio (capped to avoid extreme values)
    lift_ratio <- pmin(pmax(lift_A / (lift_B + 1e-10), 0.01), 100)
    results$lift_ratio[i] <- lift_ratio
    
    # Lift difference (more intuitive)
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
#' @param group_names Character vector of actual group names (default: c("A", "B"))
#' @return Data frame with confidence measures
compute_confidence_measures <- function(patterns, group_A_seqs, group_B_seqs, group_names = c("A", "B")) {
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
  
  # Create dynamic column names using actual group names
  confidence_col_A <- paste0("confidence_", group_names[1])
  confidence_col_B <- paste0("confidence_", group_names[2])
  
  results <- data.frame(
    pattern = patterns,
    stringsAsFactors = FALSE
  )
  
  # Add columns with actual group names
  results[[confidence_col_A]] <- numeric(length(patterns))
  results[[confidence_col_B]] <- numeric(length(patterns))
  results$confidence_diff <- numeric(length(patterns))
  results$confidence_ratio <- numeric(length(patterns))
  results$balanced_confidence <- numeric(length(patterns))
  
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
      next
    }
    
    # Confidence = P(pattern | group)
    conf_A <- count_A / count_total
    conf_B <- count_B / count_total
    
    results[[confidence_col_A]][i] <- conf_A
    results[[confidence_col_B]][i] <- conf_B
    results$confidence_diff[i] <- conf_A - conf_B
    
    # Confidence ratio with small constant
    results$confidence_ratio[i] <- (conf_A + 0.001) / (conf_B + 0.001)
    
    # Balanced confidence (accounts for group sizes)
    expected_A <- n_A / (n_A + n_B)
    expected_B <- n_B / (n_A + n_B)
    
    results$balanced_confidence[i] <- (conf_A - expected_A) - (conf_B - expected_B)
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
#' @param group_names Character vector of actual group names (default: c("A", "B"))
#' @return Data frame with effect size measures
compute_effect_sizes <- function(patterns, group_A_seqs, group_B_seqs, group_names = c("A", "B")) {
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
    
    # Cohen's h (effect size for difference in proportions)
    if (p_A == 0 && p_B == 0) {
      cohens_h <- 0
    } else {
      # Use arcsine transformation
      h_A <- 2 * asin(sqrt(pmax(0, pmin(1, p_A))))
      h_B <- 2 * asin(sqrt(pmax(0, pmin(1, p_B))))
      cohens_h <- h_A - h_B
    }
    
    # Cohen's d (standardized mean difference)
    # For binary data, use pooled standard deviation
    pooled_p <- (count_A + count_B) / (n_A + n_B)
    pooled_var <- pooled_p * (1 - pooled_p)
    if (pooled_var > 0) {
      cohens_d <- (p_A - p_B) / sqrt(pooled_var)
    } else {
      cohens_d <- 0
    }
    
    # Create contingency table for chi-square based measures
    # Present/Absent in each group
    cont_table <- matrix(c(
      count_A, n_A - count_A,      # Group A: present, absent
      count_B, n_B - count_B       # Group B: present, absent
    ), nrow = 2, byrow = TRUE)
    
    # Cramer's V
    total_n <- n_A + n_B
    if (total_n > 0) {
      chi_sq <- suppressWarnings(chisq.test(cont_table, correct = FALSE))
      if (!is.na(chi_sq$statistic)) {
        cramers_v <- sqrt(chi_sq$statistic / (total_n * (min(nrow(cont_table), ncol(cont_table)) - 1)))
      } else {
        cramers_v <- 0
      }
    } else {
      cramers_v <- 0
    }
    
    # Phi coefficient (for 2x2 tables)
    if (total_n > 0) {
      a <- cont_table[1,1]  # Group A, present
      b <- cont_table[1,2]  # Group A, absent
      c <- cont_table[2,1]  # Group B, present
      d <- cont_table[2,2]  # Group B, absent
      
      denom <- sqrt((a + b) * (c + d) * (a + c) * (b + d))
      if (denom > 0) {
        phi_coefficient <- (a * d - b * c) / denom
      } else {
        phi_coefficient <- 0
      }
    } else {
      phi_coefficient <- 0
    }
    
    # Standardized difference
    pooled_se <- sqrt(pooled_p * (1 - pooled_p) * (1/n_A + 1/n_B))
    
    results$cohens_h[i] <- cohens_h
    results$cohens_d[i] <- cohens_d
    results$cramers_v[i] <- cramers_v
    results$phi_coefficient[i] <- phi_coefficient
    results$standardized_diff[i] <- ifelse(pooled_se == 0, 0, (p_A - p_B) / pooled_se)
  }
  
  # Sort by Cohen's h
  results <- results[order(abs(results$cohens_h), decreasing = TRUE), ]
  rownames(results) <- NULL
  
  return(results)
}

# ==============================================================================
# MULTI-GROUP ANALYSIS FUNCTION
# ==============================================================================

#' Comprehensive Pattern Analysis for Multiple Groups
#'
#' Compute pattern analysis across multiple groups in sequential data.
#' Supports both data.frame input and group_tna objects from the tna package.
#'
#' @param data Data frame with sequences in wide format OR a group_tna object
#' @param group_col Column name or index containing group information (default: "Group").
#'   Ignored if data is a group_tna object.
#' @param min_length Minimum pattern length to analyze (default: 2)
#' @param max_length Maximum pattern length to analyze (default: 5)
#' @param min_frequency Minimum frequency required to include a pattern (default: 2)
#' @param measures Character vector of measures to compute (default: "support")
#' @return List containing analysis results
analyze_patterns_multi <- function(data, group_col = "Group", min_length = 2, max_length = 5,
                                   min_frequency = 2, measures = c("support")) {
  
  # =====================================================================
  # INPUT VALIDATION AND GROUP_TNA SUPPORT
  # =====================================================================
  
  # Check if input is a group_tna object
  if (is_group_tna(data)) {
    cat("Detected group_tna object, converting to tnaExtras format...\n")
    converted <- convert_group_tna(data)
    data <- converted$data
    group_col <- converted$group_col
    group_info <- converted$group_info
    
    cat("Successfully converted group_tna object:\n")
    cat("  Label:", group_info$label %||% "Unknown", "\n")
    cat("  Groups:", paste(group_info$levels, collapse = ", "), "\n")
    cat("  Total sequences:", nrow(data), "\n")
  } else {
    group_info <- NULL
  }
  
  # Input validation
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame or group_tna object")
  }
  if (nrow(data) == 0) {
    stop("'data' cannot be empty")
  }
  
  # Prepare data
  cat("Preparing sequence data...\n")
  prepared_data <- prepare_sequence_data(data, group_col, min_length)
  
  all_sequences <- prepared_data$sequences
  groups <- prepared_data$groups
  group_names <- prepared_data$group_names
  
  # Split sequences by group
  group_sequences <- list()
  for (g in group_names) {
    group_sequences[[g]] <- all_sequences[groups == g]
    cat(sprintf("Group %s: %d sequences\n", g, length(group_sequences[[g]])))
  }
  
  # Extract all unique patterns
  cat("Extracting patterns...\n")
  all_patterns <- character(0)
  for (length in min_length:max_length) {
    for (seq in all_sequences) {
      if (nchar(seq) > 0) {
        patterns <- extract_ngrams(seq, length)
        all_patterns <- c(all_patterns, patterns)
      }
    }
  }
  
  # Filter by frequency
  pattern_counts <- table(all_patterns)
  frequent_patterns <- names(pattern_counts[pattern_counts >= min_frequency])
  
  if (length(frequent_patterns) == 0) {
    max_freq <- if(length(pattern_counts) > 0) max(pattern_counts) else 0
    stop("No patterns meet the minimum frequency threshold of ", min_frequency, 
         ". Maximum pattern frequency found: ", max_freq, 
         ". Try reducing min_frequency to ", max(1, max_freq), " or lower.")
  }
  
  cat(sprintf("Found %d patterns meeting frequency threshold\n", length(frequent_patterns)))
  
  # Compute measures
  results <- list()
  
  valid_measures <- c("support")  # Can extend with other multi-group measures
  if (!all(measures %in% valid_measures)) {
    stop("measures must be one or more of: ", paste(valid_measures, collapse = ", "))
  }
  
  if ("support" %in% measures) {
    cat("Computing multi-group support measures...\n")
    results$support <- compute_support_measures_multi(frequent_patterns, group_sequences, group_names)
  }
  
  # Add metadata
  group_sizes <- sapply(group_sequences, length)
  names(group_sizes) <- group_names
  
  results$metadata <- list(
    group_names = group_names,
    n_groups = length(group_names),
    group_sizes = group_sizes,
    n_patterns = length(frequent_patterns),
    min_length = min_length,
    max_length = max_length,
    min_frequency = min_frequency,
    group_tna_info = group_info  # Store original group_tna metadata
  )
  
  cat("Multi-group analysis complete!\n")
  
  class(results) <- "pattern_analysis_multi"
  return(results)
}

# ==============================================================================
# MAIN ANALYSIS FUNCTION (keeping for backward compatibility)
# ==============================================================================

#' Comprehensive Pattern Analysis - Two Groups
#'
#' Compute comprehensive pattern analysis between two groups in sequential data.
#' This function provides detailed statistical measures for pattern comparison.
#' Supports both data.frame input and group_tna objects from the tna package.
#'
#' @param data Data frame with sequences in wide format OR a group_tna object
#' @param group_col Column name or index containing group information (default: "Group").
#'   Ignored if data is a group_tna object.
#' @param min_length Minimum pattern length to analyze (default: 2)
#' @param max_length Maximum pattern length to analyze (default: 5)
#' @param min_frequency Minimum frequency required to include a pattern (default: 2). 
#'   If too high, may result in "No patterns meet threshold" error. Try reducing if this occurs.
#' @param measures Character vector of measures to compute (default: all)
#' @param statistical Whether to include statistical testing (default: FALSE)
#' @param correction Multiple comparison correction method (default: "bonferroni")
#' @return List containing analysis results
analyze_patterns <- function(data, group_col = "Group", min_length = 2, max_length = 5,
                            min_frequency = 2, measures = c("support", "lift", "confidence", "effect_size"), statistical = FALSE,
                            correction = "bonferroni") {
  
  # =====================================================================
  # INPUT VALIDATION AND GROUP_TNA SUPPORT
  # =====================================================================
  
  # Check if input is a group_tna object
  if (is_group_tna(data)) {
    cat("Detected group_tna object, converting to tnaExtras format...\n")
    converted <- convert_group_tna(data)
    data <- converted$data
    group_col <- converted$group_col
    group_info <- converted$group_info
    
    cat("Successfully converted group_tna object:\n")
    cat("  Label:", group_info$label %||% "Unknown", "\n")
    cat("  Groups:", paste(group_info$levels, collapse = ", "), "\n")
    cat("  Total sequences:", nrow(data), "\n")
  } else {
    group_info <- NULL
  }
  
  # Input validation
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame or group_tna object")
  }
  if (nrow(data) == 0) {
    stop("'data' cannot be empty")
  }
  
  # Prepare data
  cat("Preparing sequence data...\n")
  prepared_data <- prepare_sequence_data(data, group_col, min_length)
  
  all_sequences <- prepared_data$sequences
  groups <- prepared_data$groups
  group_names <- prepared_data$group_names
  
  # Check for multiple groups and suggest alternative
  if (length(group_names) > 2) {
    warning("Detected ", length(group_names), " groups. Consider using analyze_patterns_multi() for better multi-group analysis.")
    cat("Note: Using first 2 groups for traditional analysis: ", group_names[1], " and ", group_names[2], "\n")
    # Use only first 2 groups
    keep_groups <- groups %in% group_names[1:2]
    all_sequences <- all_sequences[keep_groups]
    groups <- groups[keep_groups]
    group_names <- group_names[1:2]
  }
  
  # Split sequences by group
  group_A_sequences <- all_sequences[groups == group_names[1]]
  group_B_sequences <- all_sequences[groups == group_names[2]]
  
  cat(sprintf("Group %s: %d sequences\n", group_names[1], length(group_A_sequences)))
  cat(sprintf("Group %s: %d sequences\n", group_names[2], length(group_B_sequences)))
  
  # Extract all unique patterns
  cat("Extracting patterns...\n")
  all_patterns <- character(0)
  for (length in min_length:max_length) {
    for (seq in all_sequences) {
      if (nchar(seq) > 0) {
        patterns <- extract_ngrams(seq, length)
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
  
  # Handle "all" parameter
  if (length(measures) == 1 && measures[1] == "all") {
    measures <- valid_measures
  }
  
  if (!all(measures %in% valid_measures)) {
    stop("measures must be one or more of: ", paste(valid_measures, collapse = ", "))
  }
  
  if ("support" %in% measures) {
    cat("Computing support measures...\n")
    results$support <- compute_support_measures(frequent_patterns, group_A_sequences, group_B_sequences, group_names)
  }
  
  if ("lift" %in% measures) {
    cat("Computing lift measures...\n")
    results$lift <- compute_lift_measures(frequent_patterns, group_A_sequences, group_B_sequences, group_names)
  }
  
  if ("confidence" %in% measures) {
    cat("Computing confidence measures...\n")
    results$confidence <- compute_confidence_measures(frequent_patterns, group_A_sequences, group_B_sequences, group_names)
  }
  
  if ("effect_size" %in% measures) {
    cat("Computing effect size measures...\n")
    results$effect_size <- compute_effect_sizes(frequent_patterns, group_A_sequences, group_B_sequences, group_names)
  }
  
  # Add metadata
  results$metadata <- list(
    group_names = group_names,
    n_sequences = c(length(group_A_sequences), length(group_B_sequences)),
    n_patterns = length(frequent_patterns),
    min_length = min_length,
    max_length = max_length,
    min_frequency = min_frequency,
    group_tna_info = group_info  # Store original group_tna metadata
  )
  
  cat("Analysis complete!\n")
  
  class(results) <- "pattern_analysis"
  return(results)
}

# ==============================================================================
# PRINT AND SUMMARY METHODS
# ==============================================================================

#' Print method for pattern_analysis_multi objects
#' @param x pattern_analysis_multi object
#' @param ... additional arguments
print.pattern_analysis_multi <- function(x, ...) {
  cat("Multi-Group Pattern Analysis Results\n")
  cat("====================================\n\n")
  
  if (!is.null(x$metadata)) {
    cat("Groups:", paste(x$metadata$group_names, collapse = ", "), "\n")
    cat("Group sizes:\n")
    for (i in seq_along(x$metadata$group_names)) {
      cat(sprintf("  %s: %d sequences\n", x$metadata$group_names[i], x$metadata$group_sizes[i]))
    }
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

#' Summary method for pattern_analysis_multi objects
#' @param object pattern_analysis_multi object
#' @param measure which measure to summarize
#' @param top_n number of top patterns to show
#' @param ... additional arguments
summary.pattern_analysis_multi <- function(object, measure = NULL, top_n = 10, ...) {
  available_measures <- names(object)[names(object) != "metadata"]
  
  if (is.null(measure)) {
    measure <- available_measures[1]
  }
  
  if (!measure %in% available_measures) {
    stop("Measure '", measure, "' not found. Available: ", 
         paste(available_measures, collapse = ", "))
  }
  
  cat("Multi-Group Pattern Analysis Summary -", toupper(measure), "Measures\n")
  cat(paste(rep("=", 60), collapse = ""), "\n\n")
  
  data <- object[[measure]]
  n_show <- min(top_n, nrow(data))
  
  cat("Top", n_show, "patterns by discrimination:\n")
  print(head(data, n_show))
  
  if (nrow(data) > n_show) {
    cat(sprintf("\n... and %d more patterns\n", nrow(data) - n_show))
  }
}

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
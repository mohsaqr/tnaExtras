# ==============================================================================
# SEQUENCE PATTERN EXPLORATION TOOLKIT
# ==============================================================================
#
# Unified function for exploring sequence patterns in wide-format data.
# Supports n-gram extraction, abstract patterns, and full sequence analysis.
#
# Main function: explore_patterns()
# Helper functions: significant_patterns(), filter_patterns(), top_sequences()
#
# ==============================================================================

#' Explore Sequence Patterns
#'
#' Unified function for discovering and analyzing patterns in sequential data.
#' Supports n-gram extraction, abstract structural patterns (returns, repetitions,
#' oscillations, progressions), and full sequence frequency analysis.
#'
#' @param data Data frame with sequences in wide format (rows are sequences,
#'   columns are time points). Can also be a group_tna object.
#' @param type Type of pattern analysis to perform:
#'   \itemize{
#'     \item \code{"ngrams"} - Extract n-grams of specified lengths (default)
#'     \item \code{"abstract"} - Detect abstract patterns (returns, repetitions,
#'       oscillations, progressions)
#'     \item \code{"full"} - Analyze complete sequence frequencies
#'   }
#' @param min_length Minimum pattern length / n-gram size (default: 1)
#' @param max_length Maximum pattern length / n-gram size (default: 5)
#' @param min_gap Minimum gap size for abstract return patterns (default: 1)
#' @param max_gap Maximum gap size for abstract return patterns (default: 5)
#' @param min_support Minimum support threshold (default: 0.01)
#' @param min_count Minimum count threshold (default: 2)
#' @param test_significance Whether to perform significance tests (default: TRUE)
#' @param correction Multiple testing correction method (default: "fdr").
#'   Options: "bonferroni", "holm", "hochberg", "hommel", "BH", "BY", "fdr", "none"
#' @param alpha Significance level (default: 0.05)
#' @param start_state Filter patterns starting with this state (default: NULL)
#' @param end_state Filter patterns ending with this state (default: NULL)
#' @param contains_state Filter patterns containing this state (default: NULL)
#' @param verbose Whether to print progress messages (default: TRUE)
#'
#' @return An object of class "sequence_patterns" containing:
#'   \item{patterns}{Data frame with pattern statistics}
#'   \item{summary}{Summary statistics}
#'   \item{parameters}{Analysis parameters}
#'
#' @examples
#' \dontrun{
#' # Load regulation data
#' data <- tna::group_regulation
#' seq_data <- data[, -1]
#'
#' # Basic n-gram analysis
#' results <- explore_patterns(seq_data)
#' print(results)
#'
#' # Filter patterns starting with "plan"
#' results <- explore_patterns(seq_data, start_state = "plan")
#'
#' # Abstract pattern detection
#' abstract <- explore_patterns(seq_data, type = "abstract")
#'
#' # Full sequence analysis
#' full <- explore_patterns(seq_data, type = "full")
#' }
#'
#' @export
explore_patterns <- function(data,
                             type = "ngrams",
                             min_length = 1,
                             max_length = 5,
                             min_gap = 1,
                             max_gap = 5,
                             min_support = 0.01,
                             min_count = 2,
                             test_significance = TRUE,
                             correction = "fdr",
                             alpha = 0.05,
                             start_state = NULL,
                             end_state = NULL,
                             contains_state = NULL,
                             verbose = TRUE) {
  
  # Validate type
  valid_types <- c("ngrams", "abstract", "full")
  if (!type %in% valid_types) {
    stop("type must be one of: ", paste(valid_types, collapse = ", "))
  }
  
  # Validate common parameters
  validate_sequence_params(list(
    min_length = min_length, max_length = max_length,
    min_support = min_support, min_count = min_count,
    correction = correction
  ))
  
  # Prepare data using shared utility
  seq_data <- prepare_sequence_input(data, verbose)
  
  # Dispatch to appropriate analysis function
  patterns_df <- switch(type,
    "ngrams" = extract_ngram_patterns(seq_data, min_length, max_length, 
                                       min_support, min_count, verbose),
    "abstract" = extract_abstract_patterns(seq_data, min_gap, max_gap,
                                           min_support, min_count, verbose),
    "full" = extract_full_sequences(seq_data, min_support, min_count, verbose)
  )
  
  # Apply state filters using unified utility
  if (!is.null(start_state) || !is.null(end_state) || !is.null(contains_state)) {
    if (verbose) cat("Applying state filters...\n")
    patterns_df <- filter_patterns_by_text(patterns_df, start_state, end_state, 
                                           contains_state, text_col = "pattern")
    if (verbose) cat("  Patterns after filtering:", nrow(patterns_df), "\n")
  }
  
  # Compute significance using unified utility
  if (test_significance && nrow(patterns_df) > 0) {
    if (verbose) cat("Computing significance statistics...\n")
    patterns_df <- compute_significance_stats(patterns_df, seq_data$n_sequences,
                                              seq_data$n_states, correction, alpha)
  }
  
  # Standardize columns
  patterns_df <- standardize_columns(patterns_df, include_schema = FALSE)
  
  # Build result object
  result <- list(
    patterns = patterns_df,
    summary = list(
      n_sequences = seq_data$n_sequences,
      n_states = seq_data$n_states,
      n_patterns = nrow(patterns_df),
      n_significant = if (test_significance && "significant" %in% names(patterns_df))
        sum(patterns_df$significant, na.rm = TRUE) else NA,
      type = type,
      all_states = seq_data$all_states
    ),
    parameters = list(
      type = type,
      min_length = min_length, max_length = max_length,
      min_gap = min_gap, max_gap = max_gap,
      min_support = min_support, min_count = min_count,
      test_significance = test_significance,
      correction = correction, alpha = alpha,
      start_state = start_state, end_state = end_state, contains_state = contains_state
    )
  )
  
  class(result) <- "sequence_patterns"
  
  if (verbose) {
    cat("\n=== Analysis Complete ===\n")
    cat("Type:", type, "\n")
    cat("Patterns found:", nrow(patterns_df), "\n")
    if (test_significance && "significant" %in% names(patterns_df)) {
      cat("Significant patterns:", sum(patterns_df$significant, na.rm = TRUE), "\n")
    }
  }
  
  return(result)
}

# ==============================================================================
# PATTERN EXTRACTION FUNCTIONS
# ==============================================================================

#' Extract n-gram patterns
#' @keywords internal
extract_ngram_patterns <- function(seq_data, min_length, max_length,
                                   min_support, min_count, verbose) {
  
  if (verbose) cat("Extracting n-gram patterns...\n")
  
  sequences <- seq_data$sequences
  n_sequences <- seq_data$n_sequences
  
  all_ngrams <- list()
  ngram_seqs <- list()
  
  for (n in min_length:max_length) {
    if (verbose) cat("  Extracting", n, "-grams...\n")
    
    for (seq_idx in seq_along(sequences)) {
      seq <- sequences[[seq_idx]]
      if (length(seq) < n) next
      
      for (start in 1:(length(seq) - n + 1)) {
        ngram <- paste(seq[start:(start + n - 1)], collapse = "->")
        
        if (is.null(all_ngrams[[ngram]])) {
          all_ngrams[[ngram]] <- 0
          ngram_seqs[[ngram]] <- c()
        }
        all_ngrams[[ngram]] <- all_ngrams[[ngram]] + 1
        ngram_seqs[[ngram]] <- unique(c(ngram_seqs[[ngram]], seq_idx))
      }
    }
  }
  
  if (length(all_ngrams) == 0) {
    return(create_empty_patterns_df())
  }
  
  # Build data frame using aggregation utility
  patterns_df <- aggregate_instances(
    rep(names(all_ngrams), unlist(all_ngrams)), # Reconstruct raw list
    ngram_seqs, 
    n_sequences, 
    seq_data$n_states,
    test_significance = FALSE # Stats added later in main function
  )
  
  filter_by_thresholds(patterns_df, min_support, min_count)
}

#' Extract abstract patterns (returns, repetitions, oscillations, progressions)
#' @keywords internal
extract_abstract_patterns <- function(seq_data, min_gap, max_gap,
                                      min_support, min_count, verbose) {
  
  if (verbose) cat("Extracting abstract patterns...\n")
  
  sequences <- seq_data$sequences
  n_sequences <- seq_data$n_sequences
  all_states <- seq_data$all_states
  
  patterns_list <- list()
  
  # Helper to check threshold and add pattern
  add_if_valid <- function(name, count, seqs_n, type, len) {
    if (count >= min_count && seqs_n / n_sequences >= min_support) {
      patterns_list[[length(patterns_list) + 1]] <<- list(
        pattern = name, type = type, length = len,
        count = count, sequences_containing = seqs_n
      )
    }
  }
  
  # 1. Returns (A->*->A)
  if (verbose) cat("  Detecting returns...\n")
  for (state in all_states) {
    for (gap in min_gap:max_gap) {
      count <- 0; seqs_containing <- 0
      for (seq in sequences) {
        positions <- which(seq == state)
        if (length(positions) < 2) next
        seq_matches <- 0
        for (i in 1:(length(positions) - 1)) {
          if (any(positions[(i+1):length(positions)] - positions[i] - 1 == gap)) seq_matches <- seq_matches + 1
        }
        if (seq_matches > 0) { count <- count + seq_matches; seqs_containing <- seqs_containing + 1 }
      }
      add_if_valid(paste0(state, "->(*", gap, ")->", state), count, seqs_containing, "return", gap + 2)
    }
  }
  
  # 2. Repetitions (A->A, A->A->A)
  if (verbose) cat("  Detecting repetitions...\n")
  for (state in all_states) {
    for (rep_len in 2:min(5, max_gap + 2)) {
      count <- 0; seqs_containing <- 0
      for (seq in sequences) {
        if (length(seq) < rep_len) next
        seq_matches <- 0
        for (i in 1:(length(seq) - rep_len + 1)) {
          if (all(seq[i:(i + rep_len - 1)] == state)) seq_matches <- seq_matches + 1
        }
        if (seq_matches > 0) { count <- count + seq_matches; seqs_containing <- seqs_containing + 1 }
      }
      add_if_valid(paste(rep(state, rep_len), collapse = "->"), count, seqs_containing, "repetition", rep_len)
    }
  }
  
  # 3. Oscillations (A->B->A->B)
  if (verbose) cat("  Detecting oscillations...\n")
  for (s1 in all_states) {
    for (s2 in all_states) {
      if (s1 >= s2) next
      count <- 0; seqs_containing <- 0
      for (seq in sequences) {
        if (length(seq) < 4) next
        seq_matches <- 0
        for (i in 1:(length(seq) - 3)) {
          if (seq[i] == s1 && seq[i+1] == s2 && seq[i+2] == s1 && seq[i+3] == s2) seq_matches <- seq_matches + 1
        }
        if (seq_matches > 0) { count <- count + seq_matches; seqs_containing <- seqs_containing + 1 }
      }
      add_if_valid(paste(c(s1, s2, s1, s2), collapse = "->"), count, seqs_containing, "oscillation", 4)
    }
  }
  
  # 4. Progressions (unique state chains)
  if (verbose) cat("  Detecting progressions...\n")
  prog_counts <- list(); prog_seqs <- list()
  
  for (seq_idx in seq_along(sequences)) {
    seq <- sequences[[seq_idx]]
    i <- 1
    while (i <= length(seq)) {
      j <- i
      while (j < length(seq) && seq[j + 1] != seq[j] && !seq[j + 1] %in% seq[i:j]) j <- j + 1
      if (j - i + 1 >= 3) {
        prog <- paste(seq[i:j], collapse = "->")
        if (is.null(prog_counts[[prog]])) { prog_counts[[prog]] <- 0; prog_seqs[[prog]] <- c() }
        prog_counts[[prog]] <- prog_counts[[prog]] + 1
        prog_seqs[[prog]] <- unique(c(prog_seqs[[prog]], seq_idx))
      }
      i <- j + 1
    }
  }
  
  for (prog in names(prog_counts)) {
    add_if_valid(prog, prog_counts[[prog]], length(prog_seqs[[prog]]), "progression", length(strsplit(prog, "->")[[1]]))
  }
  
  if (length(patterns_list) == 0) return(create_empty_patterns_df())
  
  # Build data frame
  patterns_df <- data.frame(
    pattern = sapply(patterns_list, `[[`, "pattern"),
    length = sapply(patterns_list, `[[`, "length"),
    count = sapply(patterns_list, `[[`, "count"),
    sequences_containing = sapply(patterns_list, `[[`, "sequences_containing"),
    stringsAsFactors = FALSE
  )
  
  patterns_df$support <- round(patterns_df$sequences_containing / n_sequences, 4)
  patterns_df$proportion <- round(patterns_df$count / sum(patterns_df$count), 4)
  
  # Sort and return
  patterns_df <- patterns_df[order(patterns_df$count, decreasing = TRUE), ]
  rownames(patterns_df) <- NULL
  
  return(patterns_df)
}

#' Extract full sequence frequencies
#' @keywords internal
extract_full_sequences <- function(seq_data, min_support, min_count, verbose) {
  
  if (verbose) cat("Extracting full sequence frequencies...\n")
  
  sequences <- seq_data$sequences
  n_sequences <- seq_data$n_sequences
  
  # Convert to strings
  seq_strings <- sapply(sequences, function(s) paste(s, collapse = "->"))
  
  # Count frequencies
  seq_table <- table(seq_strings)
  
  patterns_df <- data.frame(
    pattern = names(seq_table),
    length = sapply(strsplit(names(seq_table), "->"), length),
    count = as.integer(seq_table),
    sequences_containing = as.integer(seq_table),
    stringsAsFactors = FALSE
  )
  
  patterns_df$support <- round(patterns_df$count / n_sequences, 4)
  patterns_df$proportion <- round(patterns_df$count / sum(patterns_df$count), 4)
  
  filter_by_thresholds(patterns_df, min_support, min_count)
}

# ==============================================================================
# S3 METHODS
# ==============================================================================

#' Print method for sequence_patterns objects
#' @export
print.sequence_patterns <- function(x, ...) {
  cat("Sequence Pattern Analysis Results\n")
  cat("=================================\n\n")
  
  cat("Analysis Type:", x$parameters$type, "\n")
  cat("Sequences:", x$summary$n_sequences, "\n")
  cat("Unique States:", x$summary$n_states, "\n")
  cat("Patterns Found:", x$summary$n_patterns, "\n")
  if (!is.na(x$summary$n_significant)) {
    cat("Significant Patterns:", x$summary$n_significant, "\n")
  }
  cat("\n")
  
  print_patterns_summary(x$patterns, "Top Patterns", n_show = 15)
  
  invisible(x)
}

#' Summary method for sequence_patterns objects
#' @export
summary.sequence_patterns <- function(object, ...) {
  cat("Sequence Pattern Analysis - Detailed Summary\n")
  cat("=============================================\n\n")
  
  cat("ANALYSIS CONFIGURATION\n")
  cat(paste(rep("-", 50), collapse = ""), "\n")
  cat("Type:", object$parameters$type, "\n")
  cat("Correction:", object$parameters$correction, "\n\n")
  
  cat("RESULTS\n")
  cat(paste(rep("-", 50), collapse = ""), "\n")
  cat("Patterns found:", object$summary$n_patterns, "\n")
  if (!is.na(object$summary$n_significant)) {
    cat("Significant:", object$summary$n_significant, "\n")
  }
  cat("\n")
  
  if (nrow(object$patterns) > 0) {
    cat("FULL PATTERN STATISTICS\n")
    cat(paste(rep("-", 50), collapse = ""), "\n")
    display_cols <- c("pattern", "length", "count", "sequences_containing", 
                      "support", "lift", "p_adjusted", "significant")
    print(head(object$patterns[, display_cols[display_cols %in% names(object$patterns)], drop = FALSE], 30), row.names = FALSE)
  }
  
  invisible(object)
}

#' Plot method for sequence_patterns objects
#' @export
plot.sequence_patterns <- function(x, type = "patterns", top_n = 20,
                                   show_support = TRUE, col = NULL, ...) {
  
  if (nrow(x$patterns) == 0) {
    message("No patterns to plot")
    return(invisible(NULL))
  }
  
  data <- head(x$patterns[order(x$patterns$count, decreasing = TRUE), ], top_n)
  
  # Truncate long labels
  labels <- sapply(data$pattern, function(s) {
    if (nchar(s) > 40) paste0(substr(s, 1, 37), "...") else s
  })
  
  # Set colors
  if (is.null(col)) {
    col <- if ("significant" %in% names(data)) ifelse(data$significant, "steelblue", "gray70") else "steelblue"
  }
  
  graphics::par(mar = c(12, 4, 4, 2))
  bp <- graphics::barplot(
    data$count, names.arg = labels, las = 2,
    main = paste("Top", min(top_n, nrow(data)), "Patterns"),
    ylab = "Count", col = col, ...
  )
  
  if (show_support && "support" %in% names(data)) {
    graphics::text(x = bp, y = data$count, labels = sprintf("%.1f%%", data$support * 100), pos = 3, cex = 0.7)
  }
  
  invisible(NULL)
}

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

#' Extract significant patterns
#' @export
significant_patterns <- function(x, sort_by = "p_adjusted") {
  if (!inherits(x, "sequence_patterns")) stop("x must be a sequence_patterns object")
  if (!"significant" %in% names(x$patterns)) stop("No significance testing was performed")
  
  sig <- x$patterns[x$patterns$significant, ]
  if (nrow(sig) == 0) {
    message("No significant patterns found")
    return(sig)
  }
  if (sort_by %in% names(sig)) sig <- sig[order(sig[[sort_by]]), ]
  return(sig)
}

#' Filter patterns by criteria
#' @export
filter_patterns <- function(x, min_support = NULL, min_count = NULL,
                            pattern_length = NULL, significant_only = FALSE,
                            start_state = NULL, end_state = NULL,
                            contains_state = NULL) {
  if (!inherits(x, "sequence_patterns")) stop("x must be a sequence_patterns object")
  
  patterns <- x$patterns
  patterns <- filter_by_thresholds(patterns, min_support, min_count)
  
  if (!is.null(pattern_length)) patterns <- patterns[patterns$length %in% pattern_length, ]
  if (significant_only && "significant" %in% names(patterns)) patterns <- patterns[patterns$significant, ]
  
  filter_patterns_by_text(patterns, start_state, end_state, contains_state)
}

#' Get top sequences/patterns
#' @export
top_sequences <- function(x, top_n = 10, min_support = NULL, significant_only = FALSE) {
  head(filter_patterns(x, min_support, significant_only = significant_only), top_n)
}

# ==============================================================================
# BACKWARD COMPATIBILITY
# ==============================================================================

#' @rdname explore_patterns
#' @export
explore_sequence_patterns <- function(data, ...) {
  .Deprecated("explore_patterns")
  explore_patterns(data, ...)
}

#' @rdname explore_patterns
#' @export
detect_abstract_patterns <- function(data, patterns = "all", ...) {
  .Deprecated("explore_patterns", msg = "detect_abstract_patterns is deprecated. Use explore_patterns(..., type = 'abstract') instead.")
  explore_patterns(data, type = "abstract", ...)
}

cat("Sequence pattern toolkit loaded.\n")

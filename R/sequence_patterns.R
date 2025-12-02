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
  
  # Validate parameters
  if (min_length < 1) stop("min_length must be at least 1")
  if (max_length < min_length) stop("max_length must be >= min_length")
  if (min_support < 0 || min_support > 1) stop("min_support must be between 0 and 1")
  if (min_count < 1) stop("min_count must be at least 1")
  
  valid_corrections <- c("bonferroni", "holm", "hochberg", "hommel", "BH", "BY", "fdr", "none")
  if (!correction %in% valid_corrections) {
    stop("correction must be one of: ", paste(valid_corrections, collapse = ", "))
  }
  
  # Prepare data using shared utility
  seq_data <- prepare_sequence_data(data, verbose)
  
  # Dispatch to appropriate analysis function
  patterns_df <- switch(type,
    "ngrams" = extract_ngram_patterns(seq_data, min_length, max_length, 
                                       min_support, min_count, verbose),
    "abstract" = extract_abstract_patterns(seq_data, min_gap, max_gap,
                                           min_support, min_count, verbose),
    "full" = extract_full_sequences(seq_data, min_support, min_count, verbose)
  )
  
  # Apply state filters
  if (!is.null(start_state) || !is.null(end_state) || !is.null(contains_state)) {
    if (verbose) cat("Applying state filters...\n")
    patterns_df <- filter_by_state(patterns_df, start_state, end_state, 
                                   contains_state, pattern_col = "pattern",
                                   separator = "->")
    if (verbose) cat("  Patterns after filtering:", nrow(patterns_df), "\n")
  }
  
  # Compute significance
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
      min_length = min_length,
      max_length = max_length,
      min_gap = min_gap,
      max_gap = max_gap,
      min_support = min_support,
      min_count = min_count,
      test_significance = test_significance,
      correction = correction,
      alpha = alpha,
      start_state = start_state,
      end_state = end_state,
      contains_state = contains_state
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
  ngram_seqs <- list()  # Track which sequences contain each n-gram
  
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
  
  # Build data frame
  patterns_df <- data.frame(
    pattern = names(all_ngrams),
    length = sapply(strsplit(names(all_ngrams), "->"), length),
    count = unlist(all_ngrams),
    sequences_containing = sapply(names(all_ngrams), function(p) length(ngram_seqs[[p]])),
    stringsAsFactors = FALSE
  )
  
  patterns_df$support <- round(patterns_df$sequences_containing / n_sequences, 4)
  patterns_df$proportion <- round(patterns_df$count / sum(patterns_df$count), 4)
  
  # Filter by thresholds
  patterns_df <- filter_by_thresholds(patterns_df, min_support, min_count)
  
  # Sort by count
  patterns_df <- patterns_df[order(patterns_df$count, decreasing = TRUE), ]
  rownames(patterns_df) <- NULL
  
  if (verbose) cat("  Patterns meeting thresholds:", nrow(patterns_df), "\n")
  
  return(patterns_df)
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
  
  # 1. Returns (A->*->A)
  if (verbose) cat("  Detecting returns...\n")
  for (state in all_states) {
    for (gap in min_gap:max_gap) {
      pattern_name <- paste0(state, "->(*", gap, ")->", state)
      count <- 0
      seqs_containing <- 0
      
      for (seq in sequences) {
        positions <- which(seq == state)
        if (length(positions) < 2) next
        
        seq_count <- 0
        for (i in 1:(length(positions) - 1)) {
          for (j in (i + 1):length(positions)) {
            if (positions[j] - positions[i] - 1 == gap) {
              seq_count <- seq_count + 1
            }
          }
        }
        
        if (seq_count > 0) {
          count <- count + seq_count
          seqs_containing <- seqs_containing + 1
        }
      }
      
      if (count >= min_count && seqs_containing / n_sequences >= min_support) {
        patterns_list[[length(patterns_list) + 1]] <- list(
          pattern = pattern_name,
          type = "return",
          length = gap + 2,
          count = count,
          sequences_containing = seqs_containing
        )
      }
    }
  }
  
  # 2. Repetitions (A->A, A->A->A)
  if (verbose) cat("  Detecting repetitions...\n")
  for (state in all_states) {
    for (rep_len in 2:min(5, max_gap + 2)) {
      pattern_name <- paste(rep(state, rep_len), collapse = "->")
      count <- 0
      seqs_containing <- 0
      
      for (seq in sequences) {
        if (length(seq) < rep_len) next
        
        seq_count <- 0
        for (i in 1:(length(seq) - rep_len + 1)) {
          if (all(seq[i:(i + rep_len - 1)] == state)) {
            seq_count <- seq_count + 1
          }
        }
        
        if (seq_count > 0) {
          count <- count + seq_count
          seqs_containing <- seqs_containing + 1
        }
      }
      
      if (count >= min_count && seqs_containing / n_sequences >= min_support) {
        patterns_list[[length(patterns_list) + 1]] <- list(
          pattern = pattern_name,
          type = "repetition",
          length = rep_len,
          count = count,
          sequences_containing = seqs_containing
        )
      }
    }
  }
  
  # 3. Oscillations (A->B->A->B)
  if (verbose) cat("  Detecting oscillations...\n")
  for (state1 in all_states) {
    for (state2 in all_states) {
      if (state1 >= state2) next  # Avoid duplicates
      
      pattern_name <- paste(c(state1, state2, state1, state2), collapse = "->")
      count <- 0
      seqs_containing <- 0
      
      for (seq in sequences) {
        if (length(seq) < 4) next
        
        seq_count <- 0
        for (i in 1:(length(seq) - 3)) {
          if (seq[i] == state1 && seq[i + 1] == state2 &&
              seq[i + 2] == state1 && seq[i + 3] == state2) {
            seq_count <- seq_count + 1
          }
        }
        
        if (seq_count > 0) {
          count <- count + seq_count
          seqs_containing <- seqs_containing + 1
        }
      }
      
      if (count >= min_count && seqs_containing / n_sequences >= min_support) {
        patterns_list[[length(patterns_list) + 1]] <- list(
          pattern = pattern_name,
          type = "oscillation",
          length = 4,
          count = count,
          sequences_containing = seqs_containing
        )
      }
    }
  }
  
  # 4. Progressions (unique state chains)
  if (verbose) cat("  Detecting progressions...\n")
  prog_counts <- list()
  prog_seqs <- list()
  
  for (seq_idx in seq_along(sequences)) {
    seq <- sequences[[seq_idx]]
    
    # Find runs of unique states
    i <- 1
    while (i <= length(seq)) {
      j <- i
      while (j < length(seq) && seq[j + 1] != seq[j] && !seq[j + 1] %in% seq[i:j]) {
        j <- j + 1
      }
      
      if (j - i + 1 >= 3) {  # At least 3 unique states
        prog <- paste(seq[i:j], collapse = "->")
        if (is.null(prog_counts[[prog]])) {
          prog_counts[[prog]] <- 0
          prog_seqs[[prog]] <- c()
        }
        prog_counts[[prog]] <- prog_counts[[prog]] + 1
        prog_seqs[[prog]] <- unique(c(prog_seqs[[prog]], seq_idx))
      }
      
      i <- j + 1
    }
  }
  
  for (prog in names(prog_counts)) {
    if (prog_counts[[prog]] >= min_count && 
        length(prog_seqs[[prog]]) / n_sequences >= min_support) {
      patterns_list[[length(patterns_list) + 1]] <- list(
        pattern = prog,
        type = "progression",
        length = length(strsplit(prog, "->")[[1]]),
        count = prog_counts[[prog]],
        sequences_containing = length(prog_seqs[[prog]])
      )
    }
  }
  
  if (length(patterns_list) == 0) {
    return(create_empty_patterns_df())
  }
  
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
  
  # Sort by count
  patterns_df <- patterns_df[order(patterns_df$count, decreasing = TRUE), ]
  rownames(patterns_df) <- NULL
  
  if (verbose) cat("  Abstract patterns found:", nrow(patterns_df), "\n")
  
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
    sequences_containing = as.integer(seq_table),  # Each unique sequence counted once
    stringsAsFactors = FALSE
  )
  
  patterns_df$support <- round(patterns_df$count / n_sequences, 4)
  patterns_df$proportion <- round(patterns_df$count / sum(patterns_df$count), 4)
  
  # Filter by thresholds
  patterns_df <- filter_by_thresholds(patterns_df, min_support, min_count)
  
  # Sort by count
  patterns_df <- patterns_df[order(patterns_df$count, decreasing = TRUE), ]
  rownames(patterns_df) <- NULL
  
  if (verbose) cat("  Unique sequences:", nrow(patterns_df), "\n")
  
  return(patterns_df)
}

# ==============================================================================
# S3 METHODS
# ==============================================================================

#' Print method for sequence_patterns objects
#' @param x sequence_patterns object
#' @param ... additional arguments
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
  
  # Show filters if applied
  if (!is.null(x$parameters$start_state) || 
      !is.null(x$parameters$end_state) || 
      !is.null(x$parameters$contains_state)) {
    cat("Filters Applied:\n")
    if (!is.null(x$parameters$start_state)) 
      cat("  start_state:", x$parameters$start_state, "\n")
    if (!is.null(x$parameters$end_state)) 
      cat("  end_state:", x$parameters$end_state, "\n")
    if (!is.null(x$parameters$contains_state)) 
      cat("  contains_state:", x$parameters$contains_state, "\n")
    cat("\n")
  }
  
  # Print top patterns
  print_patterns_summary(x$patterns, "Top Patterns", n_show = 15)
  
  cat("Access $patterns for full data frame\n")
  
  invisible(x)
}

#' Summary method for sequence_patterns objects
#' @param object sequence_patterns object
#' @param ... additional arguments
#' @export
summary.sequence_patterns <- function(object, ...) {
  cat("Sequence Pattern Analysis - Detailed Summary\n")
  cat("=============================================\n\n")
  
  cat("ANALYSIS CONFIGURATION\n")
  cat(paste(rep("-", 50), collapse = ""), "\n")
  cat("Type:", object$parameters$type, "\n")
  cat("Length range:", object$parameters$min_length, "-", object$parameters$max_length, "\n")
  cat("Min support:", object$parameters$min_support, "\n")
  cat("Min count:", object$parameters$min_count, "\n")
  cat("Correction:", object$parameters$correction, "\n")
  cat("Alpha:", object$parameters$alpha, "\n\n")
  
  cat("DATA SUMMARY\n")
  cat(paste(rep("-", 50), collapse = ""), "\n")
  cat("Sequences:", object$summary$n_sequences, "\n")
  cat("Unique states:", object$summary$n_states, "\n")
  cat("States:", paste(head(object$summary$all_states, 10), collapse = ", "))
  if (object$summary$n_states > 10) cat(" ...")
  cat("\n\n")
  
  cat("RESULTS\n")
  cat(paste(rep("-", 50), collapse = ""), "\n")
  cat("Patterns found:", object$summary$n_patterns, "\n")
  if (!is.na(object$summary$n_significant)) {
    cat("Significant:", object$summary$n_significant, 
        sprintf("(%.1f%%)\n", 100 * object$summary$n_significant / max(1, object$summary$n_patterns)))
  }
  cat("\n")
  
  if (nrow(object$patterns) > 0) {
    cat("FULL PATTERN STATISTICS\n")
    cat(paste(rep("-", 50), collapse = ""), "\n")
    
    display_cols <- c("pattern", "length", "count", "sequences_containing", 
                      "support", "lift", "p_adjusted", "significant")
    display_cols <- display_cols[display_cols %in% names(object$patterns)]
    
    print(head(object$patterns[, display_cols, drop = FALSE], 30), row.names = FALSE)
    
    if (nrow(object$patterns) > 30) {
      cat(sprintf("\n... and %d more patterns\n", nrow(object$patterns) - 30))
    }
  }
  
  invisible(object)
}

#' Plot method for sequence_patterns objects
#' @param x sequence_patterns object
#' @param type Type of plot: "patterns" or "bars"
#' @param top_n Number of top items to plot
#' @param show_support Show support percentages on bars
#' @param col Bar color(s)
#' @param ... Additional arguments passed to barplot
#' @export
plot.sequence_patterns <- function(x, type = "patterns", top_n = 20,
                                   show_support = TRUE, col = NULL, ...) {
  
  if (nrow(x$patterns) == 0) {
    message("No patterns to plot")
    return(invisible(NULL))
  }
  
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par))
  
  data <- head(x$patterns[order(x$patterns$count, decreasing = TRUE), ], top_n)
  
  # Truncate long labels
  labels <- sapply(data$pattern, function(s) {
    if (nchar(s) > 40) paste0(substr(s, 1, 37), "...") else s
  })
  
  # Set colors based on significance
  if (is.null(col)) {
    if ("significant" %in% names(data)) {
      col <- ifelse(data$significant, "steelblue", "gray70")
    } else {
      col <- "steelblue"
    }
  }
  
  graphics::par(mar = c(12, 4, 4, 2))
  bp <- graphics::barplot(
    data$count,
    names.arg = labels,
    las = 2,
    main = paste("Top", min(top_n, nrow(data)), "Patterns"),
    ylab = "Count",
    col = col,
    ...
  )
  
  # Add support labels
  if (show_support && "support" %in% names(data)) {
    graphics::text(
      x = bp,
      y = data$count,
      labels = sprintf("%.1f%%", data$support * 100),
      pos = 3,
      cex = 0.7
    )
  }
  
  # Add legend if showing significance
  if ("significant" %in% names(data) && any(data$significant) && is.null(list(...)$col)) {
    graphics::legend("topright",
                     legend = c("Significant", "Not significant"),
                     fill = c("steelblue", "gray70"),
                     bty = "n")
  }
  
  invisible(NULL)
}

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

#' Extract significant patterns
#' @param x sequence_patterns object
#' @param sort_by Column to sort by
#' @return Data frame of significant patterns
#' @export
significant_patterns <- function(x, sort_by = "p_adjusted") {
  if (!inherits(x, "sequence_patterns")) {
    stop("x must be a sequence_patterns object")
  }
  
  if (!"significant" %in% names(x$patterns)) {
    stop("No significance testing was performed")
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

#' Filter patterns by criteria
#' @param x sequence_patterns object
#' @param min_support Minimum support
#' @param min_count Minimum count
#' @param pattern_length Pattern length(s) to include
#' @param significant_only Only return significant patterns
#' @param start_state Filter by starting state
#' @param end_state Filter by ending state
#' @param contains_state Filter by containing state
#' @return Filtered data frame
#' @export
filter_patterns <- function(x, min_support = NULL, min_count = NULL,
                            pattern_length = NULL, significant_only = FALSE,
                            start_state = NULL, end_state = NULL,
                            contains_state = NULL) {
  
  if (!inherits(x, "sequence_patterns")) {
    stop("x must be a sequence_patterns object")
  }
  
  patterns <- x$patterns
  
  # Threshold filters
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
  
  # State filters
  patterns <- filter_by_state(patterns, start_state, end_state, contains_state)
  
  return(patterns)
}

#' Get top sequences/patterns
#' @param x sequence_patterns object
#' @param top_n Number to return
#' @param min_support Minimum support
#' @param significant_only Only significant patterns
#' @return Data frame of top patterns
#' @export
top_sequences <- function(x, top_n = 10, min_support = NULL, significant_only = FALSE) {
  
  if (!inherits(x, "sequence_patterns")) {
    stop("x must be a sequence_patterns object")
  }
  
  patterns <- x$patterns
  
  if (!is.null(min_support)) {
    patterns <- patterns[patterns$support >= min_support, ]
  }
  
  if (significant_only) {
    if (!"significant" %in% names(patterns)) {
      stop("No significance testing was performed")
    }
    patterns <- patterns[patterns$significant, ]
  }
  
  patterns <- patterns[order(patterns$count, decreasing = TRUE), ]
  
  return(head(patterns, top_n))
}

# ==============================================================================
# BACKWARD COMPATIBILITY
# ==============================================================================

#' @rdname explore_patterns
#' @export
explore_sequence_patterns <- function(...) {
  .Deprecated("explore_patterns")
  explore_patterns(...)
}

cat("Sequence pattern toolkit loaded.\n")
cat("Main function: explore_patterns()\n")
cat("Helper functions: significant_patterns(), filter_patterns(), top_sequences()\n")

# ==============================================================================
# SEQUENCE MOTIF ANALYSIS TOOLKIT
# ==============================================================================
#
# Advanced functions for discovering structural patterns, gap-constrained 
# patterns, and meta-paths in sequential data.
#
# Functions:
# - detect_abstract_patterns(): Find returns, oscillations, progressions
# - find_gapped_patterns(): Find patterns with wildcards (A-*-B)
# - find_meta_paths(): Discover meta-paths across node types
#
# ==============================================================================

# ==============================================================================
# ABSTRACT PATTERN DETECTION
# ==============================================================================

#' Detect Abstract Structural Patterns in Sequences
#'
#' Discovers structural patterns in sequential data including returns (A→*→A),
#' repetitions (A→A), oscillations (A→B→A→B), and progressions (unique state chains).
#' These patterns reveal underlying behavioral dynamics in the sequences.
#'
#' @param data Data frame with sequences in wide format (rows are sequences,
#'   columns are time points). Can also be a group_tna object.
#' @param patterns Character vector of pattern types to detect. Options:
#'   \itemize{
#'     \item \code{"returns"} - State returns after gap (A→*→A)
#'     \item \code{"repetitions"} - Consecutive repeats (A→A, A→A→A)
#'     \item \code{"oscillations"} - Alternating states (A→B→A→B)
#'     \item \code{"progressions"} - Unique state chains without repeats
#'     \item \code{"all"} - Detect all pattern types (default)
#'   }
#' @param min_gap Minimum gap size for return patterns (default: 1)
#' @param max_gap Maximum gap size for return patterns (default: 5)
#' @param min_support Minimum support threshold (default: 0.01)
#' @param min_count Minimum count threshold (default: 2)
#' @param test_significance Whether to perform significance tests (default: TRUE)
#' @param correction Multiple testing correction method (default: "fdr")
#' @param alpha Significance level (default: 0.05)
#' @param verbose Whether to print progress messages (default: TRUE)
#'
#' @return An object of class "abstract_patterns" containing:
#'   \item{returns}{Data frame of return patterns with frequencies and significance}
#'   \item{repetitions}{Data frame of repetition patterns}
#'   \item{oscillations}{Data frame of oscillation patterns}
#'   \item{progressions}{Data frame of progression patterns}
#'   \item{summary}{Summary statistics}
#'   \item{parameters}{Analysis parameters}
#'
#' @examples
#' \dontrun{
#' # Load regulation data from tna package
#' data <- tna::group_regulation
#' seq_data <- data[, -1]  # Remove group column
#'
#' # Detect all abstract patterns
#' patterns <- detect_abstract_patterns(seq_data)
#' print(patterns)
#'
#' # Detect specific pattern types
#' returns <- detect_abstract_patterns(seq_data, patterns = "returns", max_gap = 3)
#' oscillations <- detect_abstract_patterns(seq_data, patterns = "oscillations")
#' }
#'
#' @export
detect_abstract_patterns <- function(data,
                                    patterns = "all",
                                    min_gap = 1,
                                    max_gap = 5,
                                    min_support = 0.01,
                                    min_count = 2,
                                    test_significance = TRUE,
                                    correction = "fdr",
                                    alpha = 0.05,
                                    verbose = TRUE) {
  
  # Input validation
  if (is_group_tna(data)) {
    if (verbose) cat("Detected group_tna object, converting...\n")
    converted <- convert_group_tna(data)
    data <- converted$data
    group_col <- converted$group_col
    if (!is.null(group_col) && group_col %in% names(data)) {
      data <- data[, names(data) != group_col, drop = FALSE]
    }
  }
  
  data <- as.data.frame(data)
  
  if (nrow(data) < 2) stop("data must contain at least 2 sequences")
  
  # Handle "all" patterns
  if (length(patterns) == 1 && patterns == "all") {
    patterns <- c("returns", "repetitions", "oscillations", "progressions")
  }
  
  valid_patterns <- c("returns", "repetitions", "oscillations", "progressions")
  if (!all(patterns %in% valid_patterns)) {
    stop("patterns must be one or more of: ", paste(valid_patterns, collapse = ", "))
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
  
  if (n_sequences == 0) stop("No valid sequences found")
  
  all_states <- unique(unlist(sequences))
  n_states <- length(all_states)
  
  if (verbose) {
    cat("  Sequences:", n_sequences, "\n")
    cat("  Unique states:", n_states, "\n")
  }
  
  results <- list()
  
  # Detect returns (A→*→A)
  if ("returns" %in% patterns) {
    if (verbose) cat("Detecting return patterns...\n")
    results$returns <- detect_returns(sequences, all_states, min_gap, max_gap,
                                     n_sequences, min_support, min_count,
                                     test_significance, correction, alpha)
  }
  
  # Detect repetitions (A→A)
  if ("repetitions" %in% patterns) {
    if (verbose) cat("Detecting repetition patterns...\n")
    results$repetitions <- detect_repetitions(sequences, all_states, n_sequences,
                                             min_support, min_count,
                                             test_significance, correction, alpha)
  }
  
  # Detect oscillations (A→B→A→B)
  if ("oscillations" %in% patterns) {
    if (verbose) cat("Detecting oscillation patterns...\n")
    results$oscillations <- detect_oscillations(sequences, all_states, n_sequences,
                                               min_support, min_count,
                                               test_significance, correction, alpha)
  }
  
  # Detect progressions (unique chains)
  if ("progressions" %in% patterns) {
    if (verbose) cat("Detecting progression patterns...\n")
    results$progressions <- detect_progressions(sequences, n_sequences,
                                               min_support, min_count,
                                               test_significance, correction, alpha)
  }
  
  # Summary
  results$summary <- list(
    n_sequences = n_sequences,
    n_states = n_states,
    patterns_detected = patterns,
    n_returns = if (!is.null(results$returns)) nrow(results$returns) else 0,
    n_repetitions = if (!is.null(results$repetitions)) nrow(results$repetitions) else 0,
    n_oscillations = if (!is.null(results$oscillations)) nrow(results$oscillations) else 0,
    n_progressions = if (!is.null(results$progressions)) nrow(results$progressions) else 0
  )
  
  results$parameters <- list(
    min_gap = min_gap,
    max_gap = max_gap,
    min_support = min_support,
    min_count = min_count,
    correction = correction,
    alpha = alpha,
    all_states = all_states
  )
  
  class(results) <- "abstract_patterns"
  
  if (verbose) cat("Abstract pattern detection complete.\n")
  
  return(results)
}

#' Detect return patterns (A→*→A)
#' @keywords internal
detect_returns <- function(sequences, all_states, min_gap, max_gap,
                          n_sequences, min_support, min_count,
                          test_significance, correction, alpha) {
  
  returns_list <- list()
  
  for (state in all_states) {
    for (gap in min_gap:max_gap) {
      pattern_name <- paste0(state, "->(*", gap, ")->", state)
      count <- 0
      seqs_containing <- 0
      
      for (seq in sequences) {
        positions <- which(seq == state)
        if (length(positions) >= 2) {
          found_in_seq <- FALSE
          for (i in 1:(length(positions) - 1)) {
            actual_gap <- positions[i + 1] - positions[i] - 1
            if (actual_gap == gap) {
              count <- count + 1
              found_in_seq <- TRUE
            }
          }
          if (found_in_seq) seqs_containing <- seqs_containing + 1
        }
      }
      
      support <- seqs_containing / n_sequences
      
      if (count >= min_count && support >= min_support) {
        returns_list[[pattern_name]] <- data.frame(
          pattern = pattern_name,
          state = state,
          gap = gap,
          count = count,
          sequences_containing = seqs_containing,
          support = support,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  
  if (length(returns_list) == 0) {
    return(data.frame(
      pattern = character(0), state = character(0), gap = integer(0),
      count = integer(0), sequences_containing = integer(0), support = numeric(0),
      stringsAsFactors = FALSE
    ))
  }
  
  returns_df <- do.call(rbind, returns_list)
  rownames(returns_df) <- NULL
  
  # Significance testing
  if (test_significance && nrow(returns_df) > 0) {
    n_states <- length(all_states)
    returns_df$expected_prob <- (1 / n_states) * (1 / n_states)
    returns_df$p_value <- sapply(1:nrow(returns_df), function(i) {
      stats::binom.test(returns_df$sequences_containing[i], n_sequences,
                       returns_df$expected_prob[i], alternative = "greater")$p.value
    })
    returns_df$p_adjusted <- stats::p.adjust(returns_df$p_value, method = correction)
    returns_df$significant <- returns_df$p_adjusted < alpha
    returns_df$lift <- returns_df$support / returns_df$expected_prob
  }
  
  returns_df <- returns_df[order(returns_df$count, decreasing = TRUE), ]
  rownames(returns_df) <- NULL
  
  return(returns_df)
}

#' Detect repetition patterns (A→A)
#' @keywords internal
detect_repetitions <- function(sequences, all_states, n_sequences,
                              min_support, min_count,
                              test_significance, correction, alpha) {
  
  reps_list <- list()
  
  for (state in all_states) {
    for (length in 2:4) {
      count <- 0
      seqs_containing <- 0
      
      for (seq in sequences) {
        rle_result <- rle(seq)
        matches <- sum(rle_result$values == state & rle_result$lengths >= length)
        if (matches > 0) {
          count <- count + matches
          seqs_containing <- seqs_containing + 1
        }
      }
      
      support <- seqs_containing / n_sequences
      
      if (count >= min_count && support >= min_support) {
        pattern_name <- paste(rep(state, length), collapse = "->")
        reps_list[[paste0(state, "_", length)]] <- data.frame(
          pattern = pattern_name,
          state = state,
          length = length,
          count = count,
          sequences_containing = seqs_containing,
          support = support,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  
  if (length(reps_list) == 0) {
    return(data.frame(
      pattern = character(0), state = character(0), length = integer(0),
      count = integer(0), sequences_containing = integer(0), support = numeric(0),
      stringsAsFactors = FALSE
    ))
  }
  
  reps_df <- do.call(rbind, reps_list)
  rownames(reps_df) <- NULL
  
  if (test_significance && nrow(reps_df) > 0) {
    n_states <- length(all_states)
    reps_df$expected_prob <- (1 / n_states) ^ reps_df$length
    reps_df$p_value <- sapply(1:nrow(reps_df), function(i) {
      stats::binom.test(reps_df$sequences_containing[i], n_sequences,
                       reps_df$expected_prob[i], alternative = "greater")$p.value
    })
    reps_df$p_adjusted <- stats::p.adjust(reps_df$p_value, method = correction)
    reps_df$significant <- reps_df$p_adjusted < alpha
    reps_df$lift <- reps_df$support / reps_df$expected_prob
  }
  
  reps_df <- reps_df[order(reps_df$count, decreasing = TRUE), ]
  rownames(reps_df) <- NULL
  
  return(reps_df)
}

#' Detect oscillation patterns (A→B→A→B)
#' @keywords internal
detect_oscillations <- function(sequences, all_states, n_sequences,
                               min_support, min_count,
                               test_significance, correction, alpha) {
  
  osc_list <- list()
  n_states <- length(all_states)
  
  for (i in 1:(n_states - 1)) {
    for (j in (i + 1):n_states) {
      state_a <- all_states[i]
      state_b <- all_states[j]
      
      count <- 0
      seqs_containing <- 0
      
      for (seq in sequences) {
        # Look for A-B-A-B pattern
        found <- FALSE
        if (length(seq) >= 4) {
          for (k in 1:(length(seq) - 3)) {
            if (seq[k] == state_a && seq[k + 1] == state_b &&
                seq[k + 2] == state_a && seq[k + 3] == state_b) {
              count <- count + 1
              found <- TRUE
            }
          }
        }
        if (found) seqs_containing <- seqs_containing + 1
      }
      
      support <- seqs_containing / n_sequences
      
      if (count >= min_count && support >= min_support) {
        pattern_name <- paste(state_a, state_b, state_a, state_b, sep = "->")
        osc_list[[paste0(state_a, "_", state_b)]] <- data.frame(
          pattern = pattern_name,
          state_a = state_a,
          state_b = state_b,
          count = count,
          sequences_containing = seqs_containing,
          support = support,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  
  if (length(osc_list) == 0) {
    return(data.frame(
      pattern = character(0), state_a = character(0), state_b = character(0),
      count = integer(0), sequences_containing = integer(0), support = numeric(0),
      stringsAsFactors = FALSE
    ))
  }
  
  osc_df <- do.call(rbind, osc_list)
  rownames(osc_df) <- NULL
  
  if (test_significance && nrow(osc_df) > 0) {
    osc_df$expected_prob <- (1 / n_states) ^ 4
    osc_df$p_value <- sapply(1:nrow(osc_df), function(i) {
      stats::binom.test(osc_df$sequences_containing[i], n_sequences,
                       osc_df$expected_prob[i], alternative = "greater")$p.value
    })
    osc_df$p_adjusted <- stats::p.adjust(osc_df$p_value, method = correction)
    osc_df$significant <- osc_df$p_adjusted < alpha
    osc_df$lift <- osc_df$support / osc_df$expected_prob
  }
  
  osc_df <- osc_df[order(osc_df$count, decreasing = TRUE), ]
  rownames(osc_df) <- NULL
  
  return(osc_df)
}

#' Detect progression patterns (unique chains)
#' @keywords internal
detect_progressions <- function(sequences, n_sequences,
                               min_support, min_count,
                               test_significance, correction, alpha) {
  
  progressions <- list()
  
  for (seq in sequences) {
    if (length(seq) < 3) next
    
    # Find maximal progressions (no repeated states)
    for (start in 1:(length(seq) - 2)) {
      seen <- c(seq[start])
      end <- start
      
      for (pos in (start + 1):length(seq)) {
        if (seq[pos] %in% seen) break
        seen <- c(seen, seq[pos])
        end <- pos
      }
      
      if (length(seen) >= 3) {
        prog <- paste(seen, collapse = "->")
        if (is.null(progressions[[prog]])) {
          progressions[[prog]] <- list(count = 0, seqs = c())
        }
        progressions[[prog]]$count <- progressions[[prog]]$count + 1
        progressions[[prog]]$seqs <- unique(c(progressions[[prog]]$seqs, 
                                              which(sapply(sequences, identical, seq))))
      }
    }
  }
  
  if (length(progressions) == 0) {
    return(data.frame(
      pattern = character(0), length = integer(0),
      count = integer(0), sequences_containing = integer(0), support = numeric(0),
      stringsAsFactors = FALSE
    ))
  }
  
  prog_df <- data.frame(
    pattern = names(progressions),
    length = sapply(names(progressions), function(p) length(strsplit(p, "->")[[1]])),
    count = sapply(progressions, function(p) p$count),
    sequences_containing = sapply(progressions, function(p) length(p$seqs)),
    stringsAsFactors = FALSE
  )
  
  prog_df$support <- prog_df$sequences_containing / n_sequences
  
  prog_df <- prog_df[prog_df$count >= min_count & prog_df$support >= min_support, ]
  
  if (nrow(prog_df) == 0) {
    return(data.frame(
      pattern = character(0), length = integer(0),
      count = integer(0), sequences_containing = integer(0), support = numeric(0),
      stringsAsFactors = FALSE
    ))
  }
  
  if (test_significance && nrow(prog_df) > 0) {
    n_states <- length(unique(unlist(sequences)))
    prog_df$expected_prob <- sapply(prog_df$length, function(l) {
      prod((n_states:(n_states - l + 1)) / n_states^l)
    })
    prog_df$p_value <- sapply(1:nrow(prog_df), function(i) {
      stats::binom.test(prog_df$sequences_containing[i], n_sequences,
                       prog_df$expected_prob[i], alternative = "greater")$p.value
    })
    prog_df$p_adjusted <- stats::p.adjust(prog_df$p_value, method = correction)
    prog_df$significant <- prog_df$p_adjusted < alpha
    prog_df$lift <- prog_df$support / prog_df$expected_prob
  }
  
  prog_df <- prog_df[order(prog_df$count, decreasing = TRUE), ]
  rownames(prog_df) <- NULL
  
  return(prog_df)
}

# ==============================================================================
# GAP-CONSTRAINED PATTERN FINDING
# ==============================================================================

#' Find Gap-Constrained Patterns in Sequences
#'
#' Discovers patterns with wildcards/gaps in sequential data. Supports single
#' wildcards (*) matching any single state and multi-wildcards (**) matching
#' any number of states.
#'
#' @param data Data frame with sequences in wide format (rows are sequences,
#'   columns are time points). Can also be a group_tna object.
#' @param pattern Pattern to search for with wildcards. Use:
#'   \itemize{
#'     \item \code{*} for single wildcard (matches exactly one state)
#'     \item \code{**} for multi-wildcard (matches zero or more states)
#'   }
#'   If NULL, discovers frequent gapped patterns automatically.
#' @param min_gap Minimum gap size for discovered patterns (default: 1)
#' @param max_gap Maximum gap size for discovered patterns (default: 3)
#' @param min_support Minimum support threshold (default: 0.01)
#' @param min_count Minimum count threshold (default: 2)
#' @param test_significance Whether to perform significance tests (default: TRUE)
#' @param correction Multiple testing correction method (default: "fdr")
#' @param alpha Significance level (default: 0.05)
#' @param verbose Whether to print progress messages (default: TRUE)
#'
#' @return An object of class "gapped_patterns" containing:
#'   \item{patterns}{Data frame of gapped patterns with frequencies and significance}
#'   \item{instances}{Data frame of actual pattern instances found}
#'   \item{summary}{Summary statistics}
#'   \item{parameters}{Analysis parameters}
#'
#' @examples
#' \dontrun{
#' # Load regulation data
#' data <- tna::group_regulation
#' seq_data <- data[, -1]
#'
#' # Search for specific gapped pattern
#' results <- find_gapped_patterns(seq_data, pattern = "plan-*-consensus")
#'
#' # Multi-gap pattern
#' results <- find_gapped_patterns(seq_data, pattern = "plan-**-plan")
#'
#' # Auto-discover gapped patterns
#' results <- find_gapped_patterns(seq_data, pattern = NULL, max_gap = 2)
#' }
#'
#' @export
find_gapped_patterns <- function(data,
                                 pattern = NULL,
                                 min_gap = 1,
                                 max_gap = 3,
                                 min_support = 0.01,
                                 min_count = 2,
                                 test_significance = TRUE,
                                 correction = "fdr",
                                 alpha = 0.05,
                                 verbose = TRUE) {
  
  # Input validation
  if (is_group_tna(data)) {
    if (verbose) cat("Detected group_tna object, converting...\n")
    converted <- convert_group_tna(data)
    data <- converted$data
    group_col <- converted$group_col
    if (!is.null(group_col) && group_col %in% names(data)) {
      data <- data[, names(data) != group_col, drop = FALSE]
    }
  }
  
  data <- as.data.frame(data)
  
  if (nrow(data) < 2) stop("data must contain at least 2 sequences")
  
  # Convert to sequences
  if (verbose) cat("Processing sequences...\n")
  sequences <- apply(data, 1, function(row) {
    valid <- row[!is.na(row) & row != "" & row != "NA"]
    if (length(valid) == 0) return(NULL)
    as.character(valid)
  })
  sequences <- sequences[!sapply(sequences, is.null)]
  n_sequences <- length(sequences)
  
  all_states <- unique(unlist(sequences))
  n_states <- length(all_states)
  
  if (verbose) {
    cat("  Sequences:", n_sequences, "\n")
    cat("  Unique states:", n_states, "\n")
  }
  
  if (!is.null(pattern)) {
    # Search for specific pattern
    if (verbose) cat("Searching for pattern:", pattern, "\n")
    results <- search_gapped_pattern(sequences, pattern, n_sequences, n_states,
                                    test_significance, correction, alpha)
  } else {
    # Auto-discover gapped patterns
    if (verbose) cat("Auto-discovering gapped patterns...\n")
    results <- discover_gapped_patterns(sequences, all_states, min_gap, max_gap,
                                       n_sequences, n_states, min_support, min_count,
                                       test_significance, correction, alpha, verbose)
  }
  
  results$summary <- list(
    n_sequences = n_sequences,
    n_states = n_states,
    n_patterns = if (!is.null(results$patterns)) nrow(results$patterns) else 0,
    n_significant = if (!is.null(results$patterns) && "significant" %in% names(results$patterns))
      sum(results$patterns$significant) else NA
  )
  
  results$parameters <- list(
    pattern = pattern,
    min_gap = min_gap,
    max_gap = max_gap,
    min_support = min_support,
    min_count = min_count,
    correction = correction,
    alpha = alpha,
    all_states = all_states
  )
  
  class(results) <- "gapped_patterns"
  
  if (verbose) cat("Gap-constrained pattern search complete.\n")
  
  return(results)
}

#' Search for specific gapped pattern
#' @keywords internal
search_gapped_pattern <- function(sequences, pattern, n_sequences, n_states,
                                 test_significance, correction, alpha) {
  
  parts <- strsplit(pattern, "-")[[1]]
  
  count <- 0
  seqs_containing <- 0
  instances <- list()
  
  for (seq_idx in seq_along(sequences)) {
    seq <- sequences[[seq_idx]]
    matches <- find_pattern_matches(seq, parts)
    
    if (length(matches) > 0) {
      count <- count + length(matches)
      seqs_containing <- seqs_containing + 1
      for (m in matches) {
        instances[[length(instances) + 1]] <- data.frame(
          sequence_id = seq_idx,
          instance = paste(m, collapse = "-"),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  
  support <- seqs_containing / n_sequences
  
  patterns_df <- data.frame(
    pattern = pattern,
    count = count,
    sequences_containing = seqs_containing,
    support = support,
    stringsAsFactors = FALSE
  )
  
  if (test_significance) {
    n_fixed <- sum(!parts %in% c("*", "**"))
    expected_prob <- (1 / n_states) ^ n_fixed
    patterns_df$expected_prob <- expected_prob
    patterns_df$p_value <- stats::binom.test(seqs_containing, n_sequences,
                                            expected_prob, alternative = "greater")$p.value
    patterns_df$p_adjusted <- patterns_df$p_value
    patterns_df$significant <- patterns_df$p_adjusted < alpha
    patterns_df$lift <- support / expected_prob
  }
  
  instances_df <- if (length(instances) > 0) do.call(rbind, instances) else
    data.frame(sequence_id = integer(0), instance = character(0), stringsAsFactors = FALSE)
  
  return(list(patterns = patterns_df, instances = instances_df))
}

#' Find pattern matches in a sequence
#' @keywords internal
find_pattern_matches <- function(seq, parts) {
  matches <- list()
  
  if (length(seq) < length(parts)) return(matches)
  
  # Handle ** (multi-wildcard)
  has_multi <- "**" %in% parts
  
  if (has_multi) {
    matches <- find_multi_wildcard_matches(seq, parts)
  } else {
    matches <- find_single_wildcard_matches(seq, parts)
  }
  
  return(matches)
}

#' Find matches with single wildcards only
#' @keywords internal
find_single_wildcard_matches <- function(seq, parts) {
  matches <- list()
  pattern_len <- length(parts)
  
  if (length(seq) < pattern_len) return(matches)
  
  for (start in 1:(length(seq) - pattern_len + 1)) {
    matched <- TRUE
    match_seq <- character(pattern_len)
    
    for (i in 1:pattern_len) {
      pos <- start + i - 1
      if (pos > length(seq) || is.na(seq[pos]) || seq[pos] == "") {
        matched <- FALSE
        break
      }
      
      if (parts[i] == "*") {
        match_seq[i] <- seq[pos]
      } else if (parts[i] == seq[pos]) {
        match_seq[i] <- seq[pos]
      } else {
        matched <- FALSE
        break
      }
    }
    
    if (matched) {
      matches[[length(matches) + 1]] <- match_seq
    }
  }
  
  return(matches)
}

#' Find matches with multi-wildcards
#' @keywords internal
find_multi_wildcard_matches <- function(seq, parts) {
  matches <- list()
  
  # Split by **
  segments <- list()
  current_segment <- c()
  
  for (p in parts) {
    if (p == "**") {
      if (length(current_segment) > 0) {
        segments[[length(segments) + 1]] <- current_segment
        current_segment <- c()
      }
      segments[[length(segments) + 1]] <- "**"
    } else {
      current_segment <- c(current_segment, p)
    }
  }
  if (length(current_segment) > 0) {
    segments[[length(segments) + 1]] <- current_segment
  }
  
  # Simple case: pattern-**-pattern
  if (length(segments) == 3 && identical(segments[[2]], "**")) {
    first_part <- segments[[1]]
    last_part <- segments[[3]]
    
    # Find all positions of first part
    for (start in 1:(length(seq) - length(first_part) - length(last_part))) {
      if (all(seq[start:(start + length(first_part) - 1)] == first_part |
              first_part == "*")) {
        # Find end part
        for (end in (start + length(first_part) + 1):(length(seq) - length(last_part) + 1)) {
          if (all(seq[end:(end + length(last_part) - 1)] == last_part |
                  last_part == "*")) {
            match_seq <- seq[start:(end + length(last_part) - 1)]
            matches[[length(matches) + 1]] <- match_seq
          }
        }
      }
    }
  }
  
  return(matches)
}

#' Discover gapped patterns automatically
#' @keywords internal
discover_gapped_patterns <- function(sequences, all_states, min_gap, max_gap,
                                    n_sequences, n_states, min_support, min_count,
                                    test_significance, correction, alpha, verbose) {
  
  patterns_list <- list()
  instances_list <- list()
  
  # Generate patterns: state-*-state, state-*-*-state, etc.
  for (state_a in all_states) {
    for (state_b in all_states) {
      for (gap in min_gap:max_gap) {
        pattern <- paste0(state_a, "-", paste(rep("*", gap), collapse = "-"), "-", state_b)
        parts <- strsplit(pattern, "-")[[1]]
        
        count <- 0
        seqs_containing <- 0
        
        for (seq_idx in seq_along(sequences)) {
          seq <- sequences[[seq_idx]]
          matches <- find_single_wildcard_matches(seq, parts)
          if (length(matches) > 0) {
            count <- count + length(matches)
            seqs_containing <- seqs_containing + 1
            for (m in matches) {
              instances_list[[length(instances_list) + 1]] <- data.frame(
                pattern = pattern,
                sequence_id = seq_idx,
                instance = paste(m, collapse = "-"),
                stringsAsFactors = FALSE
              )
            }
          }
        }
        
        support <- seqs_containing / n_sequences
        
        if (count >= min_count && support >= min_support) {
          patterns_list[[pattern]] <- data.frame(
            pattern = pattern,
            from = state_a,
            to = state_b,
            gap = gap,
            count = count,
            sequences_containing = seqs_containing,
            support = support,
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
  
  if (length(patterns_list) == 0) {
    return(list(
      patterns = data.frame(
        pattern = character(0), from = character(0), to = character(0),
        gap = integer(0), count = integer(0), sequences_containing = integer(0),
        support = numeric(0), stringsAsFactors = FALSE
      ),
      instances = data.frame(
        pattern = character(0), sequence_id = integer(0), instance = character(0),
        stringsAsFactors = FALSE
      )
    ))
  }
  
  patterns_df <- do.call(rbind, patterns_list)
  rownames(patterns_df) <- NULL
  
  instances_df <- if (length(instances_list) > 0) do.call(rbind, instances_list) else
    data.frame(pattern = character(0), sequence_id = integer(0), 
               instance = character(0), stringsAsFactors = FALSE)
  
  if (test_significance && nrow(patterns_df) > 0) {
    patterns_df$expected_prob <- (1 / n_states) ^ 2
    patterns_df$p_value <- sapply(1:nrow(patterns_df), function(i) {
      stats::binom.test(patterns_df$sequences_containing[i], n_sequences,
                       patterns_df$expected_prob[i], alternative = "greater")$p.value
    })
    patterns_df$p_adjusted <- stats::p.adjust(patterns_df$p_value, method = correction)
    patterns_df$significant <- patterns_df$p_adjusted < alpha
    patterns_df$lift <- patterns_df$support / patterns_df$expected_prob
  }
  
  patterns_df <- patterns_df[order(patterns_df$count, decreasing = TRUE), ]
  rownames(patterns_df) <- NULL
  
  if (verbose) cat("  Patterns found:", nrow(patterns_df), "\n")
  
  return(list(patterns = patterns_df, instances = instances_df))
}

# ==============================================================================
# META-PATH DISCOVERY
# ==============================================================================

#' Find Meta-Paths Across Node Types in Heterogeneous Sequences
#'
#' Discovers meta-paths (type-level patterns) in sequential data where states
#' are grouped into node types. This enables analysis of heterogeneous networks
#' embedded in sequences, revealing patterns like "cognitive→social→cognitive".
#'
#' @param data Data frame with sequences in wide format (rows are sequences,
#'   columns are time points). Can also be a group_tna object.
#' @param node_types Named list mapping type names to state vectors. Example:
#'   \code{list(cognitive = c("plan", "monitor"), social = c("discuss", "consensus"))}
#' @param schema Optional. Specific type schema to search for. Can be specified as:
#'   \itemize{
#'     \item A character vector: \code{c("cognitive", "social", "cognitive")}
#'     \item A string with arrows: \code{"cognitive->social->cognitive"}
#'   }
#'   If NULL (default), auto-discovers all frequent type-paths.
#'   Supports wildcards: "*" (single type) and "**" (multiple types).
#' @param min_length Minimum type-path length for auto-discovery (default: 2)
#' @param max_length Maximum type-path length for auto-discovery (default: 5)
#' @param allow_gaps Allow gaps in type-paths for auto-discovery (default: FALSE)
#' @param max_gap Maximum gap size if gaps allowed (default: 2)
#' @param min_support Minimum support threshold (default: 0.01)
#' @param min_count Minimum count threshold (default: 2)
#' @param test_significance Whether to perform significance tests (default: TRUE)
#' @param correction Multiple testing correction method (default: "fdr")
#' @param alpha Significance level (default: 0.05)
#' @param verbose Whether to print progress messages (default: TRUE)
#'
#' @return An object of class "meta_paths" containing:
#'   \item{meta_paths}{Data frame of type-level patterns with frequencies and significance}
#'   \item{instances}{Data frame of actual state sequences matching each meta-path}
#'   \item{type_transitions}{Data frame of type-to-type transition statistics}
#'   \item{type_mapping}{The node_types mapping used}
#'   \item{summary}{Summary statistics}
#'   \item{parameters}{Analysis parameters}
#'
#' @examples
#' \dontrun{
#' # Load regulation data
#' data <- tna::group_regulation
#' seq_data <- data[, -1]
#'
#' # Define node types for regulation states
#' node_types <- list(
#'   cognitive = c("plan", "monitor", "adapt"),
#'   social = c("discuss", "consensus", "coregulate", "synthesis"),
#'   emotional = c("emotion", "cohesion")
#' )
#'
#' # Auto-discover all meta-paths
#' meta <- find_meta_paths(seq_data, node_types = node_types)
#' print(meta)
#' summary(meta)
#'
#' # Search for specific schema (vector format - recommended)
#' meta_specific <- find_meta_paths(
#'   seq_data, 
#'   node_types = node_types,
#'   schema = c("cognitive", "social", "cognitive")
#' )
#'
#' # Schema as string (also supported)
#' meta_specific <- find_meta_paths(
#'   seq_data, 
#'   node_types = node_types,
#'   schema = "cognitive->social->cognitive"
#' )
#'
#' # Schema with wildcards
#' meta_return <- find_meta_paths(
#'   seq_data,
#'   node_types = node_types, 
#'   schema = c("cognitive", "**", "cognitive")
#' )
#' }
#'
#' @export
find_meta_paths <- function(data,
                           node_types,
                           schema = NULL,
                           min_length = 2,
                           max_length = 5,
                           allow_gaps = FALSE,
                           max_gap = 2,
                           min_support = 0.01,
                           min_count = 2,
                           test_significance = TRUE,
                           correction = "fdr",
                           alpha = 0.05,
                           verbose = TRUE) {
  
  # Input validation
  if (is_group_tna(data)) {
    if (verbose) cat("Detected group_tna object, converting...\n")
    converted <- convert_group_tna(data)
    data <- converted$data
    group_col <- converted$group_col
    if (!is.null(group_col) && group_col %in% names(data)) {
      data <- data[, names(data) != group_col, drop = FALSE]
    }
  }
  
  data <- as.data.frame(data)
  
  if (nrow(data) < 2) stop("data must contain at least 2 sequences")
  
  if (!is.list(node_types) || length(node_types) < 2) {
    stop("node_types must be a named list with at least 2 type definitions")
  }
  
  type_names <- names(node_types)
  if (is.null(type_names) || any(type_names == "")) {
    stop("node_types must have named elements")
  }
  
  # Convert schema to standard format if provided as vector
  schema_original <- schema
  if (!is.null(schema)) {
    if (is.character(schema) && length(schema) > 1) {
      # Vector format: c("type1", "type2", "type3")
      schema <- paste(schema, collapse = "->")
    } else if (is.character(schema) && length(schema) == 1 && !grepl("->", schema)) {
      # Single type name without arrows - treat as single element
      schema <- schema
    }
    # Otherwise assume it's already in "type1->type2" format
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
  
  # Create state-to-type mapping
  state_to_type <- create_state_type_mapping(node_types, all_states)
  
  # Check coverage
  mapped_states <- names(state_to_type)
  unmapped <- setdiff(all_states, mapped_states)
  coverage <- length(mapped_states) / length(all_states)
  
  if (verbose) {
    cat("  Sequences:", n_sequences, "\n")
    cat("  States mapped:", length(mapped_states), "/", length(all_states), 
        sprintf(" (%.1f%%)\n", coverage * 100))
    if (length(unmapped) > 0) {
      cat("  Unmapped states:", paste(head(unmapped, 10), collapse = ", "),
          if (length(unmapped) > 10) paste("... (", length(unmapped) - 10, " more)") else "", "\n")
    }
  }
  
  # Check if any states are mapped
  if (length(mapped_states) == 0) {
    stop("No states in data match the node_types definitions.\n",
         "  States in data: ", paste(head(all_states, 10), collapse = ", "),
         if (length(all_states) > 10) "..." else "", "\n",
         "  States in node_types: ", paste(head(unlist(node_types), 10), collapse = ", "),
         if (length(unlist(node_types)) > 10) "..." else "", "\n",
         "Please check that node_types contains states that exist in your data.")
  }
  
  # Convert sequences to type sequences
  type_sequences <- lapply(sequences, function(seq) {
    types <- sapply(seq, function(s) {
      if (s %in% names(state_to_type)) state_to_type[[s]] else NA
    })
    types[!is.na(types)]
  })
  
  # Count sequences with valid type mappings
  n_with_types <- sum(sapply(type_sequences, length) > 0)
  
  if (n_with_types == 0) {
    stop("No states in any sequence match the node_types definitions.\n",
         "  Sample states from data: ", paste(head(all_states, 10), collapse = ", "), "\n",
         "  States defined in node_types: ", paste(head(unlist(node_types), 10), collapse = ", "), "\n",
         "Please ensure node_types contains states that exist in your data.")
  }
  
  # Filter to sequences with minimum length
  type_sequences_filtered <- type_sequences[sapply(type_sequences, length) >= min_length]
  
  if (length(type_sequences_filtered) == 0) {
    warning("No sequences have ", min_length, " or more consecutive mapped states.\n",
            "  Sequences with at least 1 mapped state: ", n_with_types, "\n",
            "  Try reducing min_length or mapping more states.")
    
    # Use sequences with at least 1 mapped state
    type_sequences_filtered <- type_sequences[sapply(type_sequences, length) >= 1]
    
    if (length(type_sequences_filtered) == 0) {
      stop("No valid type sequences after mapping")
    }
  }
  
  type_sequences <- type_sequences_filtered
  
  if (!is.null(schema)) {
    # Search for specific schema
    if (verbose) cat("Searching for schema:", schema, "\n")
    results <- search_meta_path_schema(sequences, type_sequences, schema,
                                      state_to_type, node_types, n_sequences,
                                      test_significance, correction, alpha)
  } else {
    # Auto-discover meta-paths
    if (verbose) cat("Auto-discovering meta-paths...\n")
    results <- discover_meta_paths(sequences, type_sequences, type_names,
                                  state_to_type, min_length, max_length,
                                  n_sequences, min_support, min_count,
                                  test_significance, correction, alpha, verbose)
  }
  
  # Compute type transitions
  type_transitions <- compute_type_transitions(type_sequences, type_names,
                                              n_sequences, test_significance,
                                              correction, alpha)
  
  results$type_transitions <- type_transitions
  results$type_mapping <- node_types
  
  results$summary <- list(
    n_sequences = n_sequences,
    n_states = length(all_states),
    n_types = length(type_names),
    type_names = type_names,
    coverage = coverage,
    n_meta_paths = if (!is.null(results$meta_paths)) nrow(results$meta_paths) else 0,
    n_significant = if (!is.null(results$meta_paths) && "significant" %in% names(results$meta_paths))
      sum(results$meta_paths$significant) else NA
  )
  
  results$parameters <- list(
    schema = schema,
    min_length = min_length,
    max_length = max_length,
    allow_gaps = allow_gaps,
    max_gap = max_gap,
    min_support = min_support,
    min_count = min_count,
    correction = correction,
    alpha = alpha
  )
  
  class(results) <- "meta_paths"
  
  if (verbose) cat("Meta-path discovery complete.\n")
  
  return(results)
}

#' Create state-to-type mapping
#' @keywords internal
create_state_type_mapping <- function(node_types, all_states) {
  mapping <- list()
  
  for (type_name in names(node_types)) {
    for (state in node_types[[type_name]]) {
      if (state %in% all_states) {
        mapping[[state]] <- type_name
      }
    }
  }
  
  return(mapping)
}

#' Search for specific meta-path schema
#' @keywords internal
search_meta_path_schema <- function(sequences, type_sequences, schema,
                                   state_to_type, node_types, n_sequences,
                                   test_significance, correction, alpha) {
  
  schema_parts <- strsplit(gsub(" ", "", schema), "->")[[1]]
  n_types <- length(names(node_types))
  
  count <- 0
  seqs_containing <- 0
  instances <- list()
  
  has_multi <- "**" %in% schema_parts
  
  for (seq_idx in seq_along(sequences)) {
    seq <- sequences[[seq_idx]]
    type_seq <- type_sequences[[seq_idx]]
    
    if (length(type_seq) < length(schema_parts[schema_parts != "**"])) next
    
    if (has_multi) {
      matches <- find_meta_path_multi_matches(seq, type_seq, schema_parts, state_to_type)
    } else {
      matches <- find_meta_path_matches(seq, type_seq, schema_parts)
    }
    
    if (length(matches) > 0) {
      count <- count + length(matches)
      seqs_containing <- seqs_containing + 1
      for (m in matches) {
        instances[[length(instances) + 1]] <- data.frame(
          schema = schema,
          sequence_id = seq_idx,
          instance = paste(m, collapse = "->"),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  
  support <- seqs_containing / n_sequences
  
  meta_paths_df <- data.frame(
    schema = schema,
    count = count,
    sequences_containing = seqs_containing,
    support = support,
    stringsAsFactors = FALSE
  )
  
  if (test_significance) {
    n_fixed <- sum(!schema_parts %in% c("*", "**"))
    expected_prob <- (1 / n_types) ^ n_fixed
    meta_paths_df$expected_prob <- expected_prob
    meta_paths_df$p_value <- stats::binom.test(seqs_containing, n_sequences,
                                              expected_prob, alternative = "greater")$p.value
    meta_paths_df$p_adjusted <- meta_paths_df$p_value
    meta_paths_df$significant <- meta_paths_df$p_adjusted < alpha
    meta_paths_df$lift <- support / expected_prob
  }
  
  instances_df <- if (length(instances) > 0) do.call(rbind, instances) else
    data.frame(schema = character(0), sequence_id = integer(0), 
               instance = character(0), stringsAsFactors = FALSE)
  
  return(list(meta_paths = meta_paths_df, instances = instances_df))
}

#' Find meta-path matches (no multi-wildcards)
#' @keywords internal
find_meta_path_matches <- function(seq, type_seq, schema_parts) {
  matches <- list()
  pattern_len <- length(schema_parts)
  
  if (length(type_seq) < pattern_len) return(matches)
  
  for (start in 1:(length(type_seq) - pattern_len + 1)) {
    matched <- TRUE
    match_seq <- character(pattern_len)
    
    for (i in 1:pattern_len) {
      if (schema_parts[i] == "*") {
        match_seq[i] <- seq[start + i - 1]
      } else if (schema_parts[i] == type_seq[start + i - 1]) {
        match_seq[i] <- seq[start + i - 1]
      } else {
        matched <- FALSE
        break
      }
    }
    
    if (matched) {
      matches[[length(matches) + 1]] <- match_seq
    }
  }
  
  return(matches)
}

#' Find meta-path matches with multi-wildcards
#' @keywords internal
find_meta_path_multi_matches <- function(seq, type_seq, schema_parts, state_to_type) {
  matches <- list()
  
  # Simple case: type-**-type
  multi_pos <- which(schema_parts == "**")
  
  if (length(multi_pos) == 1 && length(schema_parts) == 3) {
    first_type <- schema_parts[1]
    last_type <- schema_parts[3]
    
    first_positions <- which(type_seq == first_type)
    last_positions <- which(type_seq == last_type)
    
    for (fp in first_positions) {
      for (lp in last_positions) {
        if (lp > fp + 1) {
          match_seq <- seq[fp:lp]
          matches[[length(matches) + 1]] <- match_seq
        }
      }
    }
  }
  
  return(matches)
}

#' Discover meta-paths automatically
#' @keywords internal
discover_meta_paths <- function(sequences, type_sequences, type_names,
                               state_to_type, min_length, max_length,
                               n_sequences, min_support, min_count,
                               test_significance, correction, alpha, verbose) {
  
  n_types <- length(type_names)
  meta_paths_list <- list()
  instances_list <- list()
  
  # Extract type n-grams
  for (len in min_length:max_length) {
    if (verbose) cat("  Extracting type", len, "-grams...\n")
    
    type_ngrams <- list()
    
    for (seq_idx in seq_along(type_sequences)) {
      type_seq <- type_sequences[[seq_idx]]
      seq <- sequences[[seq_idx]]
      
      if (length(type_seq) < len) next
      
      for (start in 1:(length(type_seq) - len + 1)) {
        type_gram <- paste(type_seq[start:(start + len - 1)], collapse = "->")
        state_instance <- paste(seq[start:(start + len - 1)], collapse = "->")
        
        if (is.null(type_ngrams[[type_gram]])) {
          type_ngrams[[type_gram]] <- list(count = 0, seqs = c(), instances = c())
        }
        type_ngrams[[type_gram]]$count <- type_ngrams[[type_gram]]$count + 1
        type_ngrams[[type_gram]]$seqs <- c(type_ngrams[[type_gram]]$seqs, seq_idx)
        type_ngrams[[type_gram]]$instances <- c(type_ngrams[[type_gram]]$instances, state_instance)
      }
    }
    
    for (schema in names(type_ngrams)) {
      info <- type_ngrams[[schema]]
      unique_seqs <- length(unique(info$seqs))
      support <- unique_seqs / n_sequences
      
      if (info$count >= min_count && support >= min_support) {
        meta_paths_list[[schema]] <- data.frame(
          schema = schema,
          length = len,
          count = info$count,
          sequences_containing = unique_seqs,
          support = support,
          stringsAsFactors = FALSE
        )
        
        # Store top instances
        inst_table <- sort(table(info$instances), decreasing = TRUE)
        top_instances <- head(names(inst_table), 5)
        for (inst in top_instances) {
          instances_list[[length(instances_list) + 1]] <- data.frame(
            schema = schema,
            instance = inst,
            count = as.integer(inst_table[inst]),
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
  
  if (length(meta_paths_list) == 0) {
    return(list(
      meta_paths = data.frame(
        schema = character(0), length = integer(0), count = integer(0),
        sequences_containing = integer(0), support = numeric(0),
        stringsAsFactors = FALSE
      ),
      instances = data.frame(
        schema = character(0), instance = character(0), count = integer(0),
        stringsAsFactors = FALSE
      )
    ))
  }
  
  meta_paths_df <- do.call(rbind, meta_paths_list)
  rownames(meta_paths_df) <- NULL
  
  instances_df <- if (length(instances_list) > 0) do.call(rbind, instances_list) else
    data.frame(schema = character(0), instance = character(0), 
               count = integer(0), stringsAsFactors = FALSE)
  
  if (test_significance && nrow(meta_paths_df) > 0) {
    meta_paths_df$expected_prob <- (1 / n_types) ^ meta_paths_df$length
    meta_paths_df$p_value <- sapply(1:nrow(meta_paths_df), function(i) {
      stats::binom.test(meta_paths_df$sequences_containing[i], n_sequences,
                       meta_paths_df$expected_prob[i], alternative = "greater")$p.value
    })
    meta_paths_df$p_adjusted <- stats::p.adjust(meta_paths_df$p_value, method = correction)
    meta_paths_df$significant <- meta_paths_df$p_adjusted < alpha
    meta_paths_df$lift <- meta_paths_df$support / meta_paths_df$expected_prob
  }
  
  meta_paths_df <- meta_paths_df[order(meta_paths_df$count, decreasing = TRUE), ]
  rownames(meta_paths_df) <- NULL
  
  if (verbose) cat("  Meta-paths found:", nrow(meta_paths_df), "\n")
  
  return(list(meta_paths = meta_paths_df, instances = instances_df))
}

#' Compute type-to-type transition statistics
#' @keywords internal
compute_type_transitions <- function(type_sequences, type_names, n_sequences,
                                    test_significance, correction, alpha) {
  
  n_types <- length(type_names)
  trans_counts <- matrix(0, nrow = n_types, ncol = n_types,
                        dimnames = list(type_names, type_names))
  
  for (type_seq in type_sequences) {
    if (length(type_seq) < 2) next
    for (i in 1:(length(type_seq) - 1)) {
      from <- type_seq[i]
      to <- type_seq[i + 1]
      if (from %in% type_names && to %in% type_names) {
        trans_counts[from, to] <- trans_counts[from, to] + 1
      }
    }
  }
  
  # Create data frame
  trans_list <- list()
  for (from in type_names) {
    row_sum <- sum(trans_counts[from, ])
    for (to in type_names) {
      count <- trans_counts[from, to]
      if (count > 0) {
        prob <- count / row_sum
        trans_list[[paste0(from, "_", to)]] <- data.frame(
          from_type = from,
          to_type = to,
          count = count,
          probability = prob,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  
  if (length(trans_list) == 0) {
    return(data.frame(
      from_type = character(0), to_type = character(0),
      count = integer(0), probability = numeric(0),
      stringsAsFactors = FALSE
    ))
  }
  
  trans_df <- do.call(rbind, trans_list)
  rownames(trans_df) <- NULL
  
  if (test_significance && nrow(trans_df) > 0) {
    total_trans <- sum(trans_counts)
    trans_df$expected <- (1 / n_types)
    trans_df$lift <- trans_df$probability / trans_df$expected
    
    trans_df$chi_sq <- sapply(1:nrow(trans_df), function(i) {
      observed <- trans_df$count[i]
      expected <- sum(trans_counts[trans_df$from_type[i], ]) / n_types
      if (expected > 0) (observed - expected)^2 / expected else 0
    })
    
    trans_df$p_value <- sapply(trans_df$chi_sq, function(x) {
      1 - stats::pchisq(x, df = 1)
    })
    trans_df$p_adjusted <- stats::p.adjust(trans_df$p_value, method = correction)
    trans_df$significant <- trans_df$p_adjusted < alpha
  }
  
  trans_df <- trans_df[order(trans_df$count, decreasing = TRUE), ]
  rownames(trans_df) <- NULL
  
  return(trans_df)
}

# ==============================================================================
# S3 METHODS
# ==============================================================================

#' @export
print.abstract_patterns <- function(x, ...) {
  cat("Abstract Pattern Analysis Results\n")
  cat("=================================\n\n")
  
  cat("Data Summary:\n")
  cat("  Sequences:", x$summary$n_sequences, "\n")
  cat("  States:", x$summary$n_states, "\n")
  cat("  Patterns detected:", paste(x$summary$patterns_detected, collapse = ", "), "\n\n")
  
  if (x$summary$n_returns > 0) {
    cat("Return Patterns (A->*->A):", x$summary$n_returns, "\n")
    top <- head(x$returns[order(x$returns$count, decreasing = TRUE), ], 5)
    for (i in 1:nrow(top)) {
      sig <- if ("significant" %in% names(top) && top$significant[i]) " ***" else ""
      cat(sprintf("  %s (count=%d, support=%.3f)%s\n",
                 top$pattern[i], top$count[i], top$support[i], sig))
    }
    cat("\n")
  }
  
  if (x$summary$n_repetitions > 0) {
    cat("Repetition Patterns (A->A):", x$summary$n_repetitions, "\n")
    top <- head(x$repetitions[order(x$repetitions$count, decreasing = TRUE), ], 5)
    for (i in 1:nrow(top)) {
      sig <- if ("significant" %in% names(top) && top$significant[i]) " ***" else ""
      cat(sprintf("  %s (count=%d, support=%.3f)%s\n",
                 top$pattern[i], top$count[i], top$support[i], sig))
    }
    cat("\n")
  }
  
  if (x$summary$n_oscillations > 0) {
    cat("Oscillation Patterns (A->B->A->B):", x$summary$n_oscillations, "\n")
    top <- head(x$oscillations[order(x$oscillations$count, decreasing = TRUE), ], 5)
    for (i in 1:nrow(top)) {
      sig <- if ("significant" %in% names(top) && top$significant[i]) " ***" else ""
      cat(sprintf("  %s (count=%d, support=%.3f)%s\n",
                 top$pattern[i], top$count[i], top$support[i], sig))
    }
    cat("\n")
  }
  
  if (x$summary$n_progressions > 0) {
    cat("Progression Patterns (unique chains):", x$summary$n_progressions, "\n")
    top <- head(x$progressions[order(x$progressions$count, decreasing = TRUE), ], 5)
    for (i in 1:nrow(top)) {
      sig <- if ("significant" %in% names(top) && top$significant[i]) " ***" else ""
      cat(sprintf("  %s (count=%d, support=%.3f)%s\n",
                 top$pattern[i], top$count[i], top$support[i], sig))
    }
  }
  
  cat("\n*** = statistically significant\n")
  invisible(x)
}

#' @export
print.gapped_patterns <- function(x, ...) {
  cat("Gap-Constrained Pattern Analysis Results\n")
  cat("=========================================\n\n")
  
  cat("Data Summary:\n")
  cat("  Sequences:", x$summary$n_sequences, "\n")
  cat("  States:", x$summary$n_states, "\n")
  cat("  Patterns found:", x$summary$n_patterns, "\n")
  if (!is.na(x$summary$n_significant)) {
    cat("  Significant:", x$summary$n_significant, "\n")
  }
  cat("\n")
  
  if (!is.null(x$patterns) && nrow(x$patterns) > 0) {
    cat("Top Gapped Patterns:\n")
    top <- head(x$patterns[order(x$patterns$count, decreasing = TRUE), ], 10)
    for (i in 1:nrow(top)) {
      sig <- if ("significant" %in% names(top) && top$significant[i]) " ***" else ""
      cat(sprintf("  %s (count=%d, support=%.3f)%s\n",
                 top$pattern[i], top$count[i], top$support[i], sig))
    }
  }
  
  cat("\n*** = statistically significant\n")
  invisible(x)
}

#' @export
print.meta_paths <- function(x, ...) {
  cat("Meta-Path Analysis Results\n")
  cat("==========================\n\n")
  
  cat("Node Types Defined:\n")
  for (type_name in names(x$type_mapping)) {
    states <- x$type_mapping[[type_name]]
    cat(sprintf("  %s: %s (%d states)\n", type_name,
               paste(head(states, 5), collapse = ", "),
               length(states)))
  }
  cat("\n")
  
  cat("Data Summary:\n")
  cat("  Sequences:", x$summary$n_sequences, "\n")
  cat("  Type coverage:", sprintf("%.1f%%", x$summary$coverage * 100), "\n")
  cat("  Meta-paths found:", x$summary$n_meta_paths, "\n")
  if (!is.na(x$summary$n_significant)) {
    cat("  Significant:", x$summary$n_significant, "\n")
  }
  cat("\n")
  
  if (!is.null(x$meta_paths) && nrow(x$meta_paths) > 0) {
    cat("Top Meta-Paths by Frequency:\n")
    top <- head(x$meta_paths[order(x$meta_paths$count, decreasing = TRUE), ], 10)
    for (i in 1:nrow(top)) {
      sig <- if ("significant" %in% names(top) && top$significant[i]) " ***" else ""
      cat(sprintf("  %s\n", top$schema[i]))
      cat(sprintf("    count=%d, support=%.3f, lift=%.1f%s\n",
                 top$count[i], top$support[i],
                 if ("lift" %in% names(top)) top$lift[i] else NA, sig))
    }
  }
  
  cat("\n*** = statistically significant (after FDR correction)\n")
  invisible(x)
}

#' @export
summary.meta_paths <- function(object, ...) {
  cat("Meta-Path Analysis - Detailed Summary\n")
  cat("=====================================\n\n")
  
  cat("NODE TYPE DEFINITIONS\n")
  cat(paste(rep("-", 40), collapse = ""), "\n")
  for (type_name in names(object$type_mapping)) {
    states <- object$type_mapping[[type_name]]
    cat(sprintf("%s:\n  %s\n", type_name, paste(states, collapse = ", ")))
  }
  cat("\n")
  
  cat("TYPE TRANSITIONS\n")
  cat(paste(rep("-", 40), collapse = ""), "\n")
  if (!is.null(object$type_transitions) && nrow(object$type_transitions) > 0) {
    top_trans <- head(object$type_transitions, 10)
    for (i in 1:nrow(top_trans)) {
      sig <- if ("significant" %in% names(top_trans) && top_trans$significant[i]) " ***" else ""
      cat(sprintf("  %s -> %s: count=%d, prob=%.3f, lift=%.2f%s\n",
                 top_trans$from_type[i], top_trans$to_type[i],
                 top_trans$count[i], top_trans$probability[i],
                 if ("lift" %in% names(top_trans)) top_trans$lift[i] else NA, sig))
    }
  }
  cat("\n")
  
  if (!is.null(object$instances) && nrow(object$instances) > 0) {
    cat("TOP INSTANCES BY SCHEMA\n")
    cat(paste(rep("-", 40), collapse = ""), "\n")
    
    schemas <- unique(object$instances$schema)
    for (schema in head(schemas, 5)) {
      cat(sprintf("%s:\n", schema))
      schema_instances <- object$instances[object$instances$schema == schema, ]
      schema_instances <- schema_instances[order(schema_instances$count, decreasing = TRUE), ]
      for (j in 1:min(3, nrow(schema_instances))) {
        cat(sprintf("    %s (count=%d)\n", 
                   schema_instances$instance[j], schema_instances$count[j]))
      }
    }
  }
  
  invisible(object)
}

cat("Sequence motif analysis toolkit loaded.\n")
cat("Functions: detect_abstract_patterns(), find_gapped_patterns(), find_meta_paths()\n")


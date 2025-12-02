# ==============================================================================
# SEQUENCE MOTIF ANALYSIS TOOLKIT
# ==============================================================================
#
# Functions for discovering gap-constrained patterns and meta-paths
# in sequential data.
#
# Main functions:
# - find_patterns(): Search for patterns with wildcards/gaps
# - find_meta_paths(): Discover meta-paths across node types
#
# ==============================================================================

# ==============================================================================
# GAP-CONSTRAINED PATTERN FINDING
# ==============================================================================

#' Find Patterns with Wildcards in Sequences
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
#'   Can be specified as vector: \code{c("A", "*", "B")} or string: \code{"A->*->B"}
#'   If NULL, discovers frequent gapped patterns automatically.
#' @param min_gap Minimum gap size for discovered patterns (default: 1)
#' @param max_gap Maximum gap size for discovered patterns (default: 3)
#' @param min_support Minimum support threshold (default: 0.01)
#' @param min_count Minimum count threshold (default: 2)
#' @param test_significance Whether to perform significance tests (default: TRUE)
#' @param correction Multiple testing correction method (default: "fdr")
#' @param alpha Significance level (default: 0.05)
#' @param start_state Filter patterns starting with this state (default: NULL)
#' @param end_state Filter patterns ending with this state (default: NULL)
#' @param contains_state Filter patterns containing this state (default: NULL)
#' @param verbose Whether to print progress messages (default: TRUE)
#'
#' @return An object of class "sequence_patterns" containing:
#'   \item{patterns}{Data frame of patterns with statistics}
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
#' results <- find_patterns(seq_data, pattern = c("plan", "*", "consensus"))
#'
#' # Multi-gap pattern
#' results <- find_patterns(seq_data, pattern = c("plan", "**", "plan"))
#'
#' # Auto-discover gapped patterns
#' results <- find_patterns(seq_data, pattern = NULL, max_gap = 2)
#' }
#'
#' @export
find_patterns <- function(data,
                          pattern = NULL,
                          min_gap = 1,
                          max_gap = 3,
                          min_support = 0.01,
                          min_count = 2,
                          test_significance = TRUE,
                          correction = "fdr",
                          alpha = 0.05,
                          start_state = NULL,
                          end_state = NULL,
                          contains_state = NULL,
                          verbose = TRUE) {
  
  validate_sequence_params(list(min_support = min_support, min_count = min_count, correction = correction))
  
  # Prepare data
  seq_data <- prepare_sequence_input(data, verbose)
  sequences <- seq_data$sequences
  n_sequences <- seq_data$n_sequences
  n_states <- seq_data$n_states
  
  # Normalize pattern string
  if (!is.null(pattern) && length(pattern) > 1) pattern <- paste(pattern, collapse = "->")
  
  if (!is.null(pattern)) {
    if (verbose) cat("Searching for pattern:", pattern, "\n")
    results <- search_gapped_pattern(sequences, pattern, n_sequences, n_states,
                                     test_significance, correction, alpha)
  } else {
    if (verbose) cat("Auto-discovering gapped patterns...\n")
    results <- discover_gapped_patterns(sequences, seq_data$all_states, min_gap, max_gap,
                                        n_sequences, n_states, min_support, min_count,
                                        test_significance, correction, alpha, verbose)
  }
  
  # Apply filters
  patterns_df <- filter_patterns_by_text(results$patterns, start_state, end_state, contains_state)
  instances_df <- filter_patterns_by_text(results$instances, start_state, end_state, contains_state, text_col = "instance")
  
  result <- list(
    patterns = standardize_columns(patterns_df),
    instances = standardize_columns(instances_df),
    summary = list(
      n_sequences = n_sequences,
      n_states = n_states,
      n_patterns = nrow(patterns_df),
      n_significant = if ("significant" %in% names(patterns_df)) sum(patterns_df$significant, na.rm = TRUE) else NA
    ),
    parameters = list(
      pattern = pattern, min_gap = min_gap, max_gap = max_gap,
      min_support = min_support, min_count = min_count,
      test_significance = test_significance, correction = correction, alpha = alpha,
      start_state = start_state, end_state = end_state, contains_state = contains_state
    )
  )
  
  class(result) <- "sequence_patterns"
  if (verbose) cat("Pattern search complete.\n")
  return(result)
}

#' Search for specific gapped pattern
#' @keywords internal
search_gapped_pattern <- function(sequences, pattern, n_sequences, n_states,
                                  test_significance, correction, alpha) {
  
  parts <- strsplit(gsub(" ", "", pattern), "->")[[1]]
  count <- 0; seqs_containing <- 0
  instance_list <- c(); instance_seqs <- list()
  
  for (seq_idx in seq_along(sequences)) {
    matches <- find_pattern_matches(sequences[[seq_idx]], parts)
    
    if (length(matches) > 0) {
      for (m in matches) {
        inst_str <- paste(m, collapse = "->")
        if (!grepl("NA", inst_str)) {
          count <- count + 1
          instance_list <- c(instance_list, inst_str)
          if (is.null(instance_seqs[[inst_str]])) instance_seqs[[inst_str]] <- c()
          instance_seqs[[inst_str]] <- unique(c(instance_seqs[[inst_str]], seq_idx))
        }
      }
      seqs_containing <- seqs_containing + 1
    }
  }
  
  support <- seqs_containing / n_sequences
  patterns_df <- data.frame(
    pattern = pattern,
    length = sum(!parts %in% c("*", "**")),
    count = count,
    sequences_containing = seqs_containing,
    support = round(support, 4),
    stringsAsFactors = FALSE
  )
  
  if (test_significance && count > 0) {
    patterns_df <- compute_significance_stats(patterns_df, n_sequences, n_states, correction, alpha)
  }
  
  # Build instances using unified aggregation
  instances_df <- aggregate_instances(instance_list, instance_seqs, n_sequences, n_states, 
                                      schema = NULL, test_significance, correction, alpha)
  # Add pattern column to instances
  if (nrow(instances_df) > 0) instances_df$pattern <- pattern
  
  return(list(patterns = patterns_df, instances = instances_df))
}

#' Unified pattern match finder
#' @keywords internal
find_pattern_matches <- function(seq, parts) {
  matches <- list()
  
  # Multi-wildcard case: A->**->B
  if ("**" %in% parts) {
    multi_pos <- which(parts == "**")
    if (length(multi_pos) == 1 && length(parts) == 3) {
      starts <- which(seq == parts[1])
      ends <- which(seq == parts[3])
      for (s in starts) {
        for (e in ends) {
          if (e > s + 1) matches[[length(matches) + 1]] <- seq[s:e]
        }
      }
    }
    return(matches)
  }
  
  # Single wildcard or exact match
  len <- length(parts)
  if (length(seq) < len) return(matches)
  
  for (start in 1:(length(seq) - len + 1)) {
    sub_seq <- seq[start:(start + len - 1)]
    matched <- TRUE
    for (i in 1:len) {
      if (parts[i] != "*" && parts[i] != sub_seq[i]) {
        matched <- FALSE; break
      }
    }
    if (matched) matches[[length(matches) + 1]] <- sub_seq
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
  
  for (s1 in all_states) {
    for (s2 in all_states) {
      for (gap in min_gap:max_gap) {
        pattern_name <- paste0(s1, "->", paste(rep("*", gap), collapse = "->"), "->", s2)
        count <- 0; seqs_containing <- 0
        pat_instances <- c(); inst_seqs <- list()
        
        for (seq_idx in seq_along(sequences)) {
          seq <- sequences[[seq_idx]]
          if (length(seq) < gap + 2) next
          
          seq_count <- 0
          for (i in 1:(length(seq) - gap - 1)) {
            if (seq[i] == s1 && seq[i + gap + 1] == s2) {
              seq_count <- seq_count + 1
              inst <- paste(seq[i:(i + gap + 1)], collapse = "->")
              pat_instances <- c(pat_instances, inst)
              if (is.null(inst_seqs[[inst]])) inst_seqs[[inst]] <- c()
              inst_seqs[[inst]] <- unique(c(inst_seqs[[inst]], seq_idx))
            }
          }
          if (seq_count > 0) { count <- count + seq_count; seqs_containing <- seqs_containing + 1 }
        }
        
        if (count >= min_count && seqs_containing / n_sequences >= min_support) {
          patterns_list[[length(patterns_list) + 1]] <- list(
            pattern = pattern_name, length = gap + 2,
            count = count, sequences_containing = seqs_containing
          )
          
          # Aggregate instances for this pattern
          inst_df <- aggregate_instances(pat_instances, inst_seqs, n_sequences, n_states,
                                         schema = NULL, FALSE) # No sig test for individual instances here to save time
          inst_df$pattern <- pattern_name
          instances_list[[length(instances_list) + 1]] <- inst_df
        }
      }
    }
  }
  
  if (length(patterns_list) == 0) {
    return(list(patterns = create_empty_patterns_df(), instances = create_empty_patterns_df()))
  }
  
  # Build patterns DF
  patterns_df <- do.call(rbind, lapply(patterns_list, as.data.frame))
  patterns_df$support <- round(patterns_df$sequences_containing / n_sequences, 4)
  
  if (test_significance) {
    patterns_df <- compute_significance_stats(patterns_df, n_sequences, n_states, correction, alpha)
  }
  
  # Build instances DF
  instances_df <- if (length(instances_list) > 0) do.call(rbind, instances_list) else create_empty_patterns_df()
  
  return(list(
    patterns = patterns_df[order(patterns_df$count, decreasing = TRUE), ],
    instances = instances_df[order(instances_df$count, decreasing = TRUE), ]
  ))
}

# ==============================================================================
# META-PATH DISCOVERY
# ==============================================================================

#' Find Meta-Paths Across Node Types
#'
#' Discovers meta-paths (type-level patterns) in sequential data where states
#' are grouped into node types.
#'
#' @export
find_meta_paths <- function(data,
                            node_types,
                            schema = NULL,
                            min_length = 2,
                            max_length = 5,
                            min_support = 0.01,
                            min_count = 2,
                            test_significance = TRUE,
                            correction = "fdr",
                            alpha = 0.05,
                            start_state = NULL,
                            end_state = NULL,
                            contains_state = NULL,
                            start_schema = NULL,
                            end_schema = NULL,
                            contains_schema = NULL,
                            verbose = TRUE) {
  
  validate_sequence_params(list(min_length = min_length, max_length = max_length,
                                min_support = min_support, min_count = min_count, correction = correction))
  
  if (!is.list(node_types) || length(node_types) < 2) stop("node_types must be a named list with at least 2 types")
  
  seq_data <- prepare_sequence_input(data, verbose)
  
  # Normalize schema
  if (!is.null(schema) && length(schema) > 1) schema <- paste(schema, collapse = "->")
  
  # Create mappings
  state_to_type <- list()
  for (type in names(node_types)) {
    for (state in node_types[[type]]) state_to_type[[state]] <- type
  }
  
  # Convert to type sequences
  type_sequences <- lapply(seq_data$sequences, function(seq) {
    types <- sapply(seq, function(s) if (s %in% names(state_to_type)) state_to_type[[s]] else NA)
    valid <- !is.na(types)
    list(types = types[valid], states = seq[valid])
  })
  
  # Filter short sequences
  type_sequences <- type_sequences[sapply(type_sequences, function(x) length(x$types)) >= min_length]
  if (length(type_sequences) == 0) stop("No valid type sequences after mapping")
  
  if (!is.null(schema)) {
    if (verbose) cat("Searching for schema:", schema, "\n")
    results <- search_meta_path_schema(seq_data$sequences, type_sequences, schema,
                                       state_to_type, node_types, seq_data$n_sequences,
                                       test_significance, correction, alpha)
  } else {
    if (verbose) cat("Auto-discovering meta-paths...\n")
    results <- discover_meta_paths(seq_data$sequences, type_sequences, names(node_types),
                                   state_to_type, min_length, max_length,
                                   seq_data$n_sequences, min_support, min_count,
                                   test_significance, correction, alpha, verbose)
  }
  
  patterns_df <- results$instances
  schemas_df <- results$meta_paths
  
  # Apply filters
  patterns_df <- filter_patterns_by_text(patterns_df, start_state, end_state, contains_state)
  schemas_df <- filter_patterns_by_text(schemas_df, start_schema, end_schema, contains_schema, text_col = "schema")
  
  type_transitions <- compute_type_transitions(type_sequences, names(node_types),
                                               seq_data$n_sequences, test_significance,
                                               correction, alpha)
  
  result <- list(
    patterns = standardize_columns(patterns_df, TRUE),
    schemas = standardize_columns(schemas_df, TRUE),
    type_transitions = type_transitions,
    type_mapping = node_types,
    summary = list(
      n_sequences = seq_data$n_sequences,
      n_patterns = nrow(patterns_df),
      n_schemas = nrow(schemas_df),
      n_significant = if ("significant" %in% names(patterns_df)) sum(patterns_df$significant, na.rm = TRUE) else NA
    ),
    parameters = list(
      schema = schema, min_length = min_length, max_length = max_length,
      start_state = start_state, start_schema = start_schema
    )
  )
  
  class(result) <- "meta_paths"
  if (verbose) cat("Meta-path discovery complete.\n")
  return(result)
}

#' Search for specific meta-path schema
#' @keywords internal
search_meta_path_schema <- function(sequences, type_sequences, schema,
                                    state_to_type, node_types, n_sequences,
                                    test_significance, correction, alpha) {
  
  parts <- strsplit(gsub(" ", "", schema), "->")[[1]]
  count <- 0; seqs_containing <- 0
  instance_list <- c(); instance_seqs <- list()
  
  for (seq_idx in seq_along(type_sequences)) {
    info <- type_sequences[[seq_idx]]
    matches <- find_pattern_matches(info$types, parts)
    
    if (length(matches) > 0) {
      # Need to extract original states for these matches. 
      # This logic is complex to generalize in find_pattern_matches, keeping specialized logic here for now
      # but simplified compared to before.
      
      # Actually, find_pattern_matches returns values, we need indices or to re-scan states.
      # For simplicity and correctness with multi-wildcards, we reuse the dedicated multi-matcher logic here inline
      # but reuse find_pattern_matches logic logic for non-wildcard or single-wildcard cases by mapping back.
      
      # Revert to specific meta-path matchers but keep them private/local or simplified
      if ("**" %in% parts) {
        # Multi-wildcard specific logic
        pos <- which(parts == "**")
        if (length(pos) == 1 && length(parts) == 3) {
          starts <- which(info$types == parts[1])
          ends <- which(info$types == parts[3])
          for (s in starts) for (e in ends) {
            if (e > s + 1) {
              inst <- paste(info$states[s:e], collapse = "->")
              if (!grepl("NA", inst)) {
                count <- count + 1
                instance_list <- c(instance_list, inst)
                if (is.null(instance_seqs[[inst]])) instance_seqs[[inst]] <- c()
                instance_seqs[[inst]] <- unique(c(instance_seqs[[inst]], seq_idx))
              }
            }
          }
        }
      } else {
        # Standard matching on TYPES, then map to STATES
        len <- length(parts)
        if (length(info$types) >= len) {
          for (start in 1:(length(info$types) - len + 1)) {
            matched <- TRUE
            for (i in 1:len) {
              if (parts[i] != "*" && parts[i] != info$types[start + i - 1]) { matched <- FALSE; break }
            }
            if (matched) {
              inst <- paste(info$states[start:(start + len - 1)], collapse = "->")
              if (!grepl("NA", inst)) {
                count <- count + 1
                instance_list <- c(instance_list, inst)
                if (is.null(instance_seqs[[inst]])) instance_seqs[[inst]] <- c()
                instance_seqs[[inst]] <- unique(c(instance_seqs[[inst]], seq_idx))
              }
            }
          }
        }
      }
    }
  }
  
  if (length(instance_seqs) > 0) seqs_containing <- length(unique(unlist(instance_seqs)))
  
  # Schema stats
  meta_paths_df <- data.frame(
    schema = schema, length = sum(!parts %in% c("*", "**")),
    count = count, sequences_containing = seqs_containing,
    support = round(seqs_containing / n_sequences, 4),
    stringsAsFactors = FALSE
  )
  
  if (test_significance && count > 0) {
    meta_paths_df <- compute_significance_stats(meta_paths_df, n_sequences, length(node_types), correction, alpha)
  }
  
  # Instances stats
  instances_df <- aggregate_instances(instance_list, instance_seqs, n_sequences, 1, # N_states not used for instance sig here
                                      schema = schema, test_significance, correction, alpha)
  
  return(list(meta_paths = meta_paths_df, instances = instances_df))
}

#' Discover meta-paths automatically
#' @keywords internal
discover_meta_paths <- function(sequences, type_sequences, type_names,
                                state_to_type, min_length, max_length,
                                n_sequences, min_support, min_count,
                                test_significance, correction, alpha, verbose) {
  
  meta_list <- list()
  inst_list <- list()
  
  for (len in min_length:max_length) {
    if (verbose) cat("  Extracting type", len, "-grams...\n")
    
    type_ngrams <- list()
    
    for (seq_idx in seq_along(type_sequences)) {
      info <- type_sequences[[seq_idx]]
      if (length(info$types) < len) next
      
      for (start in 1:(length(info$types) - len + 1)) {
        t_gram <- paste(info$types[start:(start + len - 1)], collapse = "->")
        s_gram <- paste(info$states[start:(start + len - 1)], collapse = "->")
        
        if (grepl("NA", s_gram)) next
        
        if (is.null(type_ngrams[[t_gram]])) type_ngrams[[t_gram]] <- list(count=0, seqs=c(), instances=c())
        type_ngrams[[t_gram]]$count <- type_ngrams[[t_gram]]$count + 1
        type_ngrams[[t_gram]]$seqs <- c(type_ngrams[[t_gram]]$seqs, seq_idx)
        type_ngrams[[t_gram]]$instances <- c(type_ngrams[[t_gram]]$instances, s_gram)
      }
    }
    
    for (schema in names(type_ngrams)) {
      info <- type_ngrams[[schema]]
      seqs_containing <- length(unique(info$seqs))
      
      if (info$count >= min_count && seqs_containing/n_sequences >= min_support) {
        meta_list[[length(meta_list)+1]] <- list(
          schema = schema, length = len, count = info$count, 
          sequences_containing = seqs_containing
        )
        
        # Process instances
        inst_df <- aggregate_instances(info$instances, info$seqs, n_sequences, 1,
                                       schema = schema, FALSE) # Defer stats
        inst_list[[length(inst_list)+1]] <- inst_df
      }
    }
  }
  
  if (length(meta_list) == 0) return(list(meta_paths = create_empty_patterns_df(TRUE), instances = create_empty_patterns_df(TRUE)))
  
  # Build meta_paths df
  meta_df <- do.call(rbind, lapply(meta_list, as.data.frame))
  meta_df$support <- round(meta_df$sequences_containing / n_sequences, 4)
  
  if (test_significance) {
    meta_df <- compute_significance_stats(meta_df, n_sequences, length(type_names), correction, alpha)
  }
  
  # Build instances df
  inst_df <- if (length(inst_list) > 0) do.call(rbind, inst_list) else create_empty_patterns_df(TRUE)
  if (test_significance && nrow(inst_df) > 0) {
    # Recalculate stats globally if needed, or rely on aggregate_instances if passed TRUE.
    # Here we computed basic stats, let's add significance now properly
    inst_df <- compute_significance_stats(inst_df, n_sequences, 1, correction, alpha)
  }
  
  return(list(
    meta_paths = meta_df[order(meta_df$count, decreasing = TRUE), ],
    instances = inst_df[order(inst_df$count, decreasing = TRUE), ]
  ))
}

#' Compute type-to-type transitions
#' @keywords internal
compute_type_transitions <- function(type_sequences, type_names, n_sequences,
                                     test_significance, correction, alpha) {
  
  n_types <- length(type_names)
  trans_counts <- matrix(0, nrow = n_types, ncol = n_types, dimnames = list(type_names, type_names))
  
  for (info in type_sequences) {
    t <- info$types
    if (length(t) < 2) next
    for (i in 1:(length(t) - 1)) {
      if (t[i] %in% type_names && t[i+1] %in% type_names) {
        trans_counts[t[i], t[i+1]] <- trans_counts[t[i], t[i+1]] + 1
      }
    }
  }
  
  trans_list <- list()
  for (from in type_names) {
    row_sum <- sum(trans_counts[from, ])
    for (to in type_names) {
      count <- trans_counts[from, to]
      if (count > 0) {
        trans_list[[length(trans_list)+1]] <- data.frame(
          from_type = from, to_type = to, count = count,
          probability = round(count / row_sum, 4), stringsAsFactors = FALSE
        )
      }
    }
  }
  
  if (length(trans_list) == 0) return(data.frame(from_type=character(0)))
  
  trans_df <- do.call(rbind, trans_list)
  
  if (test_significance) {
    trans_df$lift <- round(trans_df$probability / (1/n_types), 2)
    trans_df$significant <- trans_df$probability > (1/n_types) * 1.5
  }
  
  return(trans_df[order(trans_df$count, decreasing = TRUE), ])
}

# ==============================================================================
# S3 METHODS FOR META_PATHS
# ==============================================================================

#' @export
print.meta_paths <- function(x, ...) {
  cat("Meta-Path Analysis Results\n==========================\n\n")
  cat("Node Types:\n")
  for (n in names(x$type_mapping)) cat(sprintf("  %s: %s\n", n, paste(head(x$type_mapping[[n]], 5), collapse=", ")))
  cat("\nSummary:\n")
  cat("  Sequences:", x$summary$n_sequences, "\n")
  cat("  Patterns found:", x$summary$n_patterns, "\n")
  print_patterns_summary(x$patterns, "Top State-Level Patterns")
  invisible(x)
}

#' @export
summary.meta_paths <- function(object, ...) {
  print.meta_paths(object)
  if (!is.null(object$schemas) && nrow(object$schemas) > 0) {
    cat("TYPE-LEVEL SCHEMAS\n------------------\n")
    print(head(object$schemas, 10))
  }
  invisible(object)
}

# ==============================================================================
# BACKWARD COMPATIBILITY
# ==============================================================================

#' @rdname find_patterns
#' @export
find_gapped_patterns <- function(...) {
  .Deprecated("find_patterns")
  find_patterns(...)
}

cat("Sequence motif toolkit loaded.\n")

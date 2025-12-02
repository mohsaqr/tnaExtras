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
#'
#' # Filter by starting state
#' results <- find_patterns(seq_data, start_state = "plan")
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
  
  # Prepare data using shared utility
  seq_data <- prepare_sequence_data(data, verbose)
  sequences <- seq_data$sequences
  n_sequences <- seq_data$n_sequences
  n_states <- seq_data$n_states
  all_states <- seq_data$all_states
  
  # Convert pattern to standard format if provided as vector
  if (!is.null(pattern)) {
    if (is.character(pattern) && length(pattern) > 1) {
      pattern <- paste(pattern, collapse = "->")
    }
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
  
  patterns_df <- results$patterns
  instances_df <- results$instances
  
  # Apply state filters
  if (!is.null(start_state) || !is.null(end_state) || !is.null(contains_state)) {
    if (verbose) cat("Applying state filters...\n")
    patterns_df <- filter_by_state(patterns_df, start_state, end_state,
                                   contains_state, pattern_col = "pattern",
                                   separator = "->")
    # Also filter instances
    if (nrow(instances_df) > 0 && "instance" %in% names(instances_df)) {
      instances_df <- filter_by_state(instances_df, start_state, end_state,
                                      contains_state, pattern_col = "instance",
                                      separator = "->")
    }
    if (verbose) cat("  Patterns after filtering:", nrow(patterns_df), "\n")
  }
  
  # Standardize columns
  patterns_df <- standardize_columns(patterns_df, include_schema = FALSE)
  
  result <- list(
    patterns = patterns_df,
    instances = instances_df,
    summary = list(
      n_sequences = n_sequences,
      n_states = n_states,
      n_patterns = nrow(patterns_df),
      n_significant = if ("significant" %in% names(patterns_df))
        sum(patterns_df$significant, na.rm = TRUE) else NA,
      all_states = all_states
    ),
    parameters = list(
      pattern = pattern,
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
  
  if (verbose) cat("Pattern search complete.\n")
  
  return(result)
}

#' Search for specific gapped pattern
#' @keywords internal
search_gapped_pattern <- function(sequences, pattern, n_sequences, n_states,
                                  test_significance, correction, alpha) {
  
  parts <- strsplit(gsub(" ", "", pattern), "->")[[1]]
  
  count <- 0
  seqs_containing <- 0
  instance_list <- c()
  instance_seqs <- list()
  
  for (seq_idx in seq_along(sequences)) {
    seq <- sequences[[seq_idx]]
    
    if ("**" %in% parts) {
      matches <- find_multi_wildcard_matches(seq, parts)
    } else if ("*" %in% parts) {
      matches <- find_single_wildcard_matches(seq, parts)
    } else {
      matches <- find_exact_matches(seq, parts)
    }
    
    if (length(matches) > 0) {
      for (m in matches) {
        inst_str <- paste(m, collapse = "->")
        if (!grepl("NA", inst_str)) {
          count <- count + 1
          instance_list <- c(instance_list, inst_str)
          if (is.null(instance_seqs[[inst_str]])) {
            instance_seqs[[inst_str]] <- c()
          }
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
    n_fixed <- sum(!parts %in% c("*", "**"))
    expected_prob <- (1 / n_states) ^ n_fixed
    patterns_df$expected_prob <- round(expected_prob, 6)
    patterns_df$lift <- round(support / expected_prob, 2)
    patterns_df$p_value <- stats::binom.test(seqs_containing, n_sequences,
                                             expected_prob, alternative = "greater")$p.value
    patterns_df$p_adjusted <- patterns_df$p_value
    patterns_df$significant <- patterns_df$p_adjusted < alpha
  }
  
  # Build instances data frame with statistics
  if (length(instance_list) > 0) {
    inst_table <- sort(table(instance_list), decreasing = TRUE)
    instances_df <- data.frame(
      pattern = rep(pattern, length(inst_table)),
      instance = names(inst_table),
      count = as.integer(inst_table),
      sequences_containing = sapply(names(inst_table), function(p) length(instance_seqs[[p]])),
      stringsAsFactors = FALSE
    )
    instances_df$support <- round(instances_df$sequences_containing / n_sequences, 4)
    rownames(instances_df) <- NULL
  } else {
    instances_df <- data.frame(
      pattern = character(0),
      instance = character(0),
      count = integer(0),
      sequences_containing = integer(0),
      support = numeric(0),
      stringsAsFactors = FALSE
    )
  }
  
  return(list(patterns = patterns_df, instances = instances_df))
}

#' Find exact pattern matches
#' @keywords internal
find_exact_matches <- function(seq, parts) {
  matches <- list()
  pattern_len <- length(parts)
  
  if (length(seq) < pattern_len) return(matches)
  
  for (start in 1:(length(seq) - pattern_len + 1)) {
    if (all(seq[start:(start + pattern_len - 1)] == parts)) {
      matches[[length(matches) + 1]] <- seq[start:(start + pattern_len - 1)]
    }
  }
  
  return(matches)
}

#' Find single wildcard matches
#' @keywords internal
find_single_wildcard_matches <- function(seq, parts) {
  matches <- list()
  pattern_len <- length(parts)
  
  if (length(seq) < pattern_len) return(matches)
  if (length(parts) == 0) return(matches)
  
  for (start in 1:(length(seq) - pattern_len + 1)) {
    matched <- TRUE
    match_seq <- character(pattern_len)
    
    for (i in 1:pattern_len) {
      if (parts[i] == "*") {
        match_seq[i] <- seq[start + i - 1]
      } else if (parts[i] == seq[start + i - 1]) {
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

#' Find multi-wildcard matches
#' @keywords internal
find_multi_wildcard_matches <- function(seq, parts) {
  matches <- list()
  
  # Simple case: A->**->B
  multi_pos <- which(parts == "**")
  
  if (length(multi_pos) == 1 && length(parts) == 3) {
    first_state <- parts[1]
    last_state <- parts[3]
    
    first_positions <- which(seq == first_state)
    last_positions <- which(seq == last_state)
    
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

#' Discover gapped patterns automatically
#' @keywords internal
discover_gapped_patterns <- function(sequences, all_states, min_gap, max_gap,
                                     n_sequences, n_states, min_support, min_count,
                                     test_significance, correction, alpha, verbose) {
  
  patterns_list <- list()
  instances_list <- list()
  
  for (state1 in all_states) {
    for (state2 in all_states) {
      for (gap in min_gap:max_gap) {
        pattern_name <- paste0(state1, "->", paste(rep("*", gap), collapse = "->"), "->", state2)
        
        count <- 0
        seqs_containing <- 0
        pattern_instances <- c()
        inst_seqs <- list()
        
        for (seq_idx in seq_along(sequences)) {
          seq <- sequences[[seq_idx]]
          if (length(seq) < gap + 2) next
          
          seq_count <- 0
          for (i in 1:(length(seq) - gap - 1)) {
            if (seq[i] == state1 && seq[i + gap + 1] == state2) {
              seq_count <- seq_count + 1
              inst <- paste(seq[i:(i + gap + 1)], collapse = "->")
              pattern_instances <- c(pattern_instances, inst)
              if (is.null(inst_seqs[[inst]])) inst_seqs[[inst]] <- c()
              inst_seqs[[inst]] <- unique(c(inst_seqs[[inst]], seq_idx))
            }
          }
          
          if (seq_count > 0) {
            count <- count + seq_count
            seqs_containing <- seqs_containing + 1
          }
        }
        
        support <- seqs_containing / n_sequences
        
        if (count >= min_count && support >= min_support) {
          patterns_list[[length(patterns_list) + 1]] <- list(
            pattern = pattern_name,
            from = state1,
            to = state2,
            gap = gap,
            count = count,
            sequences_containing = seqs_containing
          )
          
          # Store instances
          if (length(pattern_instances) > 0) {
            inst_table <- table(pattern_instances)
            for (inst in names(inst_table)) {
              instances_list[[length(instances_list) + 1]] <- list(
                pattern = pattern_name,
                instance = inst,
                count = as.integer(inst_table[inst]),
                sequences_containing = length(inst_seqs[[inst]])
              )
            }
          }
        }
      }
    }
  }
  
  if (length(patterns_list) == 0) {
    return(list(
      patterns = create_empty_patterns_df(),
      instances = data.frame(
        pattern = character(0), instance = character(0),
        count = integer(0), sequences_containing = integer(0),
        support = numeric(0), stringsAsFactors = FALSE
      )
    ))
  }
  
  # Build patterns data frame
  patterns_df <- data.frame(
    pattern = sapply(patterns_list, `[[`, "pattern"),
    length = sapply(patterns_list, function(x) x$gap + 2),
    count = sapply(patterns_list, `[[`, "count"),
    sequences_containing = sapply(patterns_list, `[[`, "sequences_containing"),
    stringsAsFactors = FALSE
  )
  
  patterns_df$support <- round(patterns_df$sequences_containing / n_sequences, 4)
  patterns_df$proportion <- round(patterns_df$count / sum(patterns_df$count), 4)
  
  if (test_significance && nrow(patterns_df) > 0) {
    patterns_df <- compute_significance_stats(patterns_df, n_sequences, n_states,
                                              correction, alpha)
  }
  
  patterns_df <- patterns_df[order(patterns_df$count, decreasing = TRUE), ]
  rownames(patterns_df) <- NULL
  
  # Build instances data frame
  if (length(instances_list) > 0) {
    instances_df <- data.frame(
      pattern = sapply(instances_list, `[[`, "pattern"),
      instance = sapply(instances_list, `[[`, "instance"),
      count = sapply(instances_list, `[[`, "count"),
      sequences_containing = sapply(instances_list, `[[`, "sequences_containing"),
      stringsAsFactors = FALSE
    )
    instances_df$support <- round(instances_df$sequences_containing / n_sequences, 4)
    instances_df <- instances_df[order(instances_df$count, decreasing = TRUE), ]
    rownames(instances_df) <- NULL
  } else {
    instances_df <- data.frame(
      pattern = character(0), instance = character(0),
      count = integer(0), sequences_containing = integer(0),
      support = numeric(0), stringsAsFactors = FALSE
    )
  }
  
  if (verbose) cat("  Gapped patterns found:", nrow(patterns_df), "\n")
  
  return(list(patterns = patterns_df, instances = instances_df))
}

# ==============================================================================
# META-PATH DISCOVERY
# ==============================================================================

#' Find Meta-Paths Across Node Types
#'
#' Discovers meta-paths (type-level patterns) in sequential data where states
#' are grouped into node types. This enables analysis of heterogeneous networks
#' embedded in sequences, revealing patterns like "cognitive->social->cognitive".
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
#' @param min_support Minimum support threshold (default: 0.01)
#' @param min_count Minimum count threshold (default: 2)
#' @param test_significance Whether to perform significance tests (default: TRUE)
#' @param correction Multiple testing correction method (default: "fdr")
#' @param alpha Significance level (default: 0.05)
#' @param start_state Filter patterns starting with this state (default: NULL)
#' @param end_state Filter patterns ending with this state (default: NULL)
#' @param contains_state Filter patterns containing this state (default: NULL)
#' @param start_schema Filter by starting schema/type (default: NULL)
#' @param end_schema Filter by ending schema/type (default: NULL)
#' @param contains_schema Filter by schema containing this type (default: NULL)
#' @param verbose Whether to print progress messages (default: TRUE)
#'
#' @return An object of class "meta_paths" containing:
#'   \item{patterns}{Data frame of state-level patterns with statistics}
#'   \item{schemas}{Data frame of type-level patterns with statistics}
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
#' # Define node types
#' node_types <- list(
#'   cognitive = c("plan", "monitor", "adapt"),
#'   social = c("discuss", "consensus", "coregulate", "synthesis"),
#'   emotional = c("emotion", "cohesion")
#' )
#'
#' # Auto-discover all meta-paths
#' meta <- find_meta_paths(seq_data, node_types = node_types)
#' print(meta)
#'
#' # Search for specific schema
#' meta <- find_meta_paths(seq_data, node_types = node_types,
#'                         schema = c("cognitive", "social", "cognitive"))
#'
#' # Filter by starting state
#' meta <- find_meta_paths(seq_data, node_types = node_types,
#'                         start_state = "plan")
#'
#' # Filter by schema type
#' meta <- find_meta_paths(seq_data, node_types = node_types,
#'                         start_schema = "cognitive")
#' }
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
  
  # Validate node_types
  if (!is.list(node_types) || length(node_types) < 2) {
    stop("node_types must be a named list with at least 2 type definitions")
  }
  
  type_names <- names(node_types)
  if (is.null(type_names) || any(type_names == "")) {
    stop("node_types must have named elements")
  }
  
  # Prepare data using shared utility
  seq_data <- prepare_sequence_data(data, verbose)
  sequences <- seq_data$sequences
  n_sequences <- seq_data$n_sequences
  all_states <- seq_data$all_states
  
  # Convert schema to standard format
  if (!is.null(schema)) {
    if (is.character(schema) && length(schema) > 1) {
      schema <- paste(schema, collapse = "->")
    }
  }
  
  # Create state-to-type mapping
  state_to_type <- create_state_type_mapping(node_types, all_states)
  
  # Check coverage
  mapped_states <- names(state_to_type)
  unmapped <- setdiff(all_states, mapped_states)
  coverage <- length(mapped_states) / length(all_states)
  
  if (verbose) {
    cat("  States mapped:", length(mapped_states), "/", length(all_states),
        sprintf(" (%.1f%%)\n", coverage * 100))
    if (length(unmapped) > 0) {
      cat("  Unmapped:", paste(head(unmapped, 5), collapse = ", "),
          if (length(unmapped) > 5) "..." else "", "\n")
    }
  }
  
  if (length(mapped_states) == 0) {
    stop("No states in data match the node_types definitions.\n",
         "  States in data: ", paste(head(all_states, 10), collapse = ", "), "\n",
         "  States in node_types: ", paste(head(unlist(node_types), 10), collapse = ", "))
  }
  
  # Convert sequences to type sequences
  type_sequences <- lapply(sequences, function(seq) {
    types <- sapply(seq, function(s) {
      if (s %in% names(state_to_type)) state_to_type[[s]] else NA
    })
    valid_idx <- !is.na(types)
    list(types = types[valid_idx], states = seq[valid_idx])
  })
  
  # Filter sequences with minimum length
  type_sequences <- type_sequences[sapply(type_sequences, function(x) length(x$types)) >= min_length]
  
  if (length(type_sequences) == 0) {
    stop("No valid type sequences after mapping")
  }
  
  n_types <- length(type_names)
  
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
  
  patterns_df <- results$instances
  schemas_df <- results$meta_paths
  
  # Apply state-level filters
  if (!is.null(start_state) || !is.null(end_state) || !is.null(contains_state)) {
    if (verbose) cat("Applying state filters...\n")
    patterns_df <- filter_by_state(patterns_df, start_state, end_state,
                                   contains_state, pattern_col = "pattern",
                                   separator = "->")
  }
  
  # Apply schema-level filters
  if (!is.null(start_schema) || !is.null(end_schema) || !is.null(contains_schema)) {
    if (verbose) cat("Applying schema filters...\n")
    patterns_df <- filter_by_schema(patterns_df, start_schema, end_schema, contains_schema)
    schemas_df <- filter_by_schema(schemas_df, start_schema, end_schema, contains_schema)
  }
  
  # Compute type transitions
  type_transitions <- compute_type_transitions(type_sequences, type_names,
                                               n_sequences, test_significance,
                                               correction, alpha)
  
  result <- list(
    patterns = patterns_df,
    schemas = schemas_df,
    type_transitions = type_transitions,
    type_mapping = node_types,
    summary = list(
      n_sequences = n_sequences,
      n_states = length(all_states),
      n_types = n_types,
      type_names = type_names,
      coverage = coverage,
      n_patterns = nrow(patterns_df),
      n_schemas = nrow(schemas_df),
      n_significant = if ("significant" %in% names(patterns_df))
        sum(patterns_df$significant, na.rm = TRUE) else NA
    ),
    parameters = list(
      schema = schema,
      min_length = min_length,
      max_length = max_length,
      min_support = min_support,
      min_count = min_count,
      correction = correction,
      alpha = alpha,
      start_state = start_state,
      end_state = end_state,
      contains_state = contains_state,
      start_schema = start_schema,
      end_schema = end_schema,
      contains_schema = contains_schema
    )
  )
  
  class(result) <- "meta_paths"
  
  if (verbose) cat("Meta-path discovery complete.\n")
  
  return(result)
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
  n_fixed <- sum(!schema_parts %in% c("*", "**"))
  
  count <- 0
  seqs_containing <- 0
  instance_list <- c()
  instance_seqs <- list()
  
  has_multi <- "**" %in% schema_parts
  
  for (seq_idx in seq_along(type_sequences)) {
    type_info <- type_sequences[[seq_idx]]
    type_seq <- type_info$types
    state_seq <- type_info$states
    
    if (length(type_seq) < n_fixed) next
    
    if (has_multi) {
      matches <- find_meta_path_multi_matches(state_seq, type_seq, schema_parts, state_to_type)
    } else {
      matches <- find_meta_path_matches(state_seq, type_seq, schema_parts)
    }
    
    if (length(matches) > 0) {
      seq_had_valid <- FALSE
      for (m in matches) {
        inst_str <- paste(m, collapse = "->")
        if (!grepl("NA", inst_str)) {
          count <- count + 1
          instance_list <- c(instance_list, inst_str)
          seq_had_valid <- TRUE
          if (is.null(instance_seqs[[inst_str]])) instance_seqs[[inst_str]] <- c()
          instance_seqs[[inst_str]] <- unique(c(instance_seqs[[inst_str]], seq_idx))
        }
      }
      if (seq_had_valid) seqs_containing <- seqs_containing + 1
    }
  }
  
  support <- seqs_containing / n_sequences
  
  # Schema-level statistics
  meta_paths_df <- data.frame(
    schema = schema,
    length = n_fixed,
    count = count,
    sequences_containing = seqs_containing,
    support = round(support, 4),
    stringsAsFactors = FALSE
  )
  
  if (test_significance && count > 0) {
    expected_prob <- (1 / n_types) ^ n_fixed
    meta_paths_df$expected_prob <- round(expected_prob, 6)
    meta_paths_df$lift <- round(support / expected_prob, 2)
    meta_paths_df$p_value <- stats::binom.test(seqs_containing, n_sequences,
                                               expected_prob, alternative = "greater")$p.value
    meta_paths_df$p_adjusted <- meta_paths_df$p_value
    meta_paths_df$significant <- meta_paths_df$p_adjusted < alpha
  }
  
  # Instance-level statistics
  if (length(instance_list) > 0) {
    inst_table <- sort(table(instance_list), decreasing = TRUE)
    instances_df <- data.frame(
      pattern = names(inst_table),
      schema = rep(schema, length(inst_table)),
      length = sapply(strsplit(names(inst_table), "->"), length),
      count = as.integer(inst_table),
      sequences_containing = sapply(names(inst_table), function(p) length(instance_seqs[[p]])),
      stringsAsFactors = FALSE
    )
    instances_df$support <- round(instances_df$sequences_containing / n_sequences, 4)
    instances_df$proportion <- round(instances_df$count / sum(instances_df$count), 4)
    
    if (test_significance) {
      n_unique <- length(unique(instances_df$pattern))
      instances_df$expected_prob <- round(1 / max(n_unique, 1), 6)
      instances_df$lift <- round(instances_df$support / instances_df$expected_prob, 2)
      instances_df$p_value <- sapply(1:nrow(instances_df), function(i) {
        stats::binom.test(instances_df$sequences_containing[i], n_sequences,
                          instances_df$expected_prob[i], alternative = "greater")$p.value
      })
      instances_df$p_adjusted <- stats::p.adjust(instances_df$p_value, method = correction)
      instances_df$significant <- instances_df$p_adjusted < alpha
    }
    rownames(instances_df) <- NULL
  } else {
    instances_df <- create_empty_patterns_df(include_schema = TRUE)
  }
  
  return(list(meta_paths = meta_paths_df, instances = instances_df))
}

#' Find meta-path matches
#' @keywords internal
find_meta_path_matches <- function(state_seq, type_seq, schema_parts) {
  matches <- list()
  pattern_len <- length(schema_parts)
  
  if (length(type_seq) < pattern_len) return(matches)
  
  for (start in 1:(length(type_seq) - pattern_len + 1)) {
    matched <- TRUE
    for (i in 1:pattern_len) {
      if (schema_parts[i] != "*" && schema_parts[i] != type_seq[start + i - 1]) {
        matched <- FALSE
        break
      }
    }
    if (matched) {
      matches[[length(matches) + 1]] <- state_seq[start:(start + pattern_len - 1)]
    }
  }
  
  return(matches)
}

#' Find meta-path matches with multi-wildcards
#' @keywords internal
find_meta_path_multi_matches <- function(state_seq, type_seq, schema_parts, state_to_type) {
  matches <- list()
  
  multi_pos <- which(schema_parts == "**")
  
  if (length(multi_pos) == 1 && length(schema_parts) == 3) {
    first_type <- schema_parts[1]
    last_type <- schema_parts[3]
    
    first_positions <- which(type_seq == first_type)
    last_positions <- which(type_seq == last_type)
    
    for (fp in first_positions) {
      for (lp in last_positions) {
        if (lp > fp + 1) {
          matches[[length(matches) + 1]] <- state_seq[fp:lp]
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
  
  for (len in min_length:max_length) {
    if (verbose) cat("  Extracting type", len, "-grams...\n")
    
    type_ngrams <- list()
    
    for (seq_idx in seq_along(type_sequences)) {
      type_info <- type_sequences[[seq_idx]]
      type_seq <- type_info$types
      state_seq <- type_info$states
      
      if (length(type_seq) < len) next
      
      for (start in 1:(length(type_seq) - len + 1)) {
        type_gram <- paste(type_seq[start:(start + len - 1)], collapse = "->")
        state_instance <- paste(state_seq[start:(start + len - 1)], collapse = "->")
        
        if (grepl("NA", state_instance)) next
        
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
          support = round(support, 4),
          stringsAsFactors = FALSE
        )
        
        # Store instances with sequence tracking
        inst_table <- sort(table(info$instances), decreasing = TRUE)
        for (inst in names(inst_table)) {
          inst_seqs <- unique(info$seqs[info$instances == inst])
          instances_list[[length(instances_list) + 1]] <- list(
            pattern = inst,
            schema = schema,
            length = len,
            count = as.integer(inst_table[inst]),
            sequences_containing = length(inst_seqs)
          )
        }
      }
    }
  }
  
  if (length(meta_paths_list) == 0) {
    return(list(
      meta_paths = data.frame(schema = character(0), length = integer(0),
                              count = integer(0), sequences_containing = integer(0),
                              support = numeric(0), stringsAsFactors = FALSE),
      instances = create_empty_patterns_df(include_schema = TRUE)
    ))
  }
  
  meta_paths_df <- do.call(rbind, meta_paths_list)
  rownames(meta_paths_df) <- NULL
  
  # Add significance to schemas
  if (test_significance && nrow(meta_paths_df) > 0) {
    meta_paths_df$expected_prob <- round((1 / n_types) ^ meta_paths_df$length, 6)
    meta_paths_df$lift <- round(meta_paths_df$support / meta_paths_df$expected_prob, 2)
    meta_paths_df$p_value <- sapply(1:nrow(meta_paths_df), function(i) {
      stats::binom.test(meta_paths_df$sequences_containing[i], n_sequences,
                        meta_paths_df$expected_prob[i], alternative = "greater")$p.value
    })
    meta_paths_df$p_adjusted <- stats::p.adjust(meta_paths_df$p_value, method = correction)
    meta_paths_df$significant <- meta_paths_df$p_adjusted < alpha
  }
  
  meta_paths_df <- meta_paths_df[order(meta_paths_df$count, decreasing = TRUE), ]
  
  # Build instances data frame
  if (length(instances_list) > 0) {
    instances_df <- data.frame(
      pattern = sapply(instances_list, `[[`, "pattern"),
      schema = sapply(instances_list, `[[`, "schema"),
      length = sapply(instances_list, `[[`, "length"),
      count = sapply(instances_list, `[[`, "count"),
      sequences_containing = sapply(instances_list, `[[`, "sequences_containing"),
      stringsAsFactors = FALSE
    )
    instances_df$support <- round(instances_df$sequences_containing / n_sequences, 4)
    instances_df$proportion <- round(instances_df$count / sum(instances_df$count), 4)
    
    if (test_significance) {
      n_unique <- length(unique(instances_df$pattern))
      instances_df$expected_prob <- round(1 / max(n_unique, 1), 6)
      instances_df$lift <- round(instances_df$support / instances_df$expected_prob, 2)
      instances_df$p_value <- sapply(1:nrow(instances_df), function(i) {
        stats::binom.test(instances_df$sequences_containing[i], n_sequences,
                          instances_df$expected_prob[i], alternative = "greater")$p.value
      })
      instances_df$p_adjusted <- stats::p.adjust(instances_df$p_value, method = correction)
      instances_df$significant <- instances_df$p_adjusted < alpha
    }
    
    instances_df <- instances_df[order(instances_df$count, decreasing = TRUE), ]
    rownames(instances_df) <- NULL
  } else {
    instances_df <- create_empty_patterns_df(include_schema = TRUE)
  }
  
  if (verbose) cat("  Meta-paths found:", nrow(meta_paths_df), "\n")
  
  return(list(meta_paths = meta_paths_df, instances = instances_df))
}

#' Compute type-to-type transitions
#' @keywords internal
compute_type_transitions <- function(type_sequences, type_names, n_sequences,
                                     test_significance, correction, alpha) {
  
  n_types <- length(type_names)
  trans_counts <- matrix(0, nrow = n_types, ncol = n_types,
                         dimnames = list(type_names, type_names))
  
  for (type_info in type_sequences) {
    type_seq <- if (is.list(type_info)) type_info$types else type_info
    if (length(type_seq) < 2) next
    for (i in 1:(length(type_seq) - 1)) {
      from <- type_seq[i]
      to <- type_seq[i + 1]
      if (from %in% type_names && to %in% type_names) {
        trans_counts[from, to] <- trans_counts[from, to] + 1
      }
    }
  }
  
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
          probability = round(prob, 4),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  
  if (length(trans_list) == 0) {
    return(data.frame(from_type = character(0), to_type = character(0),
                      count = integer(0), probability = numeric(0),
                      stringsAsFactors = FALSE))
  }
  
  trans_df <- do.call(rbind, trans_list)
  rownames(trans_df) <- NULL
  
  if (test_significance) {
    expected_prob <- 1 / n_types
    trans_df$lift <- round(trans_df$probability / expected_prob, 2)
    trans_df$significant <- trans_df$probability > expected_prob * 1.5
  }
  
  trans_df <- trans_df[order(trans_df$count, decreasing = TRUE), ]
  
  return(trans_df)
}

# ==============================================================================
# S3 METHODS FOR META_PATHS
# ==============================================================================

#' @export
print.meta_paths <- function(x, ...) {
  cat("Meta-Path Analysis Results\n")
  cat("==========================\n\n")
  
  cat("Node Types:\n")
  for (type_name in names(x$type_mapping)) {
    states <- x$type_mapping[[type_name]]
    cat(sprintf("  %s: %s (%d states)\n", type_name,
                paste(head(states, 5), collapse = ", "),
                length(states)))
  }
  cat("\n")
  
  cat("Summary:\n")
  cat("  Sequences:", x$summary$n_sequences, "\n")
  cat("  Coverage:", sprintf("%.1f%%", x$summary$coverage * 100), "\n")
  cat("  Patterns found:", x$summary$n_patterns, "\n")
  if (!is.na(x$summary$n_significant)) {
    cat("  Significant:", x$summary$n_significant, "\n")
  }
  cat("\n")
  
  # Print state-level patterns (primary)
  print_patterns_summary(x$patterns, "Top State-Level Patterns", n_show = 15)
  
  # Print schema summary (secondary)
  if (!is.null(x$schemas) && nrow(x$schemas) > 0) {
    cat("Type-Level Schemas:\n")
    for (i in 1:min(5, nrow(x$schemas))) {
      cat(sprintf("  %s: count=%d, support=%.3f\n",
                  x$schemas$schema[i], x$schemas$count[i], x$schemas$support[i]))
    }
    cat("\n")
  }
  
  invisible(x)
}

#' @export
summary.meta_paths <- function(object, ...) {
  cat("Meta-Path Analysis - Detailed Summary\n")
  cat("=====================================\n\n")
  
  cat("NODE TYPES\n")
  cat(paste(rep("-", 50), collapse = ""), "\n")
  for (type_name in names(object$type_mapping)) {
    states <- object$type_mapping[[type_name]]
    cat(sprintf("%s: %s\n", type_name, paste(states, collapse = ", ")))
  }
  cat("\n")
  
  cat("STATE-LEVEL PATTERNS\n")
  cat(paste(rep("-", 50), collapse = ""), "\n")
  if (!is.null(object$patterns) && nrow(object$patterns) > 0) {
    display_cols <- c("pattern", "schema", "count", "sequences_containing",
                      "support", "lift", "p_adjusted", "significant")
    display_cols <- display_cols[display_cols %in% names(object$patterns)]
    print(head(object$patterns[, display_cols, drop = FALSE], 25), row.names = FALSE)
  }
  cat("\n")
  
  cat("TYPE-LEVEL SCHEMAS\n")
  cat(paste(rep("-", 50), collapse = ""), "\n")
  if (!is.null(object$schemas) && nrow(object$schemas) > 0) {
    print(object$schemas, row.names = FALSE)
  }
  cat("\n")
  
  cat("TYPE TRANSITIONS\n")
  cat(paste(rep("-", 50), collapse = ""), "\n")
  if (!is.null(object$type_transitions) && nrow(object$type_transitions) > 0) {
    print(head(object$type_transitions, 15), row.names = FALSE)
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

#' @rdname explore_patterns
#' @export
detect_abstract_patterns <- function(data, patterns = "all", ...) {
  .Deprecated("explore_patterns")
  explore_patterns(data, type = "abstract", ...)
}

cat("Sequence motif toolkit loaded.\n")
cat("Functions: find_patterns(), find_meta_paths()\n")

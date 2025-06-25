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
# INTERNAL HELPER FUNCTIONS (some might be candidates for utils.R if broadly applicable)
# ==============================================================================

#' @noRd
.extract_ngrams_from_sequence <- function(sequence_string, n_gram_length) {
  if (!nzchar(sequence_string)) return(character(0)) # Use nzchar for empty string check
  
  seq_parts <- unlist(strsplit(sequence_string, "-", fixed = TRUE)) # fixed = TRUE for performance
  # Consider if filtering single characters is always desired or should be an option
  # seq_parts <- seq_parts[nchar(seq_parts) > 1]
  
  if (length(seq_parts) < n_gram_length) return(character(0))
  
  # More efficient n-gram generation using sapply or rolling apply if performance becomes an issue
  # For now, loop is fine.
  ngrams <- character(length(seq_parts) - n_gram_length + 1)
  for (i in 1:(length(seq_parts) - n_gram_length + 1)) {
    ngrams[i] <- paste(seq_parts[i:(i + n_gram_length - 1)], collapse = "-")
  }
  return(ngrams)
}


# ==============================================================================
# PAIRWISE PATTERN METRIC CALCULATION FUNCTIONS
# These are the core computation engines for a pair of groups.
# They now accept group names to dynamically name output columns.
# ==============================================================================

#' @noRd
.compute_pairwise_support <- function(patterns_to_analyze, sequences_group1, sequences_group2, name_group1="G1", name_group2="G2") {
  if (!is.character(patterns_to_analyze) || length(patterns_to_analyze) == 0) {
    stop("patterns_to_analyze must be a non-empty character vector", call. = FALSE)
  }
  # Further input validation for sequences_group1, sequences_group2 omitted for brevity, assume valid
  
  n_seq_g1 <- length(sequences_group1)
  n_seq_g2 <- length(sequences_group2)

  if (n_seq_g1 == 0 || n_seq_g2 == 0) {
      stop("Both groups must contain sequences for support calculation.", call. = FALSE)
  }

  results_df <- data.frame(
    pattern = patterns_to_analyze,
    stringsAsFactors = FALSE
  )
  # Dynamically name columns
  col_support_g1 <- paste0("support_", name_group1)
  col_support_g2 <- paste0("support_", name_group2)

  results_df[[col_support_g1]] <- numeric(length(patterns_to_analyze))
  results_df[[col_support_g2]] <- numeric(length(patterns_to_analyze))
  results_df$support_diff <- numeric(length(patterns_to_analyze))
  results_df$support_ratio <- numeric(length(patterns_to_analyze))
  results_df$relative_support <- numeric(length(patterns_to_analyze))

  for (i in seq_along(patterns_to_analyze)) {
    current_pattern <- patterns_to_analyze[i]
    count_g1 <- sum(sapply(sequences_group1, function(s) grepl(current_pattern, s, fixed = TRUE)))
    count_g2 <- sum(sapply(sequences_group2, function(s) grepl(current_pattern, s, fixed = TRUE)))
    
    support_g1 <- count_g1 / n_seq_g1
    support_g2 <- count_g2 / n_seq_g2
    
    results_df[[col_support_g1]][i] <- support_g1
    results_df[[col_support_g2]][i] <- support_g2
    results_df$support_diff[i] <- support_g1 - support_g2
    results_df$support_ratio[i] <- (support_g1 + 1e-9) / (support_g2 + 1e-9) # Avoid division by zero
    total_support <- support_g1 + support_g2
    results_df$relative_support[i] <- ifelse(total_support == 0, 0, (support_g1 - support_g2) / total_support)
  }
  
  results_df <- results_df[order(abs(results_df$support_diff), decreasing = TRUE), ]
  rownames(results_df) <- NULL
  return(results_df)
}

#' @noRd
.compute_pairwise_lift <- function(patterns_to_analyze, sequences_group1, sequences_group2, name_group1="G1", name_group2="G2") {
    n_seq_g1 <- length(sequences_group1)
    n_seq_g2 <- length(sequences_group2)
    n_total_seq <- n_seq_g1 + n_seq_g2

    if (n_seq_g1 == 0 || n_seq_g2 == 0) stop("Both groups must contain sequences for lift calculation.", call. = FALSE)

    results_df <- data.frame(pattern = patterns_to_analyze, stringsAsFactors = FALSE)
    col_obs_g1 <- paste0("observed_", name_group1); col_obs_g2 <- paste0("observed_", name_group2)
    col_exp_g1 <- paste0("expected_", name_group1); col_exp_g2 <- paste0("expected_", name_group2)
    col_lift_g1 <- paste0("lift_", name_group1); col_lift_g2 <- paste0("lift_", name_group2)

    # Initialize columns
    for(col_name in c(col_obs_g1, col_obs_g2, col_exp_g1, col_exp_g2, col_lift_g1, col_lift_g2, "lift_ratio", "lift_difference")) {
        results_df[[col_name]] <- numeric(length(patterns_to_analyze))
    }

    for (i in seq_along(patterns_to_analyze)) {
        current_pattern <- patterns_to_analyze[i]
        count_g1 <- sum(sapply(sequences_group1, function(s) grepl(current_pattern, s, fixed = TRUE)))
        count_g2 <- sum(sapply(sequences_group2, function(s) grepl(current_pattern, s, fixed = TRUE)))
        count_total_pattern <- count_g1 + count_g2

        results_df[[col_obs_g1]][i] <- count_g1
        results_df[[col_obs_g2]][i] <- count_g2

        if (count_total_pattern == 0) next # All other metrics will be 0 or NA

        expected_g1 <- count_total_pattern * (n_seq_g1 / n_total_seq)
        expected_g2 <- count_total_pattern * (n_seq_g2 / n_total_seq)
        results_df[[col_exp_g1]][i] <- expected_g1
        results_df[[col_exp_g2]][i] <- expected_g2

        prop_g1 <- count_g1 / n_seq_g1
        prop_g2 <- count_g2 / n_seq_g2

        # Proportion of sequences containing the pattern in the total dataset
        prop_pattern_overall <- count_total_pattern / n_total_seq

        # Lift for Group 1: P(Pattern | Group1) / P(Pattern) = prop_g1 / prop_pattern_overall
        # Lift for Group 2: P(Pattern | Group2) / P(Pattern) = prop_g2 / prop_pattern_overall
        # Avoid division by zero if prop_pattern_overall is 0 (though count_total_pattern > 0 check handles this)

        lift_g1 <- if(prop_pattern_overall > 1e-9) prop_g1 / prop_pattern_overall else 0
        lift_g2 <- if(prop_pattern_overall > 1e-9) prop_g2 / prop_pattern_overall else 0

        results_df[[col_lift_g1]][i] <- lift_g1
        results_df[[col_lift_g2]][i] <- lift_g2
        results_df$lift_ratio[i] <- pmin(pmax((lift_g1 + 1e-9) / (lift_g2 + 1e-9), 0.01), 100)
        results_df$lift_difference[i] <- lift_g1 - lift_g2
    }
    results_df <- results_df[order(abs(results_df$lift_difference), decreasing = TRUE), ]
    rownames(results_df) <- NULL
    return(results_df)
}

#' @noRd
.compute_pairwise_confidence <- function(patterns_to_analyze, sequences_group1, sequences_group2, name_group1="G1", name_group2="G2") {
    n_seq_g1 <- length(sequences_group1)
    n_seq_g2 <- length(sequences_group2)
    if (n_seq_g1 == 0 || n_seq_g2 == 0) stop("Both groups must contain sequences for confidence calculation.", call. = FALSE)

    results_df <- data.frame(pattern = patterns_to_analyze, stringsAsFactors = FALSE)
    col_conf_g1 <- paste0("confidence_", name_group1); col_conf_g2 <- paste0("confidence_", name_group2)
    
    for(col_name in c(col_conf_g1, col_conf_g2, "confidence_diff", "confidence_ratio", "balanced_confidence")) {
        results_df[[col_name]] <- numeric(length(patterns_to_analyze))
    }

    for (i in seq_along(patterns_to_analyze)) {
        current_pattern <- patterns_to_analyze[i]
        count_g1 <- sum(sapply(sequences_group1, function(s) grepl(current_pattern, s, fixed = TRUE)))
        count_g2 <- sum(sapply(sequences_group2, function(s) grepl(current_pattern, s, fixed = TRUE)))
        count_total_pattern <- count_g1 + count_g2

        if (count_total_pattern == 0) next

        # Confidence = P(Group | Pattern) = Count(Pattern in Group) / Count(Pattern Total)
        conf_g1 <- count_g1 / count_total_pattern
        conf_g2 <- count_g2 / count_total_pattern

        results_df[[col_conf_g1]][i] <- conf_g1
        results_df[[col_conf_g2]][i] <- conf_g2
        results_df$confidence_diff[i] <- conf_g1 - conf_g2
        results_df$confidence_ratio[i] <- (conf_g1 + 1e-9) / (conf_g2 + 1e-9)

        # Balanced confidence (original formula: (conf_A - expected_A_prop) - (conf_B - expected_B_prop))
        # where expected_A_prop is n_A / (n_A + n_B)
        prop_overall_g1 <- n_seq_g1 / (n_seq_g1 + n_seq_g2)
        prop_overall_g2 <- n_seq_g2 / (n_seq_g1 + n_seq_g2)
        results_df$balanced_confidence[i] <- (conf_g1 - prop_overall_g1) - (conf_g2 - prop_overall_g2)
    }
    results_df <- results_df[order(abs(results_df$confidence_diff), decreasing = TRUE), ]
    rownames(results_df) <- NULL
    return(results_df)
}

#' @noRd
.compute_pairwise_effect_sizes <- function(patterns_to_analyze, sequences_group1, sequences_group2, name_group1="G1", name_group2="G2") {
    n_seq_g1 <- length(sequences_group1)
    n_seq_g2 <- length(sequences_group2)
    if (n_seq_g1 == 0 || n_seq_g2 == 0) stop("Both groups must contain sequences for effect size calculation.", call. = FALSE)

    results_df <- data.frame(pattern = patterns_to_analyze, stringsAsFactors = FALSE)
    for(col_name in c("cohens_h", "cohens_d", "cramers_v", "phi_coefficient", "standardized_diff")) {
        results_df[[col_name]] <- numeric(length(patterns_to_analyze))
    }
    
    # Add columns for group-specific counts/props if needed for interpretation, e.g. prop_G1, prop_G2
    col_prop_g1 <- paste0("prop_", name_group1)
    col_prop_g2 <- paste0("prop_", name_group2)
    results_df[[col_prop_g1]] <- numeric(length(patterns_to_analyze))
    results_df[[col_prop_g2]] <- numeric(length(patterns_to_analyze))


    for (i in seq_along(patterns_to_analyze)) {
        current_pattern <- patterns_to_analyze[i]
        count_g1 <- sum(sapply(sequences_group1, function(s) grepl(current_pattern, s, fixed = TRUE)))
        count_g2 <- sum(sapply(sequences_group2, function(s) grepl(current_pattern, s, fixed = TRUE)))

        prop_g1 <- count_g1 / n_seq_g1
        prop_g2 <- count_g2 / n_seq_g2
        results_df[[col_prop_g1]][i] <- prop_g1
        results_df[[col_prop_g2]][i] <- prop_g2

        # Cohen's h
        results_df$cohens_h[i] <- 2 * (asin(sqrt(prop_g1)) - asin(sqrt(prop_g2)))

        # Cohen's d
        pooled_prop <- (count_g1 + count_g2) / (n_seq_g1 + n_seq_g2)
        # Avoid division by zero if pooled_prop is 0 or 1 (pooled_var is 0)
        if (pooled_prop > 1e-9 && pooled_prop < (1 - 1e-9)) {
            pooled_var <- pooled_prop * (1 - pooled_prop)
            results_df$cohens_d[i] <- (prop_g1 - prop_g2) / sqrt(pooled_var)
        } else {
            results_df$cohens_d[i] <- 0
        }

        # Contingency table for Cramer's V and Phi
        # Rows: Group1, Group2; Cols: PatternPresent, PatternAbsent
        # This is for presence/absence of pattern in a *sequence*, not pattern vs other patterns
        # PatternPresent_G1 = count_g1
        # PatternAbsent_G1 = n_seq_g1 - count_g1
        # PatternPresent_G2 = count_g2
        # PatternAbsent_G2 = n_seq_g2 - count_g2
        cont_table <- matrix(c(count_g1, n_seq_g1 - count_g1, count_g2, n_seq_g2 - count_g2),
                             nrow = 2, byrow = TRUE)

        n_total_sequences <- n_seq_g1 + n_seq_g2
        chi_sq_obj <- suppressWarnings(stats::chisq.test(cont_table)) # Suppress warnings for small expected values
        chi_sq_val <- chi_sq_obj$statistic
        
        # Cramer's V (for 2x2 table, V = Phi)
        # Phi = sqrt(chi_sq / N)
        phi_val <- sqrt(chi_sq_val / n_total_sequences)
        results_df$cramers_v[i] <- phi_val
        results_df$phi_coefficient[i] <- phi_val

        # Standardized difference of proportions
        # (p1 - p2) / SE_diff_props
        # SE_diff_props = sqrt( P_pool*(1-P_pool)*(1/n1 + 1/n2) )
        if (pooled_prop > 1e-9 && pooled_prop < (1-1e-9)) {
             pooled_se_diff <- sqrt(pooled_prop * (1 - pooled_prop) * (1/n_seq_g1 + 1/n_seq_g2))
             results_df$standardized_diff[i] <- ifelse(pooled_se_diff == 0, 0, (prop_g1 - prop_g2) / pooled_se_diff)
        } else {
             results_df$standardized_diff[i] <- 0
        }
    }
    results_df <- results_df[order(abs(results_df$cohens_h), decreasing = TRUE), ]
    rownames(results_df) <- NULL
    return(results_df)
}


# ==============================================================================
# MAIN TWO-GROUP PATTERN ANALYSIS FUNCTION (Refactored)
# ==============================================================================

#' Comprehensive Pattern Analysis Between Two Groups (Pairwise Engine)
#'
#' Computes multiple difference measures (support, lift, confidence, effect sizes)
#' for patterns found in sequential data, comparing exactly two groups. This function
#' serves as the engine for pairwise comparisons within `analyze_patterns_multiple`.
#'
#' @param data A data frame with sequences in wide format.
#' @param group_col A character string specifying the name of the column in `data`
#'        that contains group assignments.
#' @param min_length Minimum pattern length to analyze (default: 2).
#' @param max_length Maximum pattern length to analyze (default: 5).
#' @param min_frequency Minimum frequency required for a pattern to be included in analysis (default: 2).
#' @param measures A character vector specifying which measures to compute.
#'        Options: "support", "lift", "confidence", "effect_size". Default is all.
#' @return A list object of class `pattern_analysis`. Each element in the list
#'         corresponds to a computed measure (e.g., `results$support`, `results$lift`)
#'         and contains a data frame of patterns and their associated statistics.
#'         The list also includes a `metadata` element with details about the analysis.
#'
#' @seealso \code{\link{analyze_patterns_multiple}} for analyzing more than two groups.
 #' @family Pattern Analysis Functions
#' @export
analyze_patterns <- function(data, group_col = "Group", min_length = 2, max_length = 5,
                             min_frequency = 2,
                             measures = c("support", "lift", "confidence", "effect_size")) {

  # --- Input Validation ---
  if (!is.data.frame(data)) stop("'data' must be a data frame.", call. = FALSE)
  if (nrow(data) == 0) stop("'data' cannot be empty.", call. = FALSE)
  if (!is.character(group_col) || length(group_col) != 1 || !group_col %in% names(data)) {
    stop("'group_col' must be a single valid column name in 'data'.", call. = FALSE)
  }
  # Other parameter validations (min_length, max_length, etc.)
  if(min_length > max_length) stop("min_length cannot be greater than max_length.", call.=FALSE)
  if(min_frequency < 1) stop("min_frequency must be at least 1.", call.=FALSE)
  
  valid_measure_options <- c("support", "lift", "confidence", "effect_size")
  if(!all(measures %in% valid_measure_options)){
      stop(paste0("Invalid measure(s) specified. Choose from: ", paste(valid_measure_options, collapse=", ")), call.=FALSE)
  }

  # --- Data Preparation ---
  cat("Preparing sequence data for two-group analysis...\n")
  
  # Identify data columns (all columns except the group column)
  data_cols <- setdiff(names(data), group_col)
  if(length(data_cols) == 0) stop("No data columns found after excluding the group column.", call.=FALSE)
  data_cols_indices <- match(data_cols, names(data))

  # Use the generalized preprocessing function from utils.R
  # This function handles conversion of rows to hyphenated sequences and filters by min_length_seq
  # Here, min_length_seq refers to the length of the original sequence, not n-gram pattern length.
  # Let's assume min_length for analyze_patterns refers to n-gram length, so sequence min_length is 1.
  preprocessed_output <- preprocess_sequences_for_analysis(data, group_col, data_cols_indices, min_length_seq = 1)

  sequences_by_group <- preprocessed_output$sequences_by_group

  # Validate that exactly two groups are present for this function
  group_names_vec <- names(sequences_by_group)
  if (length(group_names_vec) != 2) {
    stop(paste0("`analyze_patterns` is designed for exactly 2 groups. Found ", length(group_names_vec),
                " (", paste(group_names_vec, collapse=", "), "). ",
                "For multi-group analysis, use `analyze_patterns_multiple`."), call. = FALSE)
  }

  group1_name <- group_names_vec[1]
  group2_name <- group_names_vec[2]
  sequences_g1 <- sequences_by_group[[group1_name]]
  sequences_g2 <- sequences_by_group[[group2_name]]

  cat(sprintf("Group %s: %d sequences\n", group1_name, length(sequences_g1)))
  cat(sprintf("Group %s: %d sequences\n", group2_name, length(sequences_g2)))

  if(length(sequences_g1) == 0 || length(sequences_g2) == 0){
      stop("One or both groups have no valid sequences after preprocessing.", call.=FALSE)
  }

  # --- Pattern Extraction ---
  cat("Extracting n-gram patterns...\n")
  all_extracted_patterns <- character(0)
  # Iterate n-gram lengths from min_length to max_length
  for (n_gram_len in min_length:max_length) {
    # Extract n-grams from all valid sequences (across both groups)
    for (seq_str in c(sequences_g1, sequences_g2)) {
      if (nzchar(seq_str)) { # Ensure sequence string is not empty
        ngrams_from_seq <- .extract_ngrams_from_sequence(seq_str, n_gram_len)
        if (length(ngrams_from_seq) > 0) {
          all_extracted_patterns <- c(all_extracted_patterns, ngrams_from_seq)
        }
      }
    }
  }
  
  if (length(all_extracted_patterns) == 0) {
    stop(paste0("No n-gram patterns found for lengths ", min_length, " to ", max_length, "."), call. = FALSE)
  }

  # Filter patterns by minimum frequency (overall frequency across both groups)
  pattern_freq_table <- table(all_extracted_patterns)
  frequent_patterns_to_analyze <- names(pattern_freq_table[pattern_freq_table >= min_frequency])

  if (length(frequent_patterns_to_analyze) == 0) {
    stop(paste0("No patterns meet the minimum frequency threshold of ", min_frequency, "."), call. = FALSE)
  }
  cat(sprintf("Found %d unique patterns meeting frequency threshold.\n", length(frequent_patterns_to_analyze)))

  # --- Compute Measures ---
  analysis_results <- list() # Renamed

  if ("support" %in% measures) {
    cat("Computing support measures...\n")
    analysis_results$support <- .compute_pairwise_support(frequent_patterns_to_analyze, sequences_g1, sequences_g2, group1_name, group2_name)
  }
  if ("lift" %in% measures) {
    cat("Computing lift measures...\n")
    analysis_results$lift <- .compute_pairwise_lift(frequent_patterns_to_analyze, sequences_g1, sequences_g2, group1_name, group2_name)
  }
  if ("confidence" %in% measures) {
    cat("Computing confidence measures...\n")
    analysis_results$confidence <- .compute_pairwise_confidence(frequent_patterns_to_analyze, sequences_g1, sequences_g2, group1_name, group2_name)
  }
  if ("effect_size" %in% measures) {
    cat("Computing effect size measures...\n")
    analysis_results$effect_size <- .compute_pairwise_effect_sizes(frequent_patterns_to_analyze, sequences_g1, sequences_g2, group1_name, group2_name)
  }

  # --- Metadata ---
  analysis_results$metadata <- list(
    group_names = group_names_vec,
    n_sequences_per_group = stats::setNames(c(length(sequences_g1), length(sequences_g2)), group_names_vec),
    n_patterns_analyzed = length(frequent_patterns_to_analyze),
    pattern_length_range = c(min_length, max_length),
    min_pattern_frequency = min_frequency
  )

  cat("Pattern analysis complete for the two groups.\n")
  class(analysis_results) <- "pattern_analysis"
  return(analysis_results)
}


# ==============================================================================
# PRINT AND SUMMARY METHODS (for pattern_analysis object)
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

# ==============================================================================
# MULTI-GROUP PATTERN ANALYSIS FUNCTION
# ==============================================================================

#' Comprehensive Pattern Analysis Across Multiple Groups
#'
#' Performs pattern analysis (support, lift, confidence, effect sizes) for every
#' unique pair of groups when more than two groups are present in the data.
#' It utilizes the `analyze_patterns` function for each pairwise comparison.
#'
#' @param data A data frame with sequences in wide format. Each row represents a
#'   sequence, and columns represent time points or states. Must include a group column.
#' @param group_col A character string specifying the name of the column in `data`
#'   that contains group assignments.
#' @param min_length Minimum pattern (n-gram) length to analyze (default: 2).
#' @param max_length Maximum pattern (n-gram) length to analyze (default: 5).
#' @param min_frequency Minimum overall frequency required for a pattern to be
#'   included in the analysis (default: 2).
#' @param measures A character vector specifying which measures to compute for each
#'   pairwise comparison. Options: "support", "lift", "confidence", "effect_size".
#'   Default is all.
#'
#' @return An object of class `analyze_patterns_multiple_result`. This object is a
#'   list where each element is the result of an `analyze_patterns` call for a
#'   unique pair of groups. The list elements are named in the format "GroupA_vs_GroupB".
#'   If only two groups are provided in the input data, it effectively calls
#'   `analyze_patterns` once, and its result is wrapped in a list for consistency.
#'
#' @family Pattern Analysis Functions
#' @seealso \code{\link{analyze_patterns}} for the underlying pairwise analysis logic.
#' @export
#' @examples
#' \dontrun{
#'   # Create dummy data with 3 groups
#'   pattern_data_multi <- data.frame(
#'     T1 = sample(letters[1:4], 45, replace = TRUE),
#'     T2 = sample(letters[1:4], 45, replace = TRUE),
#'     T3 = sample(letters[1:4], 45, replace = TRUE),
#'     Category = sample(c("Alpha", "Beta", "Gamma"), 45, replace = TRUE)
#'   )
#'
#'   # Perform multi-group pattern analysis
#'   multi_pattern_results <- analyze_patterns_multiple(
#'     pattern_data_multi,
#'     group_col = "Category"
#'   )
#'   print(multi_pattern_results)
#'   summary(multi_pattern_results)
#'
#'   # Access results for a specific pair, e.g., Alpha vs Beta:
#'   # results_alpha_vs_beta <- multi_pattern_results[["Alpha_vs_Beta"]]
#'   # summary(results_alpha_vs_beta, measure="support")
#' }
analyze_patterns_multiple <- function(data, group_col = "Group", min_length = 2, max_length = 5,
                                    min_frequency = 2,
                                    measures = c("support", "lift", "confidence", "effect_size")) {

  # --- Input Validation ---
  if (!is.data.frame(data)) stop("'data' must be a data frame.", call. = FALSE)
  if (nrow(data) == 0) stop("'data' cannot be empty.", call. = FALSE)
  if (!is.character(group_col) || length(group_col) != 1 || !group_col %in% names(data)) {
    stop("'group_col' must be a single valid column name in 'data'.", call. = FALSE)
  }
  # Further validation for min_length, max_length, etc., is handled by analyze_patterns

  # --- Data Preparation & Group Identification ---
  # Validate the group column and get unique group levels
  # Use a simplified validation here just to get group levels, as analyze_patterns will do more.
  group_vector <- data[[group_col]]
  if(!is.factor(group_vector)) group_vector <- as.factor(group_vector)

  if (any(is.na(group_vector))) {
    stop("Group column '", group_col, "' contains NA values. Please remove or impute them.", call. = FALSE)
  }
  if (any(levels(group_vector) == "")){
      stop("Group column '", group_col, "' contains empty string as a group level. Please provide valid group names.", call. = FALSE)
  }

  actual_group_levels <- levels(droplevels(group_vector)) # Get current levels in data
  num_groups <- length(actual_group_levels)

  if (num_groups < 2) {
    stop("At least two groups are required for pattern analysis. Found ", num_groups, " groups.", call. = FALSE)
  }

  all_pairwise_pattern_results <- list()

  if (num_groups == 2) {
    cat("Exactly two groups found. Performing a single pairwise pattern analysis...\n")
    # Call analyze_patterns directly with the full dataset and group_col
    pair_result <- analyze_patterns(
      data = data,
      group_col = group_col,
      min_length = min_length, max_length = max_length,
      min_frequency = min_frequency, measures = measures
    )
    pair_name <- paste(actual_group_levels[1], "vs", actual_group_levels[2], sep = "_")
    all_pairwise_pattern_results[[pair_name]] <- pair_result
  } else {
    cat(sprintf("%d groups found. Performing pattern analysis for all unique pairs...\n", num_groups))
    group_pairs <- generate_group_pairs(actual_group_levels) # from utils.R

    for (pair in group_pairs) {
      group1_name <- pair[1]
      group2_name <- pair[2]
      pair_name <- paste(group1_name, "vs", group2_name, sep = "_")
      cat(sprintf("\n--- Analyzing patterns for: %s ---\n", pair_name))

      # Subset the original data to include only the current pair of groups
      current_pair_data <- data[data[[group_col]] %in% c(group1_name, group2_name), , drop = FALSE]
      # Ensure the group column in the subsetted data only contains the two current group levels
      current_pair_data[[group_col]] <- factor(current_pair_data[[group_col]], levels = c(group1_name, group2_name))


      pair_result <- tryCatch({
        analyze_patterns(
          data = current_pair_data,
          group_col = group_col, # Pass the original group column name
          min_length = min_length, max_length = max_length,
          min_frequency = min_frequency, measures = measures
        )
      }, error = function(e) {
        warning(sprintf("Error during pattern analysis of %s: %s. Skipping this pair.", pair_name, e$message), call. = FALSE)
        return(NULL) # Store NULL if a pair fails
      })

      all_pairwise_pattern_results[[pair_name]] <- pair_result
    }
  }

  class(all_pairwise_pattern_results) <- "analyze_patterns_multiple_result"
  cat("\nMulti-group pattern analysis complete.\n")
  return(all_pairwise_pattern_results)
}


#' Print Method for analyze_patterns_multiple_result Objects
#'
#' Provides a concise overview of the multi-group pattern analysis results.
#' It lists each pairwise analysis performed and a high-level summary for each.
#'
#' @param x An object of class `analyze_patterns_multiple_result`.
#' @param ... Additional arguments (unused).
#' @return Invisibly returns the input object `x`.
#' @export
print.analyze_patterns_multiple_result <- function(x, ...) {
  cat("Multi-Group Pattern Analysis Results\n")
  cat("====================================\n")
  num_comparisons <- length(x)
  cat("Number of pairwise pattern analyses performed:", num_comparisons, "\n\n")

  if (num_comparisons == 0) {
    cat("No analyses were performed or results available.\n")
    return(invisible(x))
  }

  for (i in seq_along(x)) {
    pair_name <- names(x)[i]
    result_item <- x[[i]] # This is a 'pattern_analysis' object

    cat(paste0("--- Analysis for: ", pair_name, " ---\n"))
    if (is.null(result_item)) {
        cat("  Result for this pair is NULL (an error likely occurred).\n\n")
    } else if (inherits(result_item, "pattern_analysis")) {
      cat("  Groups compared:", paste(result_item$metadata$group_names, collapse = " vs "), "\n")
      cat("  Patterns analyzed meeting frequency criteria:", result_item$metadata$n_patterns_analyzed, "\n")
      cat("  Available measures computed:", paste(names(result_item)[!names(result_item) %in% "metadata"], collapse = ", "), "\n")

      # Glimpse of top support pattern if available
      if("support" %in% names(result_item) && is.data.frame(result_item$support) && nrow(result_item$support) > 0){
          top_support_pattern <- result_item$support$pattern[1]
          cat("  Top pattern by support difference: ", top_support_pattern, "\n")
      }
      cat("\n")
    } else {
      cat("  Result for this pair is not a 'pattern_analysis' object.\n\n")
    }
  }
  cat("Use summary(object) for more details or access individual pairwise results using [[name]].\n")
  invisible(x)
}

#' Summary Method for analyze_patterns_multiple_result Objects
#'
#' @param object An `analyze_patterns_multiple_result` object.
#' @param measure_to_summarize For each pairwise comparison, which specific measure
#'   (e.g., "support", "lift") should be summarized if `show_detailed_summaries` is TRUE.
#'   If NULL, the first available measure is summarized.
#' @param top_n_detailed For detailed summaries, how many top patterns to show per measure (passed to `summary.pattern_analysis`).
#' @param show_detailed_summaries Logical, whether to call `summary()` on each
#'   individual `pattern_analysis` object within the results (default TRUE).
#' @param ... Additional arguments passed to `summary.pattern_analysis` for each pairwise summary.
#' @return Invisibly returns the input object `object`.
#' @export
summary.analyze_patterns_multiple_result <- function(object,
                                                   measure_to_summarize = NULL,
                                                   top_n_detailed = 5,
                                                   show_detailed_summaries = TRUE,
                                                   ...) {
  cat("Summary of Multi-Group Pattern Analysis\n")
  cat("=======================================\n")
  num_comparisons <- length(object)
  cat("Total pairwise analyses performed:", num_comparisons, "\n")

  if (num_comparisons == 0) {
    cat("No analyses to summarize.\n")
    return(invisible(object))
  }

  # Overall summary (e.g., list of pairs, number of patterns in each)
  cat("\nOverview of Pairwise Analyses:\n")
  for (pair_name in names(object)) {
    result_item <- object[[pair_name]]
    if (!is.null(result_item) && inherits(result_item, "pattern_analysis")) {
      cat(sprintf("  - %s: %d patterns analyzed, measures computed: %s\n",
                  pair_name,
                  result_item$metadata$n_patterns_analyzed,
                  paste(names(result_item)[!names(result_item) %in% "metadata"], collapse = ", ")))
    } else {
      cat(sprintf("  - %s: Analysis failed or result is invalid.\n", pair_name))
    }
  }

  if (show_detailed_summaries) {
    cat("\n\n--- Detailed Summaries for Each Pairwise Analysis ---\n")
    for (pair_name in names(object)) {
      cat(paste0("\n--- Summary for: ", pair_name, " ---\n"))
      result_item <- object[[pair_name]]
      if (!is.null(result_item) && inherits(result_item, "pattern_analysis")) {
        # Call summary.pattern_analysis for this item
        # If measure_to_summarize is NULL, summary.pattern_analysis picks the first one.
        summary(result_item, measure = measure_to_summarize, top_n = top_n_detailed)
      } else {
        cat("  No detailed summary available (result was NULL or invalid).\n")
      }
    }
  } else {
    cat("\nSet show_detailed_summaries = TRUE to see detailed tables for each pair.\n")
  }
  invisible(object)
}
# ==============================================================================
# SEQUENCE COMPARISON TOOLKIT
# ==============================================================================
#
# Advanced sequence comparison functions with statistical testing capabilities.
# Includes chi-square tests, Fisher's exact tests, and comprehensive pattern
# analysis with visualization support.
#
# ==============================================================================

# Ensure utils are available if this file is sourced directly for testing,
# though NAMESPACE should handle this in package context.
# No direct source() calls needed if using standard package structure.

# ==============================================================================
# INTERNAL HELPER FUNCTIONS FOR SEQUENCE COMPARISON
# (Moved from being nested inside compare_sequences)
# ==============================================================================

#' @noRd
.extract_subsequences_from_df <- function(sequences_df, n_gram_length) {
  if (n_gram_length < 1) return(character(0))
  if (nrow(sequences_df) == 0) return(character(0))

  all_subsequences <- character(0)
  for (i in seq_len(nrow(sequences_df))) {
    seq_row_vector <- unlist(sequences_df[i, ])
    # Filter out NA, empty strings, and single characters (potential group markers if not handled upstream)
    # Assuming upstream preprocessing (e.g. in compare_sequences or its callers) handles NA conversion
    # and that data passed here is sequence states.
    clean_seq_vector <- seq_row_vector[!is.na(seq_row_vector) & nzchar(as.character(seq_row_vector))]
    # Optional: Add check for nchar > 1 if single characters are not valid states
    # clean_seq_vector <- clean_seq_vector[nchar(as.character(clean_seq_vector)) > 1]


    if (length(clean_seq_vector) >= n_gram_length) {
      for (j in seq_len(length(clean_seq_vector) - n_gram_length + 1)) {
        subseq <- paste(clean_seq_vector[j:(j + n_gram_length - 1)], collapse = "-")
        all_subsequences <- c(all_subsequences, subseq)
      }
    }
  }
  return(all_subsequences)
}

#' @noRd
.analyze_pairwise_discrimination <- function(seq_data_group1, seq_data_group2, n_gram_length, group1_name = "Group1", group2_name = "Group2") {
  subseq_g1 <- .extract_subsequences_from_df(seq_data_group1, n_gram_length)
  subseq_g2 <- .extract_subsequences_from_df(seq_data_group2, n_gram_length)

  if (length(subseq_g1) == 0 && length(subseq_g2) == 0) {
    return(list(
      n = n_gram_length,
      patterns = data.frame(
        pattern = character(0),
        freq_group1 = numeric(0),
        freq_group2 = numeric(0),
        prop_group1 = numeric(0),
        prop_group2 = numeric(0),
        prop_diff = numeric(0),
        discrimination = numeric(0),
        stringsAsFactors = FALSE
      ),
      total_group1 = 0, total_group2 = 0, n_patterns = 0
    ))
  }

  freq_g1 <- table(subseq_g1)
  freq_g2 <- table(subseq_g2)
  all_patterns_char <- unique(c(names(freq_g1), names(freq_g2)))

  patterns_df <- data.frame(
    pattern = all_patterns_char,
    freq_group1 = as.numeric(freq_g1[all_patterns_char]),
    freq_group2 = as.numeric(freq_g2[all_patterns_char]),
    stringsAsFactors = FALSE
  )
  patterns_df[is.na(patterns_df)] <- 0
  names(patterns_df) <- c("pattern", paste0("freq_", group1_name), paste0("freq_", group2_name))


  total_g1 <- sum(patterns_df[[paste0("freq_", group1_name)]])
  total_g2 <- sum(patterns_df[[paste0("freq_", group2_name)]])

  col_prop_g1 <- paste0("prop_", group1_name)
  col_prop_g2 <- paste0("prop_", group2_name)

  if (total_g1 > 0 && total_g2 > 0) {
    patterns_df[[col_prop_g1]] <- patterns_df[[paste0("freq_", group1_name)]] / total_g1
    patterns_df[[col_prop_g2]] <- patterns_df[[paste0("freq_", group2_name)]] / total_g2
    patterns_df$prop_diff <- patterns_df[[col_prop_g1]] - patterns_df[[col_prop_g2]]
    patterns_df$discrimination <- abs(patterns_df$prop_diff) * sqrt(patterns_df[[paste0("freq_", group1_name)]] + patterns_df[[paste0("freq_", group2_name)]])
    patterns_df <- patterns_df[order(patterns_df$discrimination, decreasing = TRUE), ]
  } else {
    patterns_df[[col_prop_g1]] <- 0
    patterns_df[[col_prop_g2]] <- 0
    patterns_df$prop_diff <- 0
    patterns_df$discrimination <- 0
  }

  # Standardize output column names for freq and prop if needed by downstream generic functions
  # For now, keeping them group-specific.

  return(list(
    n = n_gram_length,
    patterns = patterns_df,
    total_group1 = total_g1,
    total_group2 = total_g2,
    n_patterns = nrow(patterns_df)
  ))
}

#' @noRd
.analyze_pairwise_statistical <- function(seq_data_group1, seq_data_group2, n_gram_length,
                                       correction_method, test_method, min_expected_val,
                                       group1_name = "Group1", group2_name = "Group2") {

  subseq_g1 <- .extract_subsequences_from_df(seq_data_group1, n_gram_length)
  subseq_g2 <- .extract_subsequences_from_df(seq_data_group2, n_gram_length)

  if (length(subseq_g1) == 0 && length(subseq_g2) == 0) {
    return(list(n = n_gram_length, patterns = data.frame(), n_patterns = 0, n_significant = 0))
  }

  freq_g1 <- table(subseq_g1)
  freq_g2 <- table(subseq_g2)
  all_patterns_char <- unique(c(names(freq_g1), names(freq_g2)))

  results_df <- data.frame(
    pattern = all_patterns_char,
    freq_group1 = as.numeric(freq_g1[all_patterns_char]),
    freq_group2 = as.numeric(freq_g2[all_patterns_char]),
    stringsAsFactors = FALSE
  )
  results_df[is.na(results_df)] <- 0
  names(results_df) <- c("pattern", paste0("freq_", group1_name), paste0("freq_", group2_name))

  results_df$p_value <- NA_real_
  results_df$test_used <- NA_character_
  results_df$statistic <- NA_real_

  sum_freq_g1 <- sum(results_df[[paste0("freq_", group1_name)]])
  sum_freq_g2 <- sum(results_df[[paste0("freq_", group2_name)]])

  for (i in seq_len(nrow(results_df))) {
    # Create contingency table:
    #          PatternPresent PatternAbsent
    # Group1   a              b
    # Group2   c              d
    # a = freq of pattern in Group1
    # b = total patterns in Group1 - freq of pattern in Group1
    # c = freq of pattern in Group2
    # d = total patterns in Group2 - freq of pattern in Group2

    # This interpretation of contingency table (pattern vs other patterns within group)
    # might be different from (pattern presence in sequence vs absence in sequence).
    # The original code implies the former: sum(results$freq_A) - results$freq_A[i]
    # This compares one pattern against all *other* patterns' occurrences.
    # If the goal is to compare presence/absence of a pattern in *sequences*, the table setup is different.
    # Current setup:
    #   freq_pattern_g1, sum_other_patterns_g1
    #   freq_pattern_g2, sum_other_patterns_g2

    # Let's stick to the original interpretation for now:
    # This compares the proportion of a specific pattern relative to all observed patterns in one group vs. another.
    a <- results_df[[paste0("freq_", group1_name)]][i]
    c <- results_df[[paste0("freq_", group2_name)]][i]
    b <- sum_freq_g1 - a
    d <- sum_freq_g2 - c

    cont_table <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)

    # Use the safe_statistical_test helper
    test_output <- safe_statistical_test(cont_table, test_method, min_expected_val)
    results_df$p_value[i] <- test_output$p_value
    results_df$test_used[i] <- test_output$test_used
    results_df$statistic[i] <- test_output$statistic
  }

  if (correction_method != "none" && sum(!is.na(results_df$p_value)) > 0) {
    results_df$p_adjusted <- stats::p.adjust(results_df$p_value, method = correction_method)
  } else {
    results_df$p_adjusted <- results_df$p_value
  }

  results_df$significant <- !is.na(results_df$p_adjusted) & results_df$p_adjusted < 0.05
  results_df <- results_df[order(results_df$significant, results_df$p_value, decreasing = c(TRUE, FALSE)), ]

  return(list(
    n = n_gram_length,
    patterns = results_df,
    n_patterns = nrow(results_df),
    n_significant = sum(results_df$significant, na.rm = TRUE)
  ))
}


# ==============================================================================
# MAIN TWO-GROUP SEQUENCE COMPARISON FUNCTION (Refactored Internals)
# ==============================================================================

#' Compare Sequences Between Two Groups (Pairwise Engine)
#'
#' Performs comprehensive sequence comparison analysis between exactly two groups,
#' identifying discriminating subsequences and their statistical associations.
#' This function is the engine for pairwise comparisons and is called by
#' `compare_sequences_multiple` for multi-group scenarios.
#'
#' @param data A data frame containing sequence data in wide format.
#' @param group A character string specifying the name of the group column in `data`,
#'        or a factor/vector of group assignments.
#' @param min_length Minimum subsequence length to analyze (default: 2).
#' @param max_length Maximum subsequence length to analyze (default: 5).
#' @param top_n Number of top patterns to return and display (default: 10).
#' @param detailed Logical, whether to return detailed results by subsequence length.
#' @param statistical Logical, whether to perform statistical testing.
#' @param correction Character, multiple comparison correction method.
#' @param test_method Character, statistical test selection ("auto", "fisher", "chi.squared").
#' @param min_expected Numeric, minimum expected count for Chi-squared.
#'
#' @return A `compare_sequences` object.
#'
#' @seealso \code{\link{compare_sequences_multiple}} for comparing more than two groups.
#' @family Sequence Comparison Functions
#' @export
compare_sequences <- function(data, group, min_length = 2, max_length = 5, top_n = 10,
                             detailed = FALSE, statistical = FALSE, correction = "bonferroni",
                             test_method = "auto", min_expected = 5) {

  # --- Input Validation ---
  if (!is.data.frame(data)) stop("'data' must be a data frame.", call. = FALSE)
  if (nrow(data) == 0) stop("'data' cannot be empty.", call. = FALSE)
  # Basic parameter checks
  if (!is.numeric(min_length) || !is.numeric(max_length) || min_length < 1 || max_length < min_length) {
    stop("'min_length' and 'max_length' must be positive integers with min_length <= max_length.", call. = FALSE)
  }
  if (!is.numeric(top_n) || top_n < 1) stop("'top_n' must be a positive integer.", call. = FALSE)
  if (!is.logical(detailed) || !is.logical(statistical)) stop("'detailed' and 'statistical' must be logical values.", call. = FALSE)
  if (!correction %in% c("bonferroni", "holm", "hochberg", "BH", "BY", "none")) {
    stop("'correction' must be one of: bonferroni, holm, hochberg, BH, BY, none.", call. = FALSE)
  }
  if (!test_method %in% c("auto", "fisher", "chi.squared")) {
    stop("'test_method' must be one of: auto, fisher, chi.squared.", call. = FALSE)
  }


  # --- Data Preprocessing ---
  group_col_name <- NULL
  sequence_df <- NULL

  if (is.character(group) && length(group) == 1 && group %in% names(data)) {
    group_col_name <- group
    group_vector <- data[[group_col_name]]
    sequence_df <- data[, !names(data) %in% group_col_name, drop = FALSE]
  } else if ((is.vector(group) || is.factor(group)) && length(group) == nrow(data)) {
    group_vector <- group
    # Assume `data` does not contain the group column if `group` is provided as a vector
    # Or, if it might, we need a way to exclude it if user also names it.
    # For now, assume 'data' is purely sequence columns if 'group' is a vector.
    sequence_df <- data
  } else {
    stop("'group' must be a valid column name in 'data', or a vector/factor of the same length as nrow(data).", call. = FALSE)
  }

  # Validate group vector (min_groups = 2 for this specific function)
  group_factor <- validate_group_input(group_vector, nrow(sequence_df), min_groups = 2)

  actual_group_levels <- levels(group_factor)
  if (length(actual_group_levels) != 2) {
    stop(paste0("`compare_sequences` is designed for exactly 2 groups. Found ", length(actual_group_levels),
                " (", paste(actual_group_levels, collapse=", "), "). ",
                "For multi-group comparisons, use `compare_sequences_multiple`."), call. = FALSE)
  }

  # Split data into two groups
  # Ensure sequence_df contains only sequence columns
  seq_A_df <- sequence_df[group_factor == actual_group_levels[1], , drop = FALSE]
  seq_B_df <- sequence_df[group_factor == actual_group_levels[2], , drop = FALSE]

  # Clean data: remove empty rows and convert empty strings to NA
  seq_A_df <- seq_A_df[apply(seq_A_df, 1, function(r) !all(is.na(r) | nzchar(as.character(r))==FALSE ) ), , drop = FALSE]
  seq_B_df <- seq_B_df[apply(seq_B_df, 1, function(r) !all(is.na(r) | nzchar(as.character(r))==FALSE ) ), , drop = FALSE]

  if (nrow(seq_A_df) == 0 || nrow(seq_B_df) == 0) {
    stop("Both groups must have at least one non-empty sequence after cleaning.", call. = FALSE)
  }
  # seq_A_df[seq_A_df == ""] <- NA # Not strictly necessary if nzchar used above
  # seq_B_df[seq_B_df == ""] <- NA


  # --- Main Analysis Loop ---
  cat("Analyzing sequences between groups:", actual_group_levels[1], "and", actual_group_levels[2], "\n")
  cat("Group sizes:", nrow(seq_A_df), "vs", nrow(seq_B_df), "sequences\n")
  cat("Subsequence lengths:", min_length, "to", max_length, "\n\n")

  analysis_results_list <- list() # Renamed
  for (n_len in min_length:max_length) {
    if (statistical) {
      analysis_results_list[[paste0("length_", n_len)]] <- .analyze_pairwise_statistical(
        seq_A_df, seq_B_df, n_len, correction, test_method, min_expected,
        group1_name = actual_group_levels[1], group2_name = actual_group_levels[2]
      )
    } else {
      analysis_results_list[[paste0("length_", n_len)]] <- .analyze_pairwise_discrimination(
        seq_A_df, seq_B_df, n_len,
        group1_name = actual_group_levels[1], group2_name = actual_group_levels[2]
      )
    }
  }

  # --- Prepare Output ---
  summary_output_data <- NULL # Renamed
  if (detailed) {
    summary_tables_list <- list() # Renamed
    for (n_len in min_length:max_length) {
      result_item <- analysis_results_list[[paste0("length_", n_len)]]
      if (is.data.frame(result_item$patterns) && nrow(result_item$patterns) > 0) {
        summary_tables_list[[paste0("length_", n_len)]] <- utils::head(result_item$patterns, top_n)
      } else {
        summary_tables_list[[paste0("length_", n_len)]] <- data.frame() # Empty df for consistency
      }
    }
    summary_output_data <- summary_tables_list
  } else {
    all_patterns_df <- data.frame() # Renamed
    for (n_len in min_length:max_length) {
      result_item <- analysis_results_list[[paste0("length_", n_len)]]
      if (is.data.frame(result_item$patterns) && nrow(result_item$patterns) > 0) {
        patterns_with_len <- cbind(length = n_len, result_item$patterns) # Renamed
        all_patterns_df <- rbind(all_patterns_df, patterns_with_len)
      }
    }

    if (nrow(all_patterns_df) > 0) {
      if (statistical && "significant" %in% names(all_patterns_df) && "p_value" %in% names(all_patterns_df)) {
        all_patterns_df <- all_patterns_df[order(all_patterns_df$significant, all_patterns_df$p_value, decreasing = c(TRUE, FALSE)), ]
      } else if (!statistical && "discrimination" %in% names(all_patterns_df)) {
        all_patterns_df <- all_patterns_df[order(all_patterns_df$discrimination, decreasing = TRUE), ]
      } # else, no specific sort if columns missing, head will just take top rows
      summary_output_data <- utils::head(all_patterns_df, top_n)
    } else {
      summary_output_data <- data.frame()
    }
  }

  total_patterns_num <- sum(sapply(analysis_results_list, function(x) ifelse(is.null(x$n_patterns), 0, x$n_patterns))) # Renamed
  summary_stats_list <- NULL # Renamed
  if (statistical) {
    total_significant_num <- sum(sapply(analysis_results_list, function(x) ifelse(is.null(x$n_significant), 0, x$n_significant))) # Renamed
    summary_stats_list <- list(
      total_patterns = total_patterns_num,
      total_significant = total_significant_num,
      significance_rate = if (total_patterns_num > 0) total_significant_num / total_patterns_num else 0
    )
  } else {
    summary_stats_list <- list(
      total_patterns = total_patterns_num,
      # Naming group sizes with actual group names for clarity
      group_sizes = stats::setNames(c(nrow(seq_A_df), nrow(seq_B_df)), actual_group_levels)
    )
  }

  # --- Create Visualizations (if applicable) ---
  # The create_plots function and its helper create_heatmap are kept here for the two-group case.
  # compare_sequences_multiple will need its own strategy for visualization.

  # Local function for plotting, specific to compare_sequences output structure
  .create_cs_plots <- function(summary_plot_data, parameters_plot, groups_plot_names) {

    .create_cs_heatmap <- function(patterns_heat_df, groups_heat_names, top_n_heat, is_statistical_heat) {
      if (nrow(patterns_heat_df) == 0) return()
      plot_data_df <- utils::head(patterns_heat_df, top_n_heat) # Renamed

      col_freq_g1 <- paste0("freq_", groups_heat_names[1])
      col_freq_g2 <- paste0("freq_", groups_heat_names[2])
      col_prop_g1 <- paste0("prop_", groups_heat_names[1])
      col_prop_g2 <- paste0("prop_", groups_heat_names[2])

      residual_matrix <- NULL

      if (is_statistical_heat) {
        if (!(col_freq_g1 %in% names(plot_data_df)) || !(col_freq_g2 %in% names(plot_data_df))) {
          warning("Frequency columns (e.g., freq_GroupA) not found for heatmap. Skipping.", call. = FALSE)
          return()
        }
        total_g1 <- sum(plot_data_df[[col_freq_g1]], na.rm = TRUE)
        total_g2 <- sum(plot_data_df[[col_freq_g2]], na.rm = TRUE)

        if (total_g1 == 0 || total_g2 == 0) {
          warning("One group has zero total frequency for top patterns in plot. Heatmap may fail or be misleading.", call. = FALSE)
          if ("prop_diff" %in% names(plot_data_df)) {
            prop_d <- plot_data_df$prop_diff # Renamed
            residual_matrix <- cbind(prop_d, -prop_d)
          } else return()
        } else {
          prop_g1_vals <- plot_data_df[[col_freq_g1]] / total_g1 # Renamed
          prop_g2_vals <- plot_data_df[[col_freq_g2]] / total_g2 # Renamed
          overall_props <- (plot_data_df[[col_freq_g1]] + plot_data_df[[col_freq_g2]]) / (total_g1 + total_g2) # Renamed

          se_g1 <- sqrt(overall_props * (1 - overall_props) / total_g1) # Renamed
          se_g2 <- sqrt(overall_props * (1 - overall_props) / total_g2) # Renamed

          resid_g1 <- ifelse(se_g1 == 0, 0, (prop_g1_vals - overall_props) / se_g1) # Renamed
          resid_g2 <- ifelse(se_g2 == 0, 0, (prop_g2_vals - overall_props) / se_g2) # Renamed
          residual_matrix <- cbind(resid_g1, resid_g2)
        }
      } else { # Discrimination
        if (!("prop_diff" %in% names(plot_data_df))) {
          warning("prop_diff not found for discrimination heatmap. Skipping.", call. = FALSE)
          return()
        }
        prop_d <- plot_data_df$prop_diff # Renamed
        residual_matrix <- cbind(prop_d, -prop_d)
      }

      if(is.null(residual_matrix) || nrow(residual_matrix) == 0) return()

      rownames(residual_matrix) <- plot_data_df$pattern
      colnames(residual_matrix) <- groups_heat_names

      original_layout <- grDevices::layout(matrix(c(1, 2), nrow = 1), widths = c(4, 1))
      on.exit(grDevices::layout(original_layout), add = TRUE) # Ensure layout is reset
      graphics::par(mar = c(5, 12, 4, 1))

      max_abs_val <- max(abs(residual_matrix), na.rm = TRUE) # Renamed
      if (!is.finite(max_abs_val) || max_abs_val == 0) max_abs_val <- 1

      color_palette <- grDevices::colorRampPalette(c("blue", "white", "red"))(100) # Renamed

      graphics::image(1:ncol(residual_matrix), 1:nrow(residual_matrix), t(residual_matrix),
                      col = color_palette, breaks = seq(-max_abs_val, max_abs_val, length.out = 101),
                      xlab = "Groups", ylab = "", axes = FALSE)
      graphics::axis(1, at = 1:ncol(residual_matrix), labels = colnames(residual_matrix))
      graphics::axis(2, at = 1:nrow(residual_matrix), labels = rownames(residual_matrix), las = 2)

      graphics::par(mar = c(5, 1, 4, 3))
      legend_y_coords <- seq(0, 1, length.out = 100) # Renamed
      legend_color_matrix <- matrix(legend_y_coords, ncol = 1) # Renamed
      graphics::image(1, legend_y_coords, t(legend_color_matrix), col = color_palette, axes = FALSE, xlab = "", ylab = "")
      legend_axis_labels <- seq(-max_abs_val, max_abs_val, length.out = 5) # Renamed
      graphics::axis(4, at = seq(0, 1, length.out = 5), labels = sprintf("%.2f", legend_axis_labels), las = 2)
      graphics::text(0.5, 0.95, "Overrep.", pos = 4, cex = 0.8)
      graphics::text(0.5, 0.05, "Underrep.", pos = 4, cex = 0.8)
    } # end .create_cs_heatmap

    cat("\n=== CREATING VISUALIZATIONS ===\n")
    if (parameters_plot$detailed && is.list(summary_plot_data)) {
      for (name in names(summary_plot_data)) {
        len_num_plot <- gsub("length_", "", name) # Renamed
        pats_df_plot <- summary_plot_data[[name]] # Renamed
        if (is.data.frame(pats_df_plot) && nrow(pats_df_plot) > 0) {
          tryCatch({
            .create_cs_heatmap(pats_df_plot, groups_plot_names, parameters_plot$top_n, parameters_plot$statistical)
            cat(sprintf("✓ Created %s-length subsequences plot\n", len_num_plot))
          }, error = function(e) {
            cat(sprintf("⚠ Could not create %s-length plot: %s\n", len_num_plot, e$message))
          })
        }
      }
    } else if (is.data.frame(summary_plot_data) && nrow(summary_plot_data) > 0) {
      tryCatch({
        .create_cs_heatmap(summary_plot_data, groups_plot_names, parameters_plot$top_n, parameters_plot$statistical)
        cat("✓ Created combined discriminating patterns plot\n")
      }, error = function(e) {
        cat("⚠ Could not create combined plot:", e$message, "\n")
      })
    }
    cat("Visualization complete!\n\n")
  } # end .create_cs_plots

  if (grDevices::dev.cur() != 1 || interactive()) {
    .create_cs_plots(summary_output_data,
                    list(statistical = statistical, detailed = detailed, top_n = top_n),
                    actual_group_levels)
  } else {
    cat("Skipping plot generation as no active graphics device or session is not interactive.\n")
  }

  # --- Construct Result Object ---
  result_object <- list( # Renamed
    results = analysis_results_list,
    summary = summary_output_data,
    parameters = list(
      groups = actual_group_levels,
      min_length = min_length, max_length = max_length, top_n = top_n,
      detailed = detailed, statistical = statistical, correction = correction,
      test_method = test_method, min_expected = min_expected
    ),
    stats = summary_stats_list
  )
  class(result_object) <- "compare_sequences"

  print(result_object)
  return(result_object)
}


#' Print Method for compare_sequences Objects
#'
#' @param x A compare_sequences object
#' @param ... Additional arguments (unused)
#' @export
print.compare_sequences <- function(x, ...) {
  cat("Sequence Comparison Analysis\n")
  cat("============================\n\n")

  cat("Parameters:\n")
  if (!is.null(x$parameters$groups) && length(x$parameters$groups) > 0) {
    cat("  Compared Groups:", paste(x$parameters$groups, collapse = " vs "), "\n")
  }
  cat("  Subsequence lengths:", x$parameters$min_length, "to", x$parameters$max_length, "\n")
  cat("  Analysis type:", if (x$parameters$statistical) "Statistical testing" else "Discrimination analysis", "\n")
  cat("  Results format:", if (x$parameters$detailed) "Detailed by length" else "Combined", "\n")
  cat("  Top patterns shown:", x$parameters$top_n, "\n\n")

  cat("Summary Statistics:\n")
  cat("  Total patterns analyzed:", x$stats$total_patterns, "\n")
  if (x$parameters$statistical) {
    cat("  Total significant patterns:", x$stats$total_significant,
        sprintf(" (Rate: %.1f%%)", if(x$stats$total_patterns > 0) (x$stats$significance_rate * 100) else 0), "\n")
  } else {
    if (!is.null(x$stats$group_sizes)) {
        for(group_name in names(x$stats$group_sizes)){
            cat("  Sequences in ", group_name, ": ", x$stats$group_sizes[[group_name]], "\n", sep="")
        }
    }
  }
  cat("\n")

  cat("Results Preview (Top 5 or fewer):\n")
  if (x$parameters$detailed && is.list(x$summary)) {
    for (name in names(x$summary)) {
      len_num_p <- gsub("length_", "", name) # Renamed
      pats_df_p <- x$summary[[name]] # Renamed
      if (is.data.frame(pats_df_p) && nrow(pats_df_p) > 0 && "pattern" %in% names(pats_df_p)) {
        cat(sprintf("  %s-length (top pattern): %s\n", len_num_p, pats_df_p$pattern[1]))
      } else {
        cat(sprintf("  %s-length: No patterns to display.\n", len_num_p))
      }
    }
  } else if (is.data.frame(x$summary) && nrow(x$summary) > 0) {
    num_to_show <- min(5, nrow(x$summary))
    for (i in seq_len(num_to_show)) {
      pinfo <- x$summary[i, ] # Renamed
      main_id_col <- "pattern" # Assuming 'pattern' is always the identifier

      if (main_id_col %in% names(pinfo)) {
          id_val <- pinfo[[main_id_col]]
          details <- ""
          if (x$parameters$statistical && "p_value" %in% names(pinfo) && "test_used" %in% names(pinfo)) {
            details <- sprintf("(p=%.3f, Test: %s)", pinfo$p_value, pinfo$test_used)
          } else if (!x$parameters$statistical && "discrimination" %in% names(pinfo)) {
            details <- sprintf("(Discrimination=%.3f)", pinfo$discrimination)
          }
          cat(sprintf("  %d. %s %s\n", i, id_val, details))
      }
    }
  } else {
    cat("  No patterns found to display in preview.\n")
  }

  cat("\nUse summary(object) for detailed results tables.\n")
  cat("Use plot(object) to regenerate visualizations if applicable.\n")
}

#' Summary Method for compare_sequences Objects
#'
#' @param object A compare_sequences object
#' @param ... Additional arguments (unused)
#' @export
summary.compare_sequences <- function(object, ...) {
  # First, print the standard summary information via print()
  # This avoids duplicating parameter and high-level stats output.
  print(object) # This will call print.compare_sequences

  cat("\n--- Detailed Pattern Tables ---\n")
  if (object$parameters$detailed && is.list(object$summary)) {
    if (length(object$summary) == 0) {
        cat("No detailed summary tables available.\n")
    } else {
        for (name in names(object$summary)) {
          len_num_s <- gsub("length_", "", name) # Renamed
          cat(sprintf("\nSubsequences of length %s:\n", len_num_s))
          # cat(paste(rep("-", 25 + nchar(len_num_s)), collapse = ""), "\n") # Dynamic underline

          pats_df_s <- object$summary[[name]] # Renamed
          if (is.data.frame(pats_df_s) && nrow(pats_df_s) > 0) {
            print(pats_df_s, row.names = FALSE)
          } else {
            cat("  No patterns found for this length.\n")
          }
        }
    }
  } else if (is.data.frame(object$summary) && nrow(object$summary) > 0) {
    cat("\nTop patterns (combined across lengths):\n")
    # cat(paste(rep("-", 40), collapse = ""), "\n")
    print(object$summary, row.names = FALSE)
  } else {
    cat("\nNo summary pattern data to display.\n")
  }
}

#' Plot Method for compare_sequences Objects
#'
#' Regenerates visualizations for a `compare_sequences` object.
#'
#' @param x A `compare_sequences` object.
#' @param ... Additional arguments (unused).
#' @export
plot.compare_sequences <- function(x, ...) {
  cat("Attempting to regenerate visualizations for the two-group comparison...\n")

  if (grDevices::dev.cur() == 1 && !interactive()) {
    cat("No active graphics device or session is not interactive. Plotting skipped.\n")
    return(invisible())
  }

  # Define the heatmap function locally within plot, copying its logic
  # from the main compare_sequences function, ensuring it uses `x` (the object)
  .plot_cs_heatmap <- function(patterns_heat_df, groups_heat_names, top_n_heat, is_statistical_heat) {
      if (nrow(patterns_heat_df) == 0) return()
      plot_data_df <- utils::head(patterns_heat_df, top_n_heat)

      col_freq_g1 <- paste0("freq_", groups_heat_names[1])
      col_freq_g2 <- paste0("freq_", groups_heat_names[2])

      residual_matrix <- NULL

      if (is_statistical_heat) {
        if (!(col_freq_g1 %in% names(plot_data_df)) || !(col_freq_g2 %in% names(plot_data_df))) {
          warning("Frequency columns not found for heatmap. Skipping.", call. = FALSE); return()
        }
        total_g1 <- sum(plot_data_df[[col_freq_g1]], na.rm = TRUE)
        total_g2 <- sum(plot_data_df[[col_freq_g2]], na.rm = TRUE)

        if (total_g1 == 0 || total_g2 == 0) {
          if ("prop_diff" %in% names(plot_data_df)) {
            prop_d <- plot_data_df$prop_diff; residual_matrix <- cbind(prop_d, -prop_d)
          } else return()
        } else {
          prop_g1_vals <- plot_data_df[[col_freq_g1]] / total_g1
          prop_g2_vals <- plot_data_df[[col_freq_g2]] / total_g2
          overall_props <- (plot_data_df[[col_freq_g1]] + plot_data_df[[col_freq_g2]]) / (total_g1 + total_g2)
          se_g1 <- sqrt(overall_props * (1 - overall_props) / total_g1)
          se_g2 <- sqrt(overall_props * (1 - overall_props) / total_g2)
          resid_g1 <- ifelse(se_g1 == 0, 0, (prop_g1_vals - overall_props) / se_g1)
          resid_g2 <- ifelse(se_g2 == 0, 0, (prop_g2_vals - overall_props) / se_g2)
          residual_matrix <- cbind(resid_g1, resid_g2)
        }
      } else {
        if (!("prop_diff" %in% names(plot_data_df))) {
          warning("prop_diff not found for discrimination heatmap. Skipping.", call. = FALSE); return()
        }
        prop_d <- plot_data_df$prop_diff; residual_matrix <- cbind(prop_d, -prop_d)
      }

      if(is.null(residual_matrix) || nrow(residual_matrix) == 0) return()

      rownames(residual_matrix) <- plot_data_df$pattern
      colnames(residual_matrix) <- groups_heat_names

      original_layout <- grDevices::layout(matrix(c(1, 2), nrow = 1), widths = c(4, 1))
      on.exit(grDevices::layout(original_layout), add = TRUE)
      graphics::par(mar = c(5, 12, 4, 1))

      max_abs_val <- max(abs(residual_matrix), na.rm = TRUE)
      if (!is.finite(max_abs_val) || max_abs_val == 0) max_abs_val <- 1
      color_palette <- grDevices::colorRampPalette(c("blue", "white", "red"))(100)

      graphics::image(1:ncol(residual_matrix), 1:nrow(residual_matrix), t(residual_matrix),
                      col = color_palette, breaks = seq(-max_abs_val, max_abs_val, length.out = 101),
                      xlab = "Groups", ylab = "", axes = FALSE)
      graphics::axis(1, at = 1:ncol(residual_matrix), labels = colnames(residual_matrix))
      graphics::axis(2, at = 1:nrow(residual_matrix), labels = rownames(residual_matrix), las = 2)

      graphics::par(mar = c(5, 1, 4, 3))
      legend_y_coords <- seq(0, 1, length.out = 100)
      legend_color_matrix <- matrix(legend_y_coords, ncol = 1)
      graphics::image(1, legend_y_coords, t(legend_color_matrix), col = color_palette, axes = FALSE, xlab = "", ylab = "")
      legend_axis_labels <- seq(-max_abs_val, max_abs_val, length.out = 5)
      graphics::axis(4, at = seq(0, 1, length.out = 5), labels = sprintf("%.2f", legend_axis_labels), las = 2)
      graphics::text(0.5, 0.95, "Overrep.", pos = 4, cex = 0.8)
      graphics::text(0.5, 0.05, "Underrep.", pos = 4, cex = 0.8)
  } # end .plot_cs_heatmap


  params <- x$parameters
  summary_data <- x$summary
  group_names_plot <- params$groups

  if (params$detailed && is.list(summary_data)) {
    cat("Generating detailed plots by subsequence length:\n")
    for (name in names(summary_data)) {
      len_num <- gsub("length_", "", name)
      pats_data <- summary_data[[name]]
      if (is.data.frame(pats_data) && nrow(pats_data) > 0) {
        cat(sprintf("  Plotting for %s-length subsequences...\n", len_num))
        tryCatch({
          .plot_cs_heatmap(pats_data, group_names_plot, params$top_n, params$statistical)
        }, error = function(e) {
          cat(sprintf("  ⚠ Could not create %s-length plot: %s\n", len_num, e$message))
        })
      }
    }
  } else if (is.data.frame(summary_data) && nrow(summary_data) > 0) {
    cat("Generating combined plot for top patterns overall:\n")
    tryCatch({
      .plot_cs_heatmap(summary_data, group_names_plot, params$top_n, params$statistical)
    }, error = function(e) {
      cat("  ⚠ Could not create combined plot:", e$message, "\n")
    })
  } else {
    cat("No summary data available to generate plots.\n")
  }
  cat("Plot recreation complete!\n")
}

# ==============================================================================
# MULTI-GROUP SEQUENCE COMPARISON FUNCTION
# ==============================================================================

#' Compare Sequences Across Multiple Groups
#'
#' Performs sequence comparison analysis for every unique pair of groups when
#' more than two groups are present in the data. It utilizes the
#' `compare_sequences` function for each pairwise comparison.
#'
#' @param data A data frame containing sequence data in wide format. Each row
#'   represents a sequence, and columns represent time points or states.
#' @param group A character string specifying the name of the column in `data` that
#'   contains group assignments, or a vector/factor of group assignments
#'   corresponding to each row in `data`.
#' @param min_length Minimum subsequence length to analyze (default: 2).
#' @param max_length Maximum subsequence length to analyze (default: 5).
#' @param top_n Number of top patterns to return from each pairwise comparison's
#'   summary (default: 10).
#' @param detailed Logical. If TRUE, each pairwise comparison returns detailed
#'   results by subsequence length. If FALSE, returns a combined summary.
#'   (default: FALSE).
#' @param statistical Logical. If TRUE, performs statistical testing for each
#'   pairwise comparison (default: FALSE).
#' @param correction Character string. Multiple comparison correction method to apply
#'   within each pairwise statistical test (default: "bonferroni").
#' @param test_method Character string. Statistical test to use ("auto", "fisher",
#'   "chi.squared") within each pairwise comparison (default: "auto").
#' @param min_expected Numeric. Minimum expected count for using Chi-squared test
#'   when `test_method = "auto"` (default: 5).
#'
#' @return An object of class `compare_sequences_multiple_result`. This object is a
#'   list where each element is the result of a `compare_sequences` call for a
#'   unique pair of groups. The list elements are named in the format "GroupA_vs_GroupB".
#'
#' @family Sequence Comparison Functions
#' @seealso \code{\link{compare_sequences}} for the underlying pairwise comparison logic.
#' @export
#' @examples
#' \dontrun{
#'   # Create dummy data with 3 groups
#'   seq_data_multi <- data.frame(
#'     T1 = sample(letters[1:3], 30, replace = TRUE),
#'     T2 = sample(letters[1:3], 30, replace = TRUE),
#'     T3 = sample(letters[1:3], 30, replace = TRUE),
#'     Group = sample(c("X", "Y", "Z"), 30, replace = TRUE)
#'   )
#'
#'   # Perform multi-group comparison
#'   multi_results <- compare_sequences_multiple(seq_data_multi, group = "Group")
#'   print(multi_results)
#'   summary(multi_results)
#'
#'   # To access results for a specific pair, e.g., X vs Y:
#'   # results_X_vs_Y <- multi_results[["X_vs_Y"]]
#' }
compare_sequences_multiple <- function(data, group, min_length = 2, max_length = 5, top_n = 10,
                                     detailed = FALSE, statistical = FALSE, correction = "bonferroni",
                                     test_method = "auto", min_expected = 5) {

  # --- Input Validation ---
  if (!is.data.frame(data)) stop("'data' must be a data frame.", call. = FALSE)
  if (nrow(data) == 0) stop("'data' cannot be empty.", call. = FALSE)

  # Extract sequence data and group vector
  group_col_name <- NULL
  sequence_df <- NULL

  if (is.character(group) && length(group) == 1 && group %in% names(data)) {
    group_col_name <- group
    group_vector <- data[[group_col_name]]
    sequence_df <- data # compare_sequences will handle group column removal internally if name is passed
  } else if ((is.vector(group) || is.factor(group)) && length(group) == nrow(data)) {
    group_vector <- group
    sequence_df <- data # Assumes data is just sequences
    # We need a temporary group column name if group is a vector, for compare_sequences
    # Or, pass the vector directly if compare_sequences is further refactored.
    # For now, let's add a temporary group column for consistency with compare_sequences's current main path.
    temp_group_col_name <- ".temp_group_col_for_cs_multiple"
    # Avoid overwriting if user has such a column
    while(temp_group_col_name %in% names(sequence_df)){ temp_group_col_name <- paste0(temp_group_col_name, "_") }
    sequence_df[[temp_group_col_name]] <- group_vector
    group_col_name <- temp_group_col_name # Use this for calls to compare_sequences

  } else {
    stop("'group' must be a valid column name in 'data', or a vector/factor of the same length as nrow(data).", call. = FALSE)
  }

  # Validate group vector (min_groups = 2, as we need at least two for any comparison)
  # `validate_group_input` returns a factor
  validated_group_factor <- validate_group_input(group_vector, nrow(sequence_df), min_groups = 2)
  actual_group_levels <- levels(validated_group_factor)
  num_groups <- length(actual_group_levels)

  if (num_groups < 2) {
    stop("At least two groups are required for comparison.", call. = FALSE) # Should be caught by validate_group_input
  }

  # If exactly two groups, just call compare_sequences directly and wrap its class for consistency, or return directly.
  # For now, let's make it consistent: always return a list, even for 2 groups, named by the pair.
  # This simplifies downstream processing of the compare_sequences_multiple_result object.

  all_pairwise_results <- list()

  if (num_groups == 2) {
    cat("Exactly two groups found. Performing a single pairwise comparison...\n")
    # Call compare_sequences with the original data and group specification
    # If group was a vector, sequence_df now has the temp group column, and group_col_name points to it.
    # If group was a name, sequence_df is original data, group_col_name is original name.

    # We need to ensure 'group' argument for compare_sequences is the name of the group column
    # within the 'data' it receives.
    current_pair_data <- sequence_df # sequence_df has the group column (original or temporary)

    pair_result <- compare_sequences(
      data = current_pair_data,
      group = group_col_name, # This must be the name of the group column in current_pair_data
      min_length = min_length, max_length = max_length, top_n = top_n,
      detailed = detailed, statistical = statistical, correction = correction,
      test_method = test_method, min_expected = min_expected
    )
    pair_name <- paste(actual_group_levels[1], "vs", actual_group_levels[2], sep = "_")
    all_pairwise_results[[pair_name]] <- pair_result
  } else {
    cat(sprintf("%d groups found. Performing pairwise comparisons for all unique pairs...\n", num_groups))
    group_pairs <- generate_group_pairs(actual_group_levels)

    for (pair in group_pairs) {
      group1_name <- pair[1]
      group2_name <- pair[2]
      pair_name <- paste(group1_name, "vs", group2_name, sep = "_")
      cat(sprintf("\n--- Comparing: %s ---\n", pair_name))

      # Subset data for the current pair
      # The 'sequence_df' already has the (potentially temporary) group column.
      # We use 'validated_group_factor' for subsetting, which is based on original group_vector.
      pair_indices <- validated_group_factor %in% c(group1_name, group2_name)
      current_pair_data <- sequence_df[pair_indices, , drop = FALSE]

      # The 'group_col_name' is the name of the column in current_pair_data to use for grouping.
      pair_result <- tryCatch({
        compare_sequences(
          data = current_pair_data,
          group = group_col_name,
          min_length = min_length, max_length = max_length, top_n = top_n,
          detailed = detailed, statistical = statistical, correction = correction,
          test_method = test_method, min_expected = min_expected
        )
      }, error = function(e) {
        warning(sprintf("Error during comparison of %s: %s. Skipping this pair.", pair_name, e$message), call. = FALSE)
        return(NULL) # Store NULL if a pair fails
      })

      all_pairwise_results[[pair_name]] <- pair_result
    }
  }

  # Clean up temporary group column if it was added
  if (!is.null(temp_group_col_name) && temp_group_col_name %in% names(sequence_df)) {
      # This cleanup is tricky if sequence_df was passed into compare_sequences and modified there.
      # However, compare_sequences takes 'data' and 'group' (name), and subsets data internally.
      # The 'data' for compare_sequences_multiple should remain untouched.
      # The temporary column was added to `sequence_df` which was a copy of input `data` or `data` itself.
      # This is fine as `sequence_df` is local to this function scope.
  }


  class(all_pairwise_results) <- "compare_sequences_multiple_result"
  cat("\nMulti-group comparison complete.\n")
  return(all_pairwise_results)
}

#' Print Method for compare_sequences_multiple_result Objects
#'
#' Provides a concise overview of the multi-group sequence comparison results.
#' It lists each pairwise comparison performed and a high-level summary for each.
#'
#' @param x An object of class `compare_sequences_multiple_result`.
#' @param ... Additional arguments (unused).
#' @return Invisibly returns the input object `x`.
#' @export
print.compare_sequences_multiple_result <- function(x, ...) {
  cat("Multi-Group Sequence Comparison Results\n")
  cat("=======================================\n")
  num_comparisons <- length(x)
  cat("Number of pairwise comparisons performed:", num_comparisons, "\n\n")

  if (num_comparisons == 0) {
    cat("No comparisons were performed or results available.\n")
    return(invisible(x))
  }

  for (i in seq_along(x)) {
    pair_name <- names(x)[i]
    result_item <- x[[i]]

    cat(paste0("--- Comparison: ", pair_name, " ---\n"))
    if (is.null(result_item)) {
        cat("  Result for this pair is NULL (an error likely occurred).\n\n")
    } else if (inherits(result_item, "compare_sequences")) {
      # Briefly summarize the result from the compare_sequences object
      cat("  Analysis type:", if (result_item$parameters$statistical) "Statistical" else "Discrimination", "\n")
      cat("  Total patterns analyzed in this pair:", result_item$stats$total_patterns, "\n")
      if (result_item$parameters$statistical) {
        cat("  Significant patterns in this pair:", result_item$stats$total_significant, "\n")
      }
      # Show top 1 pattern from summary to give a glimpse
      if(is.data.frame(result_item$summary) && nrow(result_item$summary) > 0 && "pattern" %in% names(result_item$summary)){
          cat("  Top overall pattern for this pair: ", result_item$summary$pattern[1], "\n")
      } else if (is.list(result_item$summary) && length(result_item$summary) > 0 &&
                 is.data.frame(result_item$summary[[1]]) && nrow(result_item$summary[[1]]) > 0 &&
                 "pattern" %in% names(result_item$summary[[1]])) {
          # Case for detailed results
          cat("  Top overall pattern for this pair (first length category): ", result_item$summary[[1]]$pattern[1], "\n")
      } else {
          cat("  No top patterns to display for this pair.\n")
      }
      cat("\n")
    } else {
      cat("  Result for this pair is not a 'compare_sequences' object.\n\n")
    }
  }
  cat("Use summary(object) for more details or access individual pairwise results using [[name]], e.g., results[['GroupA_vs_GroupB']].\n")
  invisible(x)
}

#' Summary Method for compare_sequences_multiple_result Objects
#'
#' Provides a detailed summary of multi-group sequence comparison results.
#' It first gives an overall count of patterns and significant patterns (if applicable)
#' across all pairwise comparisons, then optionally prints the full summary
#' for each individual pairwise comparison.
#'
#' @param object An object of class `compare_sequences_multiple_result`.
#' @param ... Additional arguments, typically passed to `summary.compare_sequences`
#'   for each pairwise summary (e.g., `top_n` for those individual summaries, though
#'   `top_n` for `compare_sequences_multiple` already determined this for storage).
#' @return Invisibly returns the input object `object`.
#' @export
summary.compare_sequences_multiple_result <- function(object, ...) {
  cat("Summary of Multi-Group Sequence Comparison\n")
  cat("========================================\n")
  num_comparisons <- length(object)
  cat("Total pairwise comparisons performed:", num_comparisons, "\n\n")

  if (num_comparisons == 0) {
    cat("No comparisons to summarize.\n")
    return(invisible(object))
  }

  # Consolidate some key stats across all comparisons
  all_total_patterns <- 0
  all_significant_patterns <- 0 # Only if statistical
  num_statistical_comparisons <- 0

  for (pair_name in names(object)) {
    result_item <- object[[pair_name]]
    if (!is.null(result_item) && inherits(result_item, "compare_sequences")) {
      all_total_patterns <- all_total_patterns + result_item$stats$total_patterns
      if (result_item$parameters$statistical) {
        all_significant_patterns <- all_significant_patterns + result_item$stats$total_significant
        num_statistical_comparisons <- num_statistical_comparisons + 1
      }
    }
  }

  cat("Overall Summary Across All Pairs:\n")
  cat("  Total patterns analyzed across all comparisons:", all_total_patterns, "\n")
  if (num_statistical_comparisons > 0) {
    cat("  Total significant patterns across all statistical comparisons:", all_significant_patterns, "\n")
  }
  cat("\nFor detailed summaries of each pairwise comparison, print the individual result objects:\n")

  for (pair_name in names(object)) {
      cat(paste0("\n--- Summary for: ", pair_name, " ---\n"))
      result_item <- object[[pair_name]]
      if (!is.null(result_item) && inherits(result_item, "compare_sequences")) {
          summary(result_item) # Calls summary.compare_sequences
      } else {
          cat("  No summary available for this pair (result was NULL or invalid).\n")
      }
  }
  invisible(object)
}
#'
#' Performs comprehensive sequence comparison analysis between two groups, identifying
#' discriminating subsequences and their statistical associations. Supports both
#' discrimination-based analysis and rigorous statistical testing.
#'
#' @param data A data frame containing sequence data in wide format, where each row
#'   represents a sequence and each column represents a time point
#' @param group Column name or index containing group information, or a vector indicating 
#'   group membership for each sequence. If data contains a group column, specify column 
#'   name (e.g., "Group") or pass as vector (e.g., data$Group)
#' @param min_length Minimum subsequence length to analyze (default: 2)
#' @param max_length Maximum subsequence length to analyze (default: 5)
#' @param top_n Number of top patterns to return and display (default: 10)
#' @param detailed Logical, whether to return detailed results by subsequence length
#'   (default: FALSE). When TRUE, results are organized by length; when FALSE,
#'   shows overall top patterns across all lengths
#' @param statistical Logical, whether to perform statistical testing (default: FALSE).
#'   When TRUE, performs Chi-squared or Fisher's Exact tests; when FALSE, uses
#'   discrimination score analysis
#' @param correction Character, multiple comparison correction method (default: "bonferroni").
#'   Supports: "bonferroni", "holm", "hochberg", "BH", "BY", "none"
#' @param test_method Character, statistical test selection (default: "auto").
#'   Options: "auto", "fisher", "chi.squared"
#' @param min_expected Numeric, minimum expected count for automatic test selection (default: 5)
#'
#' @return A compare_sequences object containing:
#' \describe{
#'   \item{results}{Main analysis results}
#'   \item{summary}{Summary table of top patterns}
#'   \item{parameters}{Analysis parameters used}
#'   \item{stats}{Summary statistics}
#' }
#'
#' @examples
#' \dontrun{
#' # Load example data
#' data(seqdata)
#' 
#' # Basic analysis using column name
#' result <- compare_sequences(seqdata, "Group")
#' print(result)
#' 
#' # Basic analysis using vector
#' result <- compare_sequences(seqdata, seqdata$Group)
#' print(result)
#' 
#' # Statistical analysis
#' result <- compare_sequences(seqdata, "Group", statistical = TRUE)
#' }
#'
#' @export
compare_sequences <- function(data, group, min_length = 2, max_length = 5, top_n = 10, 
                             detailed = FALSE, statistical = FALSE, correction = "bonferroni", 
                             test_method = "auto", min_expected = 5) {
  
  # =====================================================================
  # INPUT VALIDATION
  # =====================================================================
  
  # Validate inputs
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame", call. = FALSE)
  }
  
  if (nrow(data) == 0) {
    stop("'data' cannot be empty", call. = FALSE)
  }
  
  if (!is.numeric(min_length) || !is.numeric(max_length) || min_length < 1 || max_length < min_length) {
    stop("'min_length' and 'max_length' must be positive integers with min_length <= max_length", call. = FALSE)
  }
  
  if (!is.numeric(top_n) || top_n < 1) {
    stop("'top_n' must be a positive integer", call. = FALSE)
  }
  
  if (!is.logical(detailed) || !is.logical(statistical)) {
    stop("'detailed' and 'statistical' must be logical values", call. = FALSE)
  }
  
  if (!correction %in% c("bonferroni", "holm", "hochberg", "BH", "BY", "none")) {
    stop("'correction' must be one of: bonferroni, holm, hochberg, BH, BY, none", call. = FALSE)
  }
  
  if (!test_method %in% c("auto", "fisher", "chi.squared")) {
    stop("'test_method' must be one of: auto, fisher, chi.squared", call. = FALSE)
  }
  
  # =====================================================================
  # DATA PREPROCESSING
  # =====================================================================
  
  # Handle group parameter - can be column name or vector
  if (is.character(group) && length(group) == 1) {
    # group is a column name
    if (!group %in% names(data)) {
      stop("Group column '", group, "' not found in data", call. = FALSE)
    }
    group_vector <- data[[group]]
    # Remove group column from data
    data <- data[, !names(data) %in% group, drop = FALSE]
  } else if (is.numeric(group) && length(group) == 1) {
    # group is a column index
    if (group < 1 || group > ncol(data)) {
      stop("Group column index ", group, " is out of range", call. = FALSE)
    }
    group_vector <- data[[group]]
    # Remove group column from data
    data <- data[, -group, drop = FALSE]
  } else {
    # group is a vector
    if (length(group) != nrow(data)) {
      stop("'group' must have the same length as number of rows in 'data'", call. = FALSE)
    }
    group_vector <- group
  }
  
  # Convert group to factor and validate
  group_factor <- as.factor(group_vector)
  groups <- levels(group_factor)
  
  if (length(groups) != 2) {
    stop("Group variable must have exactly 2 levels, found: ", length(groups), call. = FALSE)
  }
  
  # Split data by groups
  seq_A <- data[group_factor == groups[1], , drop = FALSE]
  seq_B <- data[group_factor == groups[2], , drop = FALSE]
  
  # Clean data: remove empty rows and convert empty strings to NA
  seq_A <- seq_A[apply(seq_A, 1, function(x) !all(is.na(x) | x == "")), , drop = FALSE]
  seq_B <- seq_B[apply(seq_B, 1, function(x) !all(is.na(x) | x == "")), , drop = FALSE]
  
  if (nrow(seq_A) == 0 || nrow(seq_B) == 0) {
    stop("Both groups must have at least one non-empty sequence", call. = FALSE)
  }
  
  # Convert empty strings to NA
  seq_A[seq_A == ""] <- NA
  seq_B[seq_B == ""] <- NA
  
  # =====================================================================
  # CORE ANALYSIS FUNCTIONS
  # =====================================================================
  
  extract_subsequences <- function(sequences, n) {
    if (n < 1) return(character(0))
    
    subsequences <- character(0)
    for (i in seq_len(nrow(sequences))) {
      seq_row <- unlist(sequences[i, ])
      # Filter out NA, empty strings, and single characters (potential group markers)
      clean_seq <- seq_row[!is.na(seq_row) & nchar(as.character(seq_row)) > 1]
      
      if (length(clean_seq) >= n) {
        for (j in seq_len(length(clean_seq) - n + 1)) {
          subseq <- paste(clean_seq[j:(j + n - 1)], collapse = "-")
          subsequences <- c(subsequences, subseq)
        }
      }
    }
    return(subsequences)
  }
  
  analyze_discrimination <- function(seq_A, seq_B, n) {
    # Extract subsequences
    subseq_A <- extract_subsequences(seq_A, n)
    subseq_B <- extract_subsequences(seq_B, n)
    
    if (length(subseq_A) == 0 && length(subseq_B) == 0) {
      return(list(n = n, patterns = data.frame(), total_A = 0, total_B = 0, n_patterns = 0))
    }
    
    # Calculate frequencies
    freq_A <- table(subseq_A)
    freq_B <- table(subseq_B)
    all_patterns <- unique(c(names(freq_A), names(freq_B)))
    
    # Create frequency table
    patterns <- data.frame(
      pattern = all_patterns,
      freq_A = as.numeric(freq_A[all_patterns]),
      freq_B = as.numeric(freq_B[all_patterns]),
      stringsAsFactors = FALSE
    )
    patterns[is.na(patterns)] <- 0
    
    # Calculate metrics
    total_A <- sum(patterns$freq_A)
    total_B <- sum(patterns$freq_B)
    
    if (total_A > 0 && total_B > 0) {
      patterns$prop_A <- patterns$freq_A / total_A
      patterns$prop_B <- patterns$freq_B / total_B
      patterns$prop_diff <- patterns$prop_A - patterns$prop_B
      
      # Discrimination score (simple and robust)
      patterns$discrimination <- abs(patterns$prop_diff) * sqrt(patterns$freq_A + patterns$freq_B)
      
      # Sort by discrimination
      patterns <- patterns[order(patterns$discrimination, decreasing = TRUE), ]
    } else {
      patterns$prop_A <- patterns$prop_B <- patterns$prop_diff <- patterns$discrimination <- 0
    }
    
    return(list(
      n = n,
      patterns = patterns,
      total_A = total_A,
      total_B = total_B,
      n_patterns = nrow(patterns)
    ))
  }
  
  analyze_statistical <- function(seq_A, seq_B, n, correction_method, test_method, min_expected) {
    # Extract subsequences
    subseq_A <- extract_subsequences(seq_A, n)
    subseq_B <- extract_subsequences(seq_B, n)
    
    if (length(subseq_A) == 0 && length(subseq_B) == 0) {
      return(list(n = n, patterns = data.frame(), n_patterns = 0, n_significant = 0))
    }
    
    # Calculate frequencies
    freq_A <- table(subseq_A)
    freq_B <- table(subseq_B)
    all_patterns <- unique(c(names(freq_A), names(freq_B)))
    
    # Create results table
    results <- data.frame(
      pattern = all_patterns,
      freq_A = as.numeric(freq_A[all_patterns]),
      freq_B = as.numeric(freq_B[all_patterns]),
      stringsAsFactors = FALSE
    )
    results[is.na(results)] <- 0
    
    # Perform statistical tests
    results$p_value <- NA
    results$test_used <- NA
    results$statistic <- NA
    
    for (i in seq_len(nrow(results))) {
      # Create contingency table
      cont_table <- matrix(c(results$freq_A[i], results$freq_B[i],
                            sum(results$freq_A) - results$freq_A[i],
                            sum(results$freq_B) - results$freq_B[i]), 
                          nrow = 2, byrow = TRUE)
      
      # Determine test method
      use_fisher <- FALSE
      if (test_method == "fisher") {
        use_fisher <- TRUE
      } else if (test_method == "auto") {
        expected_min <- min(chisq.test(cont_table)$expected)
        use_fisher <- expected_min < min_expected
      }
      
      # Perform test
      tryCatch({
        if (use_fisher) {
          test_result <- fisher.test(cont_table)
          results$p_value[i] <- test_result$p.value
          results$test_used[i] <- "Fisher"
          results$statistic[i] <- NA
        } else {
          test_result <- chisq.test(cont_table)
          results$p_value[i] <- test_result$p.value
          results$test_used[i] <- "Chi-squared"
          results$statistic[i] <- test_result$statistic
        }
      }, error = function(e) {
        results$p_value[i] <<- NA
        results$test_used[i] <<- "Failed"
        results$statistic[i] <<- NA
      })
    }
    
    # Multiple comparison correction
    if (correction_method != "none") {
      results$p_adjusted <- p.adjust(results$p_value, method = correction_method)
    } else {
      results$p_adjusted <- results$p_value
    }
    
    results$significant <- !is.na(results$p_adjusted) & results$p_adjusted < 0.05
    
    # Sort by significance and p-value
    results <- results[order(results$significant, results$p_value, decreasing = c(TRUE, FALSE)), ]
    
    return(list(
      n = n,
      patterns = results,
      n_patterns = nrow(results),
      n_significant = sum(results$significant, na.rm = TRUE)
    ))
  }
  
  # =====================================================================
  # MAIN ANALYSIS
  # =====================================================================
  
  cat("Analyzing sequences between groups", groups[1], "and", groups[2], "\n")
  cat("Group sizes:", nrow(seq_A), "vs", nrow(seq_B), "sequences\n")
  cat("Subsequence lengths:", min_length, "to", max_length, "\n\n")
  
  # Run analysis for each subsequence length
  analysis_results <- list()
  for (n in min_length:max_length) {
    if (statistical) {
      analysis_results[[paste0("length_", n)]] <- analyze_statistical(
        seq_A, seq_B, n, correction, test_method, min_expected)
    } else {
      analysis_results[[paste0("length_", n)]] <- analyze_discrimination(seq_A, seq_B, n)
    }
  }
  
  # =====================================================================
  # PREPARE OUTPUT
  # =====================================================================
  
  # Create summary tables
  if (detailed) {
    # Separate table for each length
    summary_tables <- list()
    for (n in min_length:max_length) {
      result <- analysis_results[[paste0("length_", n)]]
      if (nrow(result$patterns) > 0) {
        summary_tables[[paste0("length_", n)]] <- head(result$patterns, top_n)
      }
    }
    summary_output <- summary_tables
  } else {
    # Combined table across all lengths
    all_patterns <- data.frame()
    for (n in min_length:max_length) {
      result <- analysis_results[[paste0("length_", n)]]
      if (nrow(result$patterns) > 0) {
        patterns_with_length <- cbind(length = n, result$patterns)
        all_patterns <- rbind(all_patterns, patterns_with_length)
      }
    }
    
    if (nrow(all_patterns) > 0) {
      if (statistical) {
        all_patterns <- all_patterns[order(all_patterns$significant, all_patterns$p_value, 
                                         decreasing = c(TRUE, FALSE)), ]
      } else {
        all_patterns <- all_patterns[order(all_patterns$discrimination, decreasing = TRUE), ]
      }
      summary_output <- head(all_patterns, top_n)
    } else {
      summary_output <- data.frame()
    }
  }
  
  # Create summary statistics
  total_patterns <- sum(sapply(analysis_results, function(x) x$n_patterns))
  if (statistical) {
    total_significant <- sum(sapply(analysis_results, function(x) x$n_significant))
    summary_stats <- list(
      total_patterns = total_patterns,
      total_significant = total_significant,
      significance_rate = if (total_patterns > 0) total_significant / total_patterns else 0
    )
  } else {
    summary_stats <- list(
      total_patterns = total_patterns,
      group_A_size = nrow(seq_A),
      group_B_size = nrow(seq_B)
    )
  }
  
  # =====================================================================
  # CREATE VISUALIZATIONS
  # =====================================================================
  
  create_plots <- function(analysis_results, summary_output, parameters, groups) {
    # Helper function for creating heatmaps
    create_heatmap <- function(patterns, groups) {
      if (nrow(patterns) == 0) return()
      
      # Create residual matrix for heatmap
      if (parameters$statistical) {
        # For statistical analysis, use adjusted standardized residuals for better visualization
        residual_data <- patterns[1:min(nrow(patterns), parameters$top_n), ]
        
        # Calculate better residuals for visualization
        total_A <- sum(residual_data$freq_A)
        total_B <- sum(residual_data$freq_B)
        total_overall <- total_A + total_B
        
        # Calculate proportions within each group
        prop_A <- residual_data$freq_A / total_A
        prop_B <- residual_data$freq_B / total_B
        
        # Use standardized proportion differences (z-score like)
        # This shows how much each pattern deviates from equal representation
        overall_prop <- (residual_data$freq_A + residual_data$freq_B) / total_overall
        
        # Calculate standard error for proportion difference
        se_A <- sqrt(overall_prop * (1 - overall_prop) / total_A)
        se_B <- sqrt(overall_prop * (1 - overall_prop) / total_B)
        
        # Standardized deviations from expected proportion
        resid_A <- (prop_A - overall_prop) / se_A
        resid_B <- (prop_B - overall_prop) / se_B
        
        residual_matrix <- cbind(resid_A, resid_B)
      } else {
        # For discrimination analysis, use proportion differences
        residual_data <- patterns[1:min(nrow(patterns), parameters$top_n), ]
        prop_diff <- residual_data$prop_diff
        residual_matrix <- cbind(prop_diff, -prop_diff)
      }
      
      rownames(residual_matrix) <- residual_data$pattern
      colnames(residual_matrix) <- groups
      
      # Create the heatmap
      # Set up layout for main plot + legend
      layout(matrix(c(1, 2), nrow = 1), widths = c(4, 1))
      par(mar = c(5, 12, 4, 1))
      
      # Color palette
      max_val <- max(abs(residual_matrix), na.rm = TRUE)
      colors <- colorRampPalette(c("blue", "white", "red"))(100)
      
      # Create heatmap
      image(1:ncol(residual_matrix), 1:nrow(residual_matrix), 
            t(residual_matrix), 
            col = colors,
            breaks = seq(-max_val, max_val, length.out = 101),
            xlab = "Groups", ylab = "",
            axes = FALSE)
      
      # Add axes
      axis(1, at = 1:ncol(residual_matrix), labels = colnames(residual_matrix))
      axis(2, at = 1:nrow(residual_matrix), labels = rownames(residual_matrix), las = 2)
      
      # Add gradient legend
      # Reset margins for legend
      par(mar = c(5, 1, 4, 3))
      
      # Create gradient legend
      legend_y <- seq(0, 1, length.out = 100)
      legend_matrix <- matrix(legend_y, ncol = 1)
      
      image(1, legend_y, t(legend_matrix), 
            col = colors, axes = FALSE, xlab = "", ylab = "")
      
      # Add legend labels
      legend_vals <- seq(-max_val, max_val, length.out = 5)
      axis(4, at = seq(0, 1, length.out = 5), labels = sprintf("%.2f", legend_vals), las = 2)
      
      # Add "Overrep." and "Underrep." labels
      text(0.5, 0.95, "Overrep.", pos = 4, cex = 0.8)
      text(0.5, 0.05, "Underrep.", pos = 4, cex = 0.8)
      
      # Reset layout
      layout(1)
    }
    
    cat("\n=== CREATING VISUALIZATIONS ===\n")
    
    if (parameters$detailed && is.list(summary_output)) {
      # Create separate plots for each length
      for (name in names(summary_output)) {
        length_num <- gsub("length_", "", name)
        patterns <- summary_output[[name]]
        
        if (nrow(patterns) > 0) {
          tryCatch({
            create_heatmap(patterns, groups)
            cat(sprintf("✓ Created %s-length subsequences plot\n", length_num))
          }, error = function(e) {
            cat(sprintf("⚠ Could not create %s-length plot: %s\n", length_num, e$message))
          })
        }
      }
    } else if (is.data.frame(summary_output) && nrow(summary_output) > 0) {
      # Create combined plot
      tryCatch({
        create_heatmap(summary_output, groups)
        cat("✓ Created combined discriminating patterns plot\n")
      }, error = function(e) {
        cat("⚠ Could not create combined plot:", e$message, "\n")
      })
    }
    
    cat("Visualization complete!\n\n")
  }
  
  # Create plots
  create_plots(analysis_results, summary_output, list(
    statistical = statistical,
    detailed = detailed,
    top_n = top_n
  ), groups)
  
  # =====================================================================
  # DISPLAY RESULTS
  # =====================================================================
  
  # Create result object first
  result <- list(
    results = analysis_results,
    summary = summary_output,
    parameters = list(
      groups = groups,
      min_length = min_length,
      max_length = max_length,
      top_n = top_n,
      detailed = detailed,
      statistical = statistical,
      correction = correction,
      test_method = test_method,
      min_expected = min_expected
    ),
    stats = summary_stats
  )
  
  class(result) <- "compare_sequences"
  
  # Print results immediately
  print(result)
  
  return(result)
}

#' Print Method for compare_sequences Objects
#'
#' @param x A compare_sequences object
#' @param ... Additional arguments (unused)
#' @export
print.compare_sequences <- function(x, ...) {
  cat("Sequence Comparison Analysis\n")
  cat("============================\n\n")
  
  # Parameters
  cat("Parameters:\n")
  cat("  Groups:", paste(x$parameters$groups, collapse = " vs "), "\n")
  cat("  Subsequence lengths:", x$parameters$min_length, "to", x$parameters$max_length, "\n")
  cat("  Analysis type:", if (x$parameters$statistical) "Statistical testing" else "Discrimination analysis", "\n")
  cat("  Results format:", if (x$parameters$detailed) "Detailed by length" else "Combined", "\n")
  cat("  Top patterns shown:", x$parameters$top_n, "\n\n")
  
  # Summary statistics
  cat("Summary:\n")
  cat("  Total patterns found:", x$stats$total_patterns, "\n")
  if (x$parameters$statistical) {
    cat("  Significant patterns:", x$stats$total_significant, 
        sprintf("(%.1f%%)", x$stats$significance_rate * 100), "\n")
  } else {
    cat("  Group A sequences:", x$stats$group_A_size, "\n")
    cat("  Group B sequences:", x$stats$group_B_size, "\n")
  }
  cat("\n")
  
  # Results preview
  if (x$parameters$detailed && is.list(x$summary)) {
    cat("Top patterns by subsequence length:\n")
    for (name in names(x$summary)) {
      length_num <- gsub("length_", "", name)
      patterns <- x$summary[[name]]
      if (nrow(patterns) > 0) {
        cat(sprintf("  %s-length: %d patterns (top: %s)\n", 
                   length_num, nrow(patterns), patterns$pattern[1]))
      }
    }
  } else if (is.data.frame(x$summary) && nrow(x$summary) > 0) {
    cat("Top patterns overall:\n")
    for (i in seq_len(min(5, nrow(x$summary)))) {
      pattern_info <- x$summary[i, ]
      if (x$parameters$statistical) {
        cat(sprintf("  %d. %s (p=%.3f, %s)\n", 
                   i, pattern_info$pattern, pattern_info$p_value, pattern_info$test_used))
      } else {
        cat(sprintf("  %d. %s (discrimination=%.3f)\n", 
                   i, pattern_info$pattern, pattern_info$discrimination))
      }
    }
  } else {
    cat("No significant patterns found.\n")
  }
  
  cat("\nUse summary() for detailed results.\n")
}

#' Summary Method for compare_sequences Objects
#'
#' @param object A compare_sequences object
#' @param ... Additional arguments (unused)
#' @export
summary.compare_sequences <- function(object, ...) {
  cat("Detailed Sequence Comparison Results\n")
  cat("====================================\n\n")
  
  # Print parameters
  print(object)
  
  # Detailed results
  if (object$parameters$detailed && is.list(object$summary)) {
    for (name in names(object$summary)) {
      length_num <- gsub("length_", "", name)
      cat(sprintf("\n%s-length subsequences:\n", length_num))
      cat(paste(rep("-", 40), collapse = ""), "\n")
      
      patterns <- object$summary[[name]]
      if (nrow(patterns) > 0) {
        print(patterns, row.names = FALSE)
      } else {
        cat("No patterns found.\n")
      }
    }
  } else if (is.data.frame(object$summary) && nrow(object$summary) > 0) {
    cat("\nTop patterns across all lengths:\n")
    cat(paste(rep("-", 50), collapse = ""), "\n")
    print(object$summary, row.names = FALSE)
  } else {
    cat("\nNo patterns found.\n")
  }
}

#' Plot Method for compare_sequences Objects
#'
#' @param x A compare_sequences object
#' @param ... Additional arguments (unused)
#' @export
plot.compare_sequences <- function(x, ...) {
  cat("Recreating visualizations...\n")
  
  # Recreate the plots using the stored data
  create_plots <- function(summary_output, parameters, groups) {
    # Helper function for creating heatmaps
    create_heatmap <- function(patterns, groups) {
      if (nrow(patterns) == 0) return()
      
      # Create residual matrix for heatmap
      if (parameters$statistical) {
        # For statistical analysis, use adjusted standardized residuals for better visualization
        residual_data <- patterns[1:min(nrow(patterns), parameters$top_n), ]
        
        # Calculate better residuals for visualization
        total_A <- sum(residual_data$freq_A)
        total_B <- sum(residual_data$freq_B)
        total_overall <- total_A + total_B
        
        # Calculate proportions within each group
        prop_A <- residual_data$freq_A / total_A
        prop_B <- residual_data$freq_B / total_B
        
        # Use standardized proportion differences (z-score like)
        # This shows how much each pattern deviates from equal representation
        overall_prop <- (residual_data$freq_A + residual_data$freq_B) / total_overall
        
        # Calculate standard error for proportion difference
        se_A <- sqrt(overall_prop * (1 - overall_prop) / total_A)
        se_B <- sqrt(overall_prop * (1 - overall_prop) / total_B)
        
        # Standardized deviations from expected proportion
        resid_A <- (prop_A - overall_prop) / se_A
        resid_B <- (prop_B - overall_prop) / se_B
        
        residual_matrix <- cbind(resid_A, resid_B)
      } else {
        # For discrimination analysis, use proportion differences
        residual_data <- patterns[1:min(nrow(patterns), parameters$top_n), ]
        prop_diff <- residual_data$prop_diff
        residual_matrix <- cbind(prop_diff, -prop_diff)
      }
      
      rownames(residual_matrix) <- residual_data$pattern
      colnames(residual_matrix) <- groups
      
      # Create the heatmap
      par(mar = c(5, 12, 4, 1))
      # Set up layout for main plot + legend
      layout(matrix(c(1, 2), nrow = 1), widths = c(4, 1))
      
      # Color palette
      max_val <- max(abs(residual_matrix), na.rm = TRUE)
      colors <- colorRampPalette(c("blue", "white", "red"))(100)
      
      # Create heatmap
      image(1:ncol(residual_matrix), 1:nrow(residual_matrix), 
            t(residual_matrix), 
            col = colors,
            breaks = seq(-max_val, max_val, length.out = 101),
            xlab = "Groups", ylab = "",
            axes = FALSE)
      
      # Add axes
      axis(1, at = 1:ncol(residual_matrix), labels = colnames(residual_matrix))
      axis(2, at = 1:nrow(residual_matrix), labels = rownames(residual_matrix), las = 2)
      
      # Add gradient legend
      # Reset margins for legend
      par(mar = c(5, 1, 4, 3))
      
      # Create gradient legend
      legend_y <- seq(0, 1, length.out = 100)
      legend_matrix <- matrix(legend_y, ncol = 1)
      
      image(1, legend_y, t(legend_matrix), 
            col = colors, axes = FALSE, xlab = "", ylab = "")
      
      # Add legend labels
      legend_vals <- seq(-max_val, max_val, length.out = 5)
      axis(4, at = seq(0, 1, length.out = 5), labels = sprintf("%.2f", legend_vals), las = 2)
      
      # Add "Overrep." and "Underrep." labels
      text(0.5, 0.95, "Overrep.", pos = 4, cex = 0.8)
      text(0.5, 0.05, "Underrep.", pos = 4, cex = 0.8)
      
      # Reset layout
      layout(1)
    }
    
    if (parameters$detailed && is.list(summary_output)) {
      # Create separate plots for each length
      for (name in names(summary_output)) {
        length_num <- gsub("length_", "", name)
        patterns <- summary_output[[name]]
        
        if (nrow(patterns) > 0) {
          tryCatch({
            create_heatmap(patterns, groups)
            cat(sprintf("✓ Created %s-length subsequences plot\n", length_num))
          }, error = function(e) {
            cat(sprintf("⚠ Could not create %s-length plot: %s\n", length_num, e$message))
          })
        }
      }
    } else if (is.data.frame(summary_output) && nrow(summary_output) > 0) {
      # Create combined plot
      tryCatch({
        create_heatmap(summary_output, groups)
        cat("✓ Created combined discriminating patterns plot\n")
      }, error = function(e) {
        cat("⚠ Could not create combined plot:", e$message, "\n")
      })
    }
  }
  
  # Create plots
  create_plots(x$summary, x$parameters, x$parameters$groups)
  cat("Plot recreation complete!\n")
} 
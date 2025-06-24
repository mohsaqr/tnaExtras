# ==============================================================================
# SEQUENCE COMPARISON TOOLKIT
# ==============================================================================
# 
# Advanced sequence comparison functions with statistical testing capabilities.
# Includes chi-square tests, Fisher's exact tests, and comprehensive pattern
# analysis with visualization support.
#
# Functions:
# - seqompare(): Main sequence comparison function with statistical testing
# - Helper functions for statistical analysis and visualization
#
# ==============================================================================

# Pattern analysis functions are available within the package

# ==============================================================================
# MAIN SEQUENCE COMPARISON FUNCTION
# ==============================================================================

#' Compare Sequences Between Groups
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
#' @param min_char_length_state Minimum character length for a state to be considered valid during subsequence extraction (default: 1, meaning all non-empty states are considered).
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
                             test_method = "auto", min_expected = 5,
                             min_char_length_state = 1) { # New parameter
  
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
  
  extract_subsequences <- function(sequences, n, min_char_length_state_local) { # Renamed to avoid conflict
    if (n < 1) return(character(0))
    
    subsequences <- character(0)
    for (i in seq_len(nrow(sequences))) {
      seq_row <- unlist(sequences[i, ])
      # Filter out NA, empty strings
      clean_seq_intermediate <- seq_row[!is.na(seq_row) & nzchar(as.character(seq_row))]

      # Filter based on min_char_length_state_local
      if (min_char_length_state_local > 0) {
        clean_seq <- clean_seq_intermediate[nchar(as.character(clean_seq_intermediate)) >= min_char_length_state_local]
      } else {
        clean_seq <- clean_seq_intermediate
      }
      
      if (length(clean_seq) >= n) {
        for (j in seq_len(length(clean_seq) - n + 1)) {
          subseq <- paste(clean_seq[j:(j + n - 1)], collapse = "-")
          subsequences <- c(subsequences, subseq)
        }
      }
    }
    return(subsequences)
  }
  
  analyze_discrimination <- function(seq_A, seq_B, n, min_char_length_state_local) {
    # Extract subsequences
    subseq_A <- extract_subsequences(seq_A, n, min_char_length_state_local)
    subseq_B <- extract_subsequences(seq_B, n, min_char_length_state_local)
    
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
  
  analyze_statistical <- function(seq_A, seq_B, n, correction_method, test_method, min_expected, min_char_length_state_local) {
    # Extract subsequences
    subseq_A <- extract_subsequences(seq_A, n, min_char_length_state_local)
    subseq_B <- extract_subsequences(seq_B, n, min_char_length_state_local)
    
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
        seq_A, seq_B, n, correction, test_method, min_expected, min_char_length_state_local = min_char_length_state)
    } else {
      analysis_results[[paste0("length_", n)]] <- analyze_discrimination(
        seq_A, seq_B, n, min_char_length_state_local = min_char_length_state)
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
        # For statistical analysis, use standardized residuals
        residual_data <- patterns[1:min(nrow(patterns), parameters$top_n), ]
        
        # Calculate standardized residuals
        # These totals are for the 'top_n' patterns being plotted
        total_freq_A_top_n <- sum(residual_data$freq_A)
        total_freq_B_top_n <- sum(residual_data$freq_B)
        grand_total_top_n <- total_freq_A_top_n + total_freq_B_top_n
        
        if (grand_total_top_n == 0) { # Avoid division by zero if no patterns
            resid_A <- rep(0, nrow(residual_data))
            resid_B <- rep(0, nrow(residual_data))
        } else {
            expected_A_list <- numeric(nrow(residual_data))
            expected_B_list <- numeric(nrow(residual_data))

            for (k_pattern in 1:nrow(residual_data)) {
                pattern_total_freq <- residual_data$freq_A[k_pattern] + residual_data$freq_B[k_pattern]
                expected_A_list[k_pattern] <- (pattern_total_freq * total_freq_A_top_n) / grand_total_top_n
                expected_B_list[k_pattern] <- (pattern_total_freq * total_freq_B_top_n) / grand_total_top_n
            }

            # Standardized residuals: (Observed - Expected) / sqrt(Expected)
            # Adding 0.5 to Expected in denominator to avoid issues with Expected = 0
            resid_A <- (residual_data$freq_A - expected_A_list) / sqrt(expected_A_list + 0.5)
            resid_B <- (residual_data$freq_B - expected_B_list) / sqrt(expected_B_list + 0.5)
        }
        residual_matrix <- cbind(resid_A, resid_B)
      } else {
        # For discrimination analysis, use proportion differences
        residual_data <- patterns[1:min(nrow(patterns), parameters$top_n), ]
        if (nrow(residual_data) > 0 && "prop_diff" %in% names(residual_data)) {
            prop_diff <- residual_data$prop_diff
            residual_matrix <- cbind(prop_diff, -prop_diff)
        } else {
            # Fallback if prop_diff is not available or no data
            residual_matrix <- matrix(0, nrow = nrow(residual_data), ncol = 2)
            if (nrow(residual_data) > 0) {
                 rownames(residual_matrix) <- residual_data$pattern
            }
        }
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
      min_expected = min_expected,
      min_char_length_state = min_char_length_state # Added to parameters
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
        # For statistical analysis, use standardized residuals
        residual_data <- patterns[1:min(nrow(patterns), parameters$top_n), ]
        
        # Calculate standardized residuals
        # These totals are for the 'top_n' patterns being plotted
        total_freq_A_top_n <- sum(residual_data$freq_A)
        total_freq_B_top_n <- sum(residual_data$freq_B)
        grand_total_top_n <- total_freq_A_top_n + total_freq_B_top_n

        if (grand_total_top_n == 0) { # Avoid division by zero if no patterns
            resid_A <- rep(0, nrow(residual_data))
            resid_B <- rep(0, nrow(residual_data))
        } else {
            expected_A_list <- numeric(nrow(residual_data))
            expected_B_list <- numeric(nrow(residual_data))

            for (k_pattern in 1:nrow(residual_data)) {
                pattern_total_freq <- residual_data$freq_A[k_pattern] + residual_data$freq_B[k_pattern]
                expected_A_list[k_pattern] <- (pattern_total_freq * total_freq_A_top_n) / grand_total_top_n
                expected_B_list[k_pattern] <- (pattern_total_freq * total_freq_B_top_n) / grand_total_top_n
            }

            # Standardized residuals: (Observed - Expected) / sqrt(Expected)
            # Adding 0.5 to Expected in denominator to avoid issues with Expected = 0
            resid_A <- (residual_data$freq_A - expected_A_list) / sqrt(expected_A_list + 0.5)
            resid_B <- (residual_data$freq_B - expected_B_list) / sqrt(expected_B_list + 0.5)
        }
        residual_matrix <- cbind(resid_A, resid_B)
      } else {
        # For discrimination analysis, use proportion differences
        residual_data <- patterns[1:min(nrow(patterns), parameters$top_n), ]
        if (nrow(residual_data) > 0 && "prop_diff" %in% names(residual_data)) {
            prop_diff <- residual_data$prop_diff
            residual_matrix <- cbind(prop_diff, -prop_diff)
        } else {
            # Fallback if prop_diff is not available or no data
            residual_matrix <- matrix(0, nrow = nrow(residual_data), ncol = 2)
            if (nrow(residual_data) > 0) {
                 rownames(residual_matrix) <- residual_data$pattern
            }
        }
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
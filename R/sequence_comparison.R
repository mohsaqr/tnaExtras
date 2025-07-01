# ==============================================================================
# SEQUENCE COMPARISON TOOLKIT
# ==============================================================================
# 
# Advanced sequence comparison functions with statistical testing capabilities.
# Includes chi-square tests, Fisher's exact tests, and comprehensive pattern
# analysis with visualization support.
#
# Functions:
# - compare_sequences(): Main sequence comparison function with statistical testing
# - compare_sequences_multi(): Multi-group sequence comparison
# - Helper functions for statistical analysis and visualization
#
# ==============================================================================

# Pattern analysis functions are available within the package

# ==============================================================================
# MULTI-GROUP SEQUENCE COMPARISON FUNCTION
# ==============================================================================

#' Compare Sequences Across Multiple Groups
#'
#' Performs comprehensive sequence comparison analysis across multiple groups, identifying
#' discriminating subsequences using support-based measures.
#' Supports both data.frame input and group_tna objects from the tna package.
#'
#' @param data A data frame containing sequence data in wide format, where each row
#'   represents a sequence and each column represents a time point OR a group_tna object
#' @param group Column name or index containing group information, or a vector indicating 
#'   group membership for each sequence. If data contains a group column, specify column 
#'   name (e.g., "Group") or pass as vector (e.g., data$Group). Ignored if data is a group_tna object.
#' @param min_length Minimum subsequence length to analyze (default: 2)
#' @param max_length Maximum subsequence length to analyze (default: 5)
#' @param top_n Number of top patterns to return and display (default: 10)
#' @param min_frequency Minimum frequency required to include a pattern (default: 2)
#'
#' @return A compare_sequences_multi object containing:
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
#' # Multi-group analysis with data.frame
#' result <- compare_sequences_multi(seqdata, "Group")
#' 
#' # Multi-group analysis with group_tna object
#' # group_tna_obj <- tna::create_group_tna(...)
#' # result <- compare_sequences_multi(group_tna_obj)
#' print(result)
#' }
#'
compare_sequences_multi_internal <- function(data, group, min_length = 2, max_length = 5, 
                                             top_n = 10, min_frequency = 2) {
  
  # =====================================================================
  # INPUT VALIDATION AND GROUP_TNA SUPPORT
  # =====================================================================
  
  # Check if input is a group_tna object
  if (is_group_tna(data)) {
    cat("Detected group_tna object, converting to tnaExtras format...\n")
    converted <- convert_group_tna(data)
    data <- converted$data
    group <- converted$group_col
    group_info <- converted$group_info
    
    cat("Successfully converted group_tna object:\n")
    cat("  Label:", group_info$label %||% "Unknown", "\n")
    cat("  Groups:", paste(group_info$levels, collapse = ", "), "\n")
    cat("  Total sequences:", nrow(data), "\n")
  } else {
    group_info <- NULL
  }
  
  # Validate inputs
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame or group_tna object", call. = FALSE)
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
  
  if (length(groups) < 2) {
    stop("Group variable must have at least 2 levels, found: ", length(groups), call. = FALSE)
  }
  
  # Split data by groups
  group_data <- list()
  for (g in groups) {
    group_seq <- data[group_factor == g, , drop = FALSE]
    # Clean data: remove empty rows and convert empty strings to NA
    group_seq <- group_seq[apply(group_seq, 1, function(x) !all(is.na(x) | x == "")), , drop = FALSE]
    if (nrow(group_seq) == 0) {
      stop("Group '", g, "' has no valid sequences", call. = FALSE)
    }
    group_seq[group_seq == ""] <- NA
    group_data[[g]] <- group_seq
  }
  
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
  
  analyze_multi_group <- function(group_data, n) {
    # Extract subsequences for each group
    group_subseqs <- list()
    for (g in names(group_data)) {
      group_subseqs[[g]] <- extract_subsequences(group_data[[g]], n)
    }
    
    # Get all unique patterns
    all_patterns <- unique(unlist(group_subseqs))
    
    if (length(all_patterns) == 0) {
      return(list(n = n, patterns = data.frame(), n_patterns = 0))
    }
    
    # Create frequency table
    patterns <- data.frame(
      pattern = all_patterns,
      stringsAsFactors = FALSE
    )
    
    # Add frequency columns for each group
    for (g in names(group_data)) {
      freq_col <- paste0("freq_", g)
      prop_col <- paste0("prop_", g)
      
      freq_table <- table(group_subseqs[[g]])
      patterns[[freq_col]] <- as.numeric(freq_table[all_patterns])
      patterns[[freq_col]][is.na(patterns[[freq_col]])] <- 0
      
      total_freq <- sum(patterns[[freq_col]])
      if (total_freq > 0) {
        patterns[[prop_col]] <- patterns[[freq_col]] / total_freq
      } else {
        patterns[[prop_col]] <- 0
      }
    }
    
    # Calculate discrimination measures
    # Use coefficient of variation (CV) of proportions as discrimination measure
    prop_cols <- paste0("prop_", names(group_data))
    if (length(prop_cols) > 1) {
      patterns$mean_prop <- apply(patterns[prop_cols], 1, mean)
      patterns$sd_prop <- apply(patterns[prop_cols], 1, sd)
      patterns$cv_prop <- ifelse(patterns$mean_prop == 0, 0, patterns$sd_prop / patterns$mean_prop)
      
      # Also calculate range
      patterns$max_prop <- apply(patterns[prop_cols], 1, max)
      patterns$min_prop <- apply(patterns[prop_cols], 1, min)
      patterns$range_prop <- patterns$max_prop - patterns$min_prop
      
      # Find dominant group
      max_indices <- apply(patterns[prop_cols], 1, which.max)
      patterns$dominant_group <- names(group_data)[max_indices]
      
      # Sort by range (most discriminating patterns first)
      patterns <- patterns[order(patterns$range_prop, decreasing = TRUE), ]
    } else {
      patterns$cv_prop <- patterns$range_prop <- 0
      patterns$dominant_group <- names(group_data)[1]
    }
    
    return(list(
      n = n,
      patterns = patterns,
      n_patterns = nrow(patterns)
    ))
  }
  
  # =====================================================================
  # MAIN ANALYSIS
  # =====================================================================
  
  cat("Analyzing sequences across", length(groups), "groups:", paste(groups, collapse = ", "), "\n")
  group_sizes <- sapply(group_data, nrow)
  cat("Group sizes:", paste(paste(groups, group_sizes, sep = ":"), collapse = ", "), "\n")
  cat("Subsequence lengths:", min_length, "to", max_length, "\n\n")
  
  # Run analysis for each subsequence length
  analysis_results <- list()
  for (n in min_length:max_length) {
    analysis_results[[paste0("length_", n)]] <- analyze_multi_group(group_data, n)
  }
  
  # =====================================================================
  # COMPILE RESULTS
  # =====================================================================
  
  # Combine results across all lengths for summary
  all_patterns <- do.call(rbind, lapply(analysis_results, function(x) {
    if (nrow(x$patterns) > 0) {
      x$patterns$length <- x$n
      return(x$patterns)
    } else {
      return(data.frame())
    }
  }))
  
  # Filter by minimum frequency
  if (nrow(all_patterns) > 0) {
    total_freq_col <- paste0("freq_", names(group_data))
    all_patterns$total_freq <- rowSums(all_patterns[total_freq_col])
    all_patterns <- all_patterns[all_patterns$total_freq >= min_frequency, ]
  }
  
  # Get top patterns
  if (nrow(all_patterns) > 0) {
    top_patterns <- head(all_patterns, top_n)
  } else {
    top_patterns <- data.frame()
  }
  
  # =====================================================================
  # SUMMARY STATISTICS
  # =====================================================================
  
  stats <- list(
    n_groups = length(groups),
    group_names = groups,
    group_sizes = group_sizes,
    total_sequences = sum(group_sizes),
    n_patterns_total = nrow(all_patterns),
    n_patterns_by_length = sapply(analysis_results, function(x) x$n_patterns),
    min_length = min_length,
    max_length = max_length,
    min_frequency = min_frequency,
    group_tna_info = group_info  # Store original group_tna metadata
  )
  
  # =====================================================================
  # RETURN RESULTS
  # =====================================================================
  
  result <- list(
    results = analysis_results,
    summary = top_patterns,
    all_patterns = all_patterns,
    stats = stats,
    parameters = list(
      min_length = min_length,
      max_length = max_length,
      top_n = top_n,
      min_frequency = min_frequency
    )
  )
  
  class(result) <- "compare_sequences_multi"
  return(result)
}

# ==============================================================================
# MAIN SEQUENCE COMPARISON FUNCTION
# ==============================================================================

#' Compare Sequences Between Groups
#'
#' Performs comprehensive sequence comparison analysis between groups. Automatically
#' detects the number of groups and uses appropriate analysis methods:
#' - For 2 groups: Full statistical analysis with tests, effect sizes, and detailed measures
#' - For 3+ groups: Multi-group discrimination analysis with support-based measures
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
#'   discrimination score analysis. Only applies to 2-group analysis
#' @param correction Character, multiple comparison correction method (default: "bonferroni").
#'   Supports all R p.adjust methods: "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
#'   Only applies to 2-group statistical analysis
#' @param test_method Character, statistical test selection (default: "auto").
#'   Options: "auto", "fisher", "chi.squared". Only applies to 2-group statistical analysis
#' @param min_expected Numeric, minimum expected count for automatic test selection (default: 5).
#'   Only applies to 2-group statistical analysis
#' @param min_frequency Numeric, minimum frequency for pattern inclusion in multi-group analysis (default: 2)
#'
#' @return For 2 groups: A compare_sequences object with statistical measures
#'         For 3+ groups: A compare_sequences_multi object with discrimination measures
#'
#' @examples
#' \dontrun{
#' # Load example data
#' data(seqdata)
#' 
#' # Automatically detects 2 groups - uses detailed statistical analysis
#' result <- compare_sequences(seqdata, "Group")
#' print(result)
#' 
#' # Automatically detects 3+ groups - uses multi-group discrimination analysis
#' result <- compare_sequences(multigroup_data, "Group")
#' print(result)
#' 
#' # Statistical analysis (only for 2 groups)
#' result <- compare_sequences(seqdata, "Group", statistical = TRUE)
#' }
#'
#' @export
compare_sequences <- function(data, group, min_length = 2, max_length = 5, top_n = 10, 
                             detailed = FALSE, statistical = FALSE, correction = "bonferroni", 
                             test_method = "auto", min_expected = 5, min_frequency = 2) {
  
  # =====================================================================
  # INPUT VALIDATION AND GROUP_TNA SUPPORT
  # =====================================================================
  
  # Check if input is a group_tna object
  if (is_group_tna(data)) {
    cat("Detected group_tna object, converting to tnaExtras format...\n")
    converted <- convert_group_tna(data)
    data <- converted$data
    group <- converted$group_col
    group_info <- converted$group_info
    
    cat("Successfully converted group_tna object:\n")
    cat("  Label:", group_info$label %||% "Unknown", "\n")
    cat("  Groups:", paste(group_info$levels, collapse = ", "), "\n")
    cat("  Total sequences:", nrow(data), "\n")
  } else {
    group_info <- NULL
  }
  
  # Validate inputs
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame or group_tna object", call. = FALSE)
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
  
  if (!correction %in% p.adjust.methods) {
    stop("'correction' must be one of: ", paste(p.adjust.methods, collapse = ", "), call. = FALSE)
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
  
  # Auto-detection logic: delegate to appropriate function based on number of groups
  if (length(groups) > 2) {
    cat("Detected", length(groups), "groups. Using multi-group analysis...\n")
    # Reconstruct the original data with group column for delegation
    original_data <- data
    original_data[[if(is.character(group) && length(group) == 1) group else "Group"]] <- group_vector
    return(compare_sequences_multi_internal(data = original_data,
                                           group = if(is.character(group) && length(group) == 1) group else "Group",
                                           min_length = min_length, max_length = max_length,
                                           top_n = top_n, min_frequency = min_frequency))
  } else if (length(groups) != 2) {
    stop("Group variable must have exactly 2 levels for two-group analysis, found: ", length(groups), call. = FALSE)
  }
  
  # Continue with 2-group analysis
  cat("Detected 2 groups. Using detailed two-group analysis...\n")
  
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
  
  analyze_discrimination <- function(seq_A, seq_B, n, group_names = c("A", "B")) {
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
    
    # Create dynamic column names using actual group names
    freq_col_A <- paste0("freq_", group_names[1])
    freq_col_B <- paste0("freq_", group_names[2])
    prop_col_A <- paste0("prop_", group_names[1])
    prop_col_B <- paste0("prop_", group_names[2])
    
    # Create frequency table
    patterns <- data.frame(
      pattern = all_patterns,
      stringsAsFactors = FALSE
    )
    
    # Add frequency columns with actual group names
    patterns[[freq_col_A]] <- as.numeric(freq_A[all_patterns])
    patterns[[freq_col_B]] <- as.numeric(freq_B[all_patterns])
    patterns[is.na(patterns)] <- 0
    
    # Calculate metrics
    total_A <- sum(patterns[[freq_col_A]])
    total_B <- sum(patterns[[freq_col_B]])
    
    if (total_A > 0 && total_B > 0) {
      patterns[[prop_col_A]] <- patterns[[freq_col_A]] / total_A
      patterns[[prop_col_B]] <- patterns[[freq_col_B]] / total_B
      patterns$prop_diff <- patterns[[prop_col_A]] - patterns[[prop_col_B]]
      
      # Discrimination score (simple and robust)
      patterns$discrimination <- abs(patterns$prop_diff) * sqrt(patterns[[freq_col_A]] + patterns[[freq_col_B]])
      
      # Sort by discrimination
      patterns <- patterns[order(patterns$discrimination, decreasing = TRUE), ]
    } else {
      patterns[[prop_col_A]] <- patterns[[prop_col_B]] <- patterns$prop_diff <- patterns$discrimination <- 0
    }
    
    return(list(
      n = n,
      patterns = patterns,
      total_A = total_A,
      total_B = total_B,
      n_patterns = nrow(patterns)
    ))
  }
  
  analyze_statistical <- function(seq_A, seq_B, n, correction_method, test_method, min_expected, group_names = c("A", "B")) {
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
    
    # Create dynamic column names using actual group names
    freq_col_A <- paste0("freq_", group_names[1])
    freq_col_B <- paste0("freq_", group_names[2])
    
    # Create results table
    results <- data.frame(
      pattern = all_patterns,
      stringsAsFactors = FALSE
    )
    
    # Add frequency columns with actual group names
    results[[freq_col_A]] <- as.numeric(freq_A[all_patterns])
    results[[freq_col_B]] <- as.numeric(freq_B[all_patterns])
    results[is.na(results)] <- 0
    
    # Perform statistical tests
    results$p_value <- NA
    results$test_used <- NA
    results$statistic <- NA
    
    for (i in seq_len(nrow(results))) {
      # Create contingency table
      cont_table <- matrix(c(results[[freq_col_A]][i], results[[freq_col_B]][i],
                            sum(results[[freq_col_A]]) - results[[freq_col_A]][i],
                            sum(results[[freq_col_B]]) - results[[freq_col_B]][i]), 
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
        seq_A, seq_B, n, correction, test_method, min_expected, groups)
    } else {
      analysis_results[[paste0("length_", n)]] <- analyze_discrimination(seq_A, seq_B, n, groups)
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
      
      # Detect frequency column names dynamically
      freq_cols <- grep("^freq_", names(patterns), value = TRUE)
      prop_cols <- grep("^prop_", names(patterns), value = TRUE)
      
      # Ensure we use exactly top_n patterns (or fewer if not available)
      n_patterns_to_show <- min(nrow(patterns), parameters$top_n)
      residual_data <- patterns[1:n_patterns_to_show, , drop = FALSE]
      
      # Create residual matrix for heatmap
      if (parameters$statistical) {
        # For statistical analysis, use adjusted standardized residuals for better visualization
        # Calculate better residuals for visualization using dynamic column names
        if (length(freq_cols) >= 2) {
          total_A <- sum(residual_data[[freq_cols[1]]])
          total_B <- sum(residual_data[[freq_cols[2]]])
          total_overall <- total_A + total_B
          
          # Calculate proportions within each group
          prop_A <- residual_data[[freq_cols[1]]] / total_A
          prop_B <- residual_data[[freq_cols[2]]] / total_B
          
          # Use standardized proportion differences (z-score like)
          # This shows how much each pattern deviates from equal representation
          overall_prop <- (residual_data[[freq_cols[1]]] + residual_data[[freq_cols[2]]]) / total_overall
          
          # Calculate standard error for proportion difference
          se_A <- sqrt(overall_prop * (1 - overall_prop) / total_A)
          se_B <- sqrt(overall_prop * (1 - overall_prop) / total_B)
          
          # Standardized deviations from expected proportion
          resid_A <- (prop_A - overall_prop) / se_A
          resid_B <- (prop_B - overall_prop) / se_B
          
          residual_matrix <- cbind(resid_A, resid_B)
        } else {
          # Fallback for insufficient data
          residual_matrix <- matrix(0, nrow = nrow(residual_data), ncol = 2)
        }
      } else {
        # For discrimination analysis, use proportion differences
        prop_diff <- residual_data$prop_diff
        residual_matrix <- cbind(prop_diff, -prop_diff)
      }
      
      rownames(residual_matrix) <- residual_data$pattern
      colnames(residual_matrix) <- groups
      
      # Calculate dynamic margins based on pattern name lengths
      max_pattern_length <- max(nchar(rownames(residual_matrix)), na.rm = TRUE)
      left_margin <- max(12, ceiling(max_pattern_length * 0.6))  # Dynamic left margin
      
      # Save current graphics parameters
      old_par <- par(no.readonly = TRUE)
      on.exit(par(old_par))
      
      # Set up layout for main plot + legend with better proportions
      layout(matrix(c(1, 2), nrow = 1), widths = c(8, 2))  # Thin but visible legend
      
      # Main heatmap plot with better margins
      par(mar = c(6, left_margin, 4, 1))
      
      # Color palette (reversed: red-white-blue)
      max_val <- max(abs(residual_matrix), na.rm = TRUE)
      if (max_val == 0) max_val <- 1  # Avoid division by zero
      colors <- colorRampPalette(c("#900C3F", "white", "#043996"))(100)
      
      # Create coordinate vectors for proper alignment
      x_coords <- seq_len(ncol(residual_matrix))
      y_coords <- seq_len(nrow(residual_matrix))
      
      # Create heatmap with proper coordinate system
      image(x = x_coords, y = y_coords, z = t(residual_matrix), 
            col = colors,
            breaks = seq(-max_val, max_val, length.out = 101),
            xlab = "Groups", ylab = "",
            axes = FALSE, 
            main = "")
      
      # Add properly aligned axes
      axis(1, at = x_coords, labels = colnames(residual_matrix), 
           tick = TRUE, line = 0)
      axis(2, at = y_coords, labels = rownames(residual_matrix), 
           las = 2, tick = TRUE, line = 0, cex.axis = 0.8)
      
      # Add grid lines for better readability
      abline(h = y_coords + 0.5, col = "lightgray", lwd = 0.5)
      abline(v = x_coords + 0.5, col = "lightgray", lwd = 0.5)
      
      # Add a box around the plot
      box()
      
      # Legend plot with margins
      par(mar = c(6, 1, 4, 3))
      
      # Create gradient legend
      legend_y <- seq(0, 1, length.out = 100)
      legend_matrix <- matrix(legend_y, ncol = 1)
      
      image(x = 1, y = legend_y, z = t(legend_matrix), 
            col = colors, axes = FALSE, xlab = "", ylab = "",
            main = "")
      
      # Add legend labels with better positioning
      legend_vals <- seq(-max_val, max_val, length.out = 5)
      legend_positions <- seq(0, 1, length.out = 5)
      axis(4, at = legend_positions, labels = sprintf("%.2f", legend_vals), 
           las = 2, cex.axis = 0.8)
      
      # Add a box around the legend
      box()
      
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
            cat(sprintf("[OK] Created %s-length subsequences plot\n", length_num))
          }, error = function(e) {
            cat(sprintf("[WARN] Could not create %s-length plot: %s\n", length_num, e$message))
          })
        }
      }
    } else if (is.data.frame(summary_output) && nrow(summary_output) > 0) {
      # Create combined plot
      tryCatch({
        create_heatmap(summary_output, groups)
        cat("[OK] Created combined discriminating patterns plot\n")
      }, error = function(e) {
        cat("[WARN] Could not create combined plot:", e$message, "\n")
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
      group_tna_info = group_info  # Store original group_tna metadata
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
      
      # Detect frequency column names dynamically
      freq_cols <- grep("^freq_", names(patterns), value = TRUE)
      prop_cols <- grep("^prop_", names(patterns), value = TRUE)
      
      # Ensure we use exactly top_n patterns (or fewer if not available)
      n_patterns_to_show <- min(nrow(patterns), parameters$top_n)
      residual_data <- patterns[1:n_patterns_to_show, , drop = FALSE]
      
      # Create residual matrix for heatmap
      if (parameters$statistical) {
        # For statistical analysis, use adjusted standardized residuals for better visualization
        # Calculate better residuals for visualization using dynamic column names
        if (length(freq_cols) >= 2) {
          total_A <- sum(residual_data[[freq_cols[1]]])
          total_B <- sum(residual_data[[freq_cols[2]]])
          total_overall <- total_A + total_B
          
          # Calculate proportions within each group
          prop_A <- residual_data[[freq_cols[1]]] / total_A
          prop_B <- residual_data[[freq_cols[2]]] / total_B
          
          # Use standardized proportion differences (z-score like)
          # This shows how much each pattern deviates from equal representation
          overall_prop <- (residual_data[[freq_cols[1]]] + residual_data[[freq_cols[2]]]) / total_overall
          
          # Calculate standard error for proportion difference
          se_A <- sqrt(overall_prop * (1 - overall_prop) / total_A)
          se_B <- sqrt(overall_prop * (1 - overall_prop) / total_B)
          
          # Standardized deviations from expected proportion
          resid_A <- (prop_A - overall_prop) / se_A
          resid_B <- (prop_B - overall_prop) / se_B
          
          residual_matrix <- cbind(resid_A, resid_B)
        } else {
          # Fallback for insufficient data
          residual_matrix <- matrix(0, nrow = nrow(residual_data), ncol = 2)
        }
      } else {
        # For discrimination analysis, use proportion differences
        prop_diff <- residual_data$prop_diff
        residual_matrix <- cbind(prop_diff, -prop_diff)
      }
      
      rownames(residual_matrix) <- residual_data$pattern
      colnames(residual_matrix) <- groups
      
      # Calculate dynamic margins based on pattern name lengths
      max_pattern_length <- max(nchar(rownames(residual_matrix)), na.rm = TRUE)
      left_margin <- max(12, ceiling(max_pattern_length * 0.6))  # Dynamic left margin
      
      # Save current graphics parameters
      old_par <- par(no.readonly = TRUE)
      on.exit(par(old_par))
      
      # Set up layout for main plot + legend with better proportions
      layout(matrix(c(1, 2), nrow = 1), widths = c(8, 2))  # Thin but visible legend
      
      # Main heatmap plot with better margins
      par(mar = c(6, left_margin, 4, 1))
      
      # Color palette (reversed: red-white-blue)
      max_val <- max(abs(residual_matrix), na.rm = TRUE)
      if (max_val == 0) max_val <- 1  # Avoid division by zero
      colors <- colorRampPalette(c("red", "white", "blue"))(100)
      
      # Create coordinate vectors for proper alignment
      x_coords <- seq_len(ncol(residual_matrix))
      y_coords <- seq_len(nrow(residual_matrix))
      
      # Create heatmap with proper coordinate system
      image(x = x_coords, y = y_coords, z = t(residual_matrix), 
            col = colors,
            breaks = seq(-max_val, max_val, length.out = 101),
            xlab = "Groups", ylab = "",
            axes = FALSE, 
            main = "")
      
      # Add properly aligned axes
      axis(1, at = x_coords, labels = colnames(residual_matrix), 
           tick = TRUE, line = 0)
      axis(2, at = y_coords, labels = rownames(residual_matrix), 
           las = 2, tick = TRUE, line = 0, cex.axis = 0.8)
      
      # Add grid lines for better readability
      abline(h = y_coords + 0.5, col = "lightgray", lwd = 0.5)
      abline(v = x_coords + 0.5, col = "lightgray", lwd = 0.5)
      
      # Add a box around the plot
      box()
      
      # Legend plot with margins
      par(mar = c(6, 1, 4, 3))
      
      # Create gradient legend
      legend_y <- seq(0, 1, length.out = 100)
      legend_matrix <- matrix(legend_y, ncol = 1)
      
      image(x = 1, y = legend_y, z = t(legend_matrix), 
            col = colors, axes = FALSE, xlab = "", ylab = "",
            main = "")
      
      # Add legend labels with better positioning
      legend_vals <- seq(-max_val, max_val, length.out = 5)
      legend_positions <- seq(0, 1, length.out = 5)
      axis(4, at = legend_positions, labels = sprintf("%.2f", legend_vals), 
           las = 2, cex.axis = 0.8)
      
      # Add a box around the legend
      box()
      
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
            cat(sprintf("[OK] Created %s-length subsequences plot\n", length_num))
          }, error = function(e) {
            cat(sprintf("[WARN] Could not create %s-length plot: %s\n", length_num, e$message))
          })
        }
      }
    } else if (is.data.frame(summary_output) && nrow(summary_output) > 0) {
      # Create combined plot
      tryCatch({
        create_heatmap(summary_output, groups)
        cat("[OK] Created combined discriminating patterns plot\n")
      }, error = function(e) {
        cat("[WARN] Could not create combined plot:", e$message, "\n")
      })
    }
  }
  
  # Create plots
  create_plots(x$summary, x$parameters, x$parameters$groups)
  cat("Plot recreation complete!\n")
}

# ==============================================================================
# PRINT AND SUMMARY METHODS FOR MULTI-GROUP COMPARISON
# ==============================================================================

#' Print method for compare_sequences_multi objects
#' @param x compare_sequences_multi object
#' @param ... additional arguments
print.compare_sequences_multi <- function(x, ...) {
  cat("Multi-Group Sequence Comparison Results\n")
  cat("=======================================\n\n")
  
  cat("Analysis Summary:\n")
  cat("  Groups:", paste(x$stats$group_names, collapse = ", "), "\n")
  cat("  Group sizes:", paste(paste(x$stats$group_names, x$stats$group_sizes, sep = ":"), collapse = ", "), "\n")
  cat("  Total sequences:", x$stats$total_sequences, "\n")
  cat("  Subsequence lengths:", x$parameters$min_length, "to", x$parameters$max_length, "\n")
  cat("  Total patterns found:", x$stats$n_patterns_total, "\n")
  cat("  Patterns by length:\n")
  for (i in seq_along(x$stats$n_patterns_by_length)) {
    len <- names(x$stats$n_patterns_by_length)[i]
    count <- x$stats$n_patterns_by_length[i]
    cat("    ", len, ":", count, "patterns\n")
  }
  
  cat("\nTop discriminating patterns:\n")
  if (nrow(x$summary) > 0) {
    n_show <- min(5, nrow(x$summary))
    show_cols <- c("pattern", "length", "range_prop", "dominant_group")
    available_cols <- intersect(show_cols, names(x$summary))
    print(head(x$summary[available_cols], n_show))
    
    if (nrow(x$summary) > n_show) {
      cat("  ... and", nrow(x$summary) - n_show, "more patterns\n")
    }
  } else {
    cat("  No patterns found meeting criteria\n")
  }
  
  cat("\nUse summary(result) for detailed analysis\n")
}

#' Summary method for compare_sequences_multi objects
#' @param object compare_sequences_multi object
#' @param length which subsequence length to summarize (default: all)
#' @param top_n number of top patterns to show (default: 10)
#' @param ... additional arguments
summary.compare_sequences_multi <- function(object, length = NULL, top_n = 10, ...) {
  cat("Multi-Group Sequence Comparison - Detailed Summary\n")
  cat("==================================================\n\n")
  
  # Overall statistics
  cat("Dataset Information:\n")
  cat("  Groups:", length(object$stats$group_names), "(", paste(object$stats$group_names, collapse = ", "), ")\n")
  cat("  Total sequences:", object$stats$total_sequences, "\n")
  cat("  Subsequence range:", object$parameters$min_length, "to", object$parameters$max_length, "\n")
  cat("  Minimum frequency:", object$parameters$min_frequency, "\n\n")
  
  # Group-specific information
  cat("Group Details:\n")
  for (i in seq_along(object$stats$group_names)) {
    g <- object$stats$group_names[i]
    size <- object$stats$group_sizes[i]
    prop <- round(size / object$stats$total_sequences * 100, 1)
    cat("  ", g, ":", size, "sequences (", prop, "%)\n")
  }
  cat("\n")
  
  # Pattern summary by length
  if (is.null(length)) {
    cat("Pattern Summary by Length:\n")
    for (len_key in names(object$results)) {
      len_num <- object$results[[len_key]]$n
      n_patterns <- object$results[[len_key]]$n_patterns
      cat("  Length", len_num, ":", n_patterns, "patterns\n")
    }
    cat("\n")
    
    # Show top patterns across all lengths
    cat("Top", min(top_n, nrow(object$summary)), "Discriminating Patterns (All Lengths):\n")
    if (nrow(object$summary) > 0) {
      print(head(object$summary, top_n))
    } else {
      cat("  No patterns found\n")
    }
    
  } else {
    # Show specific length
    len_key <- paste0("length_", length)
    if (len_key %in% names(object$results)) {
      len_data <- object$results[[len_key]]
      cat("Patterns for Length", length, ":\n")
      cat("  Total patterns:", len_data$n_patterns, "\n\n")
      
      if (len_data$n_patterns > 0) {
        n_show <- min(top_n, nrow(len_data$patterns))
        cat("Top", n_show, "patterns:\n")
        print(head(len_data$patterns, n_show))
      } else {
        cat("  No patterns found for this length\n")
      }
    } else {
      cat("Length", length, "not found in results\n")
    }
  }
}

# ==============================================================================
# ALIASES FOR BACKWARD COMPATIBILITY
# ==============================================================================

#' Multi-Group Sequence Comparison (Alias for compare_sequences)
#'
#' This function is now an alias for compare_sequences() which auto-detects
#' the number of groups. For multi-group analysis, use compare_sequences() directly.
#' 
#' @param ... All arguments passed to compare_sequences()
#' @return Same as compare_sequences()
#' @export
compare_sequences_multi <- function(...) {
  # Force multi-group behavior if not enough groups are auto-detected
  args <- list(...)
  result <- do.call(compare_sequences, args)
  return(result)
} 

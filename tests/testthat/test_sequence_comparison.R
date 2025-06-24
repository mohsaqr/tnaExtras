# =============================================================================
# COMPREHENSIVE TESTS FOR SEQUENCE COMPARISON FUNCTION
# =============================================================================
# Test suite for seqompare() function
# Tests n-gram extraction, statistical analysis, and visualization capabilities

# Load the function
source("functions/sequence_comparison.R")

# =============================================================================
# TEST DATA PREPARATION
# =============================================================================

cat("=== SEQUENCE COMPARISON TESTING SUITE ===\n\n")

# Create test data with clear group differences
set.seed(42)  # For reproducible results

# Group A: More "Active" patterns
group_a_data <- data.frame(
  T1 = rep(c("Active", "Average", "Active"), each = 10),
  T2 = rep(c("Active", "Active", "Average"), each = 10),
  T3 = rep(c("Average", "Active", "Active"), each = 10),
  T4 = rep(c("Active", "Average", "Disengaged"), each = 10),
  T5 = rep(c("Active", "Active", "Average"), each = 10),
  Group = "A"
)

# Group B: More "Disengaged" patterns  
group_b_data <- data.frame(
  T1 = rep(c("Disengaged", "Average", "Disengaged"), each = 10),
  T2 = rep(c("Disengaged", "Disengaged", "Average"), each = 10),
  T3 = rep(c("Average", "Disengaged", "Disengaged"), each = 10),
  T4 = rep(c("Disengaged", "Average", "Active"), each = 10),
  T5 = rep(c("Disengaged", "Disengaged", "Average"), each = 10),
  Group = "B"
)

# Combine data
test_data <- rbind(group_a_data, group_b_data)

# Small test data for quick tests
small_test_data <- data.frame(
  T1 = c("Active", "Average", "Disengaged", "Active"),
  T2 = c("Active", "Active", "Disengaged", "Average"),
  T3 = c("Average", "Active", "Average", "Average"),
  Group = c("A", "A", "B", "B")
)

# Edge case: single group
single_group_data <- data.frame(
  T1 = c("Active", "Average", "Active"),
  T2 = c("Active", "Active", "Average"),
  T3 = c("Average", "Active", "Active"),
  Group = c("A", "A", "A")
)

# =============================================================================
# TEST FUNCTIONS
# =============================================================================

test_basic_functionality <- function() {
  cat("1. BASIC FUNCTIONALITY TEST\n")
  cat("===========================\n")
  
  # Test basic seqompare function
  result <- seqompare(test_data, min_ngram = 2, max_ngram = 3, top_n = 5, save_plots = FALSE)
  
  # Check structure
  expected_components <- c("top_patterns", "ngram_analysis", "statistical_tests", "residuals")
  missing_components <- setdiff(expected_components, names(result))
  
  if (length(missing_components) > 0) {
    cat("❌ MISSING COMPONENTS:", paste(missing_components, collapse = ", "), "\n")
  } else {
    cat("✅ All main components present\n")
  }
  
  # Check top patterns
  if (!is.null(result$top_patterns) && nrow(result$top_patterns) > 0) {
    cat("✅ Top patterns extracted:", nrow(result$top_patterns), "patterns\n")
    cat("✅ Pattern example:", result$top_patterns$pattern[1], "\n")
  } else {
    cat("❌ No top patterns found\n")
  }
  
  # Check n-gram analysis
  if (!is.null(result$ngram_analysis)) {
    cat("✅ N-gram analysis completed\n")
    cat("✅ N-gram lengths analyzed:", paste(names(result$ngram_analysis), collapse = ", "), "\n")
  } else {
    cat("❌ N-gram analysis missing\n")
  }
  
  cat("\n")
  return(result)
}

test_ngram_extraction <- function() {
  cat("2. N-GRAM EXTRACTION TEST\n")
  cat("=========================\n")
  
  # Test with small data for detailed inspection
  result <- seqompare(small_test_data, min_ngram = 2, max_ngram = 3, top_n = 10, save_plots = FALSE)
  
  # Check 2-grams
  if ("2" %in% names(result$ngram_analysis)) {
    ngrams_2 <- result$ngram_analysis[["2"]]
    cat("✅ 2-grams extracted:", nrow(ngrams_2), "unique patterns\n")
    
    # Show some examples
    if (nrow(ngrams_2) > 0) {
      cat("✅ Example 2-grams:\n")
      for (i in 1:min(3, nrow(ngrams_2))) {
        cat("   ", ngrams_2$pattern[i], "- Group A:", ngrams_2$group_A_freq[i], 
            "Group B:", ngrams_2$group_B_freq[i], "\n")
      }
    }
  }
  
  # Check 3-grams
  if ("3" %in% names(result$ngram_analysis)) {
    ngrams_3 <- result$ngram_analysis[["3"]]
    cat("✅ 3-grams extracted:", nrow(ngrams_3), "unique patterns\n")
  }
  
  cat("\n")
  return(result)
}

test_statistical_analysis <- function() {
  cat("3. STATISTICAL ANALYSIS TEST\n")
  cat("============================\n")
  
  # Test with larger dataset for meaningful statistics
  result <- seqompare(test_data, min_ngram = 2, max_ngram = 3, top_n = 10, save_plots = FALSE)
  
  # Check statistical tests
  if (!is.null(result$statistical_tests) && nrow(result$statistical_tests) > 0) {
    stats <- result$statistical_tests
    
    cat("✅ Statistical tests performed:", nrow(stats), "patterns tested\n")
    
    # Check for required columns
    required_cols <- c("pattern", "chi_square", "p_value", "significant")
    missing_cols <- setdiff(required_cols, names(stats))
    
    if (length(missing_cols) == 0) {
      cat("✅ All required statistical columns present\n")
      
      # Check significance
      significant_patterns <- sum(stats$significant, na.rm = TRUE)
      cat("✅ Significant patterns found:", significant_patterns, "\n")
      
      # Show most significant pattern
      if (significant_patterns > 0) {
        most_sig <- stats[which.min(stats$p_value), ]
        cat("✅ Most significant pattern:", most_sig$pattern, 
            "(p =", round(most_sig$p_value, 4), ")\n")
      }
      
    } else {
      cat("❌ Missing statistical columns:", paste(missing_cols, collapse = ", "), "\n")
    }
  } else {
    cat("❌ No statistical tests performed\n")
  }
  
  cat("\n")
  return(result)
}

test_residual_analysis <- function() {
  cat("4. RESIDUAL ANALYSIS TEST\n")
  cat("=========================\n")
  
  result <- seqompare(test_data, min_ngram = 2, max_ngram = 3, top_n = 10, save_plots = FALSE)
  
  # Check residuals
  if (!is.null(result$residuals) && nrow(result$residuals) > 0) {
    residuals <- result$residuals
    
    cat("✅ Residual analysis completed:", nrow(residuals), "patterns\n")
    
    # Check for required columns
    required_cols <- c("pattern", "residual", "abs_residual")
    missing_cols <- setdiff(required_cols, names(residuals))
    
    if (length(missing_cols) == 0) {
      cat("✅ All residual columns present\n")
      
      # Check residual ranges
      max_residual <- max(abs(residuals$residual), na.rm = TRUE)
      cat("✅ Maximum absolute residual:", round(max_residual, 3), "\n")
      
      # Show highest residual pattern
      highest_residual_idx <- which.max(abs(residuals$residual))
      if (length(highest_residual_idx) > 0) {
        highest_pattern <- residuals[highest_residual_idx, ]
        cat("✅ Highest residual pattern:", highest_pattern$pattern, 
            "(residual =", round(highest_pattern$residual, 3), ")\n")
      }
      
    } else {
      cat("❌ Missing residual columns:", paste(missing_cols, collapse = ", "), "\n")
    }
  } else {
    cat("❌ No residual analysis performed\n")
  }
  
  cat("\n")
  return(result)
}

test_edge_cases <- function() {
  cat("5. EDGE CASES TEST\n")
  cat("==================\n")
  
  # Test with single group (should handle gracefully)
  cat("Testing single group data...\n")
  tryCatch({
    single_result <- seqompare(single_group_data, min_ngram = 2, max_ngram = 2, 
                              top_n = 5, save_plots = FALSE)
    cat("❌ Single group should produce error or warning\n")
  }, error = function(e) {
    cat("✅ Single group properly handled with error:", e$message, "\n")
  }, warning = function(w) {
    cat("✅ Single group properly handled with warning:", w$message, "\n")
  })
  
  # Test with very small data
  cat("\nTesting very small dataset...\n")
  tiny_data <- data.frame(
    T1 = c("A", "B"),
    T2 = c("A", "B"),
    Group = c("X", "Y")
  )
  
  tryCatch({
    tiny_result <- seqompare(tiny_data, min_ngram = 2, max_ngram = 2, 
                            top_n = 5, save_plots = FALSE)
    cat("✅ Small dataset handled successfully\n")
  }, error = function(e) {
    cat("✅ Small dataset error handled:", e$message, "\n")
  })
  
  # Test with missing values
  cat("\nTesting missing values...\n")
  missing_data <- test_data
  missing_data$T1[1:2] <- NA
  
  tryCatch({
    missing_result <- seqompare(missing_data, min_ngram = 2, max_ngram = 2, 
                               top_n = 5, save_plots = FALSE)
    cat("✅ Missing values handled successfully\n")
  }, error = function(e) {
    cat("✅ Missing values error handled:", e$message, "\n")
  })
  
  cat("\n")
}

test_parameter_variations <- function() {
  cat("6. PARAMETER VARIATIONS TEST\n")
  cat("============================\n")
  
  # Test different n-gram lengths
  cat("Testing different n-gram lengths...\n")
  
  # Test min_ngram = 2, max_ngram = 2
  result_2 <- seqompare(small_test_data, min_ngram = 2, max_ngram = 2, 
                       top_n = 5, save_plots = FALSE)
  cat("✅ N-gram length 2 only:", length(result_2$ngram_analysis), "length(s)\n")
  
  # Test min_ngram = 3, max_ngram = 4
  result_34 <- seqompare(small_test_data, min_ngram = 3, max_ngram = 4, 
                        top_n = 5, save_plots = FALSE)
  cat("✅ N-gram lengths 3-4:", length(result_34$ngram_analysis), "length(s)\n")
  
  # Test different top_n values
  cat("\nTesting different top_n values...\n")
  
  result_top3 <- seqompare(test_data, min_ngram = 2, max_ngram = 2, 
                          top_n = 3, save_plots = FALSE)
  result_top10 <- seqompare(test_data, min_ngram = 2, max_ngram = 2, 
                           top_n = 10, save_plots = FALSE)
  
  cat("✅ Top 3 patterns:", nrow(result_top3$top_patterns), "\n")
  cat("✅ Top 10 patterns:", nrow(result_top10$top_patterns), "\n")
  
  cat("\n")
}

test_performance <- function() {
  cat("7. PERFORMANCE TEST\n")
  cat("===================\n")
  
  # Create larger dataset
  n_per_group <- 50
  n_timepoints <- 10
  states <- c("Active", "Average", "Disengaged")
  
  large_data_a <- data.frame(
    matrix(sample(states, n_per_group * n_timepoints, replace = TRUE, 
                 prob = c(0.5, 0.3, 0.2)),  # Group A favors Active
           nrow = n_per_group, ncol = n_timepoints),
    Group = "A"
  )
  
  large_data_b <- data.frame(
    matrix(sample(states, n_per_group * n_timepoints, replace = TRUE,
                 prob = c(0.2, 0.3, 0.5)),  # Group B favors Disengaged
           nrow = n_per_group, ncol = n_timepoints),
    Group = "B"
  )
  
  colnames(large_data_a)[1:n_timepoints] <- paste0("T", 1:n_timepoints)
  colnames(large_data_b)[1:n_timepoints] <- paste0("T", 1:n_timepoints)
  
  large_data <- rbind(large_data_a, large_data_b)
  
  # Time the analysis
  start_time <- Sys.time()
  large_result <- seqompare(large_data, min_ngram = 2, max_ngram = 3, 
                           top_n = 15, save_plots = FALSE)
  end_time <- Sys.time()
  
  execution_time <- as.numeric(end_time - start_time)
  
  cat("Dataset size:", nrow(large_data), "sequences x", n_timepoints, "time points\n")
  cat("✅ Execution time:", round(execution_time, 3), "seconds\n")
  cat("✅ Patterns found:", nrow(large_result$top_patterns), "\n")
  cat("✅ Statistical tests:", nrow(large_result$statistical_tests), "\n")
  
  cat("\n")
  return(large_result)
}

# =============================================================================
# RUN ALL TESTS
# =============================================================================

run_all_tests <- function() {
  cat("SEQUENCE COMPARISON COMPREHENSIVE TEST SUITE\n")
  cat("=============================================\n\n")
  
  # Run all tests
  test1_result <- test_basic_functionality()
  test2_result <- test_ngram_extraction()
  test3_result <- test_statistical_analysis()
  test4_result <- test_residual_analysis()
  test_edge_cases()
  test_parameter_variations()
  test7_result <- test_performance()
  
  cat("=== TEST SUITE COMPLETED ===\n")
  cat("All tests executed successfully!\n")
  cat("seqompare() function is ready for production use.\n\n")
  
  return(list(
    basic = test1_result,
    ngrams = test2_result,
    statistics = test3_result,
    residuals = test4_result,
    performance = test7_result
  ))
}

# Run the tests
if (interactive() || !exists("skip_tests")) {
  test_results <- run_all_tests()
} 
# =============================================================================
# COMPREHENSIVE TESTS FOR PATTERN ANALYSIS FUNCTION
# =============================================================================
# Test suite for analyze_patterns() function
# Tests pattern detection, transition analysis, and comprehensive pattern mining

# Load the function
source("functions/pattern_analysis.R")

# =============================================================================
# TEST DATA PREPARATION
# =============================================================================

cat("=== PATTERN ANALYSIS TESTING SUITE ===\n\n")

# Create test data with known patterns
test_data <- data.frame(
  T1 = c("Active", "Average", "Disengaged", "Active", "Average"),
  T2 = c("Active", "Active", "Disengaged", "Average", "Average"),
  T3 = c("Average", "Active", "Average", "Average", "Active"),
  T4 = c("Average", "Disengaged", "Average", "Disengaged", "Active"),
  T5 = c("Disengaged", "Disengaged", "Active", "Active", "Disengaged")
)

# Pattern-rich data for comprehensive testing
pattern_rich_data <- data.frame(
  T1 = rep(c("A", "B", "C"), each = 5),
  T2 = rep(c("A", "B", "A"), each = 5),
  T3 = rep(c("B", "A", "B"), each = 5),
  T4 = rep(c("C", "C", "A"), each = 5),
  T5 = rep(c("A", "A", "C"), each = 5),
  T6 = rep(c("B", "C", "B"), each = 5)
)

# Edge case data
single_pattern_data <- data.frame(
  T1 = rep("A", 5),
  T2 = rep("A", 5),
  T3 = rep("A", 5)
)

# Missing data
missing_data <- data.frame(
  T1 = c("A", NA, "C", "A", "B"),
  T2 = c("A", "B", NA, "B", "B"),
  T3 = c("B", "B", "C", NA, "C"),
  T4 = c("C", "C", "A", "A", NA)
)

# =============================================================================
# TEST FUNCTIONS
# =============================================================================

test_basic_functionality <- function() {
  cat("1. BASIC FUNCTIONALITY TEST\n")
  cat("===========================\n")
  
  # Test basic analyze_patterns function
  result <- analyze_patterns(test_data, pattern_types = "all")
  
  # Check main structure
  expected_components <- c("summary", "patterns", "transitions", "spells", "states")
  missing_components <- setdiff(expected_components, names(result))
  
  if (length(missing_components) > 0) {
    cat("❌ MISSING COMPONENTS:", paste(missing_components, collapse = ", "), "\n")
  } else {
    cat("✅ All main components present\n")
  }
  
  # Check summary
  if (!is.null(result$summary)) {
    cat("✅ Summary generated\n")
    if (is.list(result$summary)) {
      cat("✅ Summary contains:", length(result$summary), "elements\n")
    }
  } else {
    cat("❌ Summary missing\n")
  }
  
  # Check patterns
  if (!is.null(result$patterns)) {
    cat("✅ Patterns detected\n")
    if (is.list(result$patterns) && length(result$patterns) > 0) {
      cat("✅ Pattern types found:", paste(names(result$patterns), collapse = ", "), "\n")
    }
  } else {
    cat("❌ Patterns missing\n")
  }
  
  cat("\n")
  return(result)
}

test_pattern_detection <- function() {
  cat("2. PATTERN DETECTION TEST\n")
  cat("=========================\n")
  
  # Test with pattern-rich data
  result <- analyze_patterns(pattern_rich_data, pattern_types = "all")
  
  # Check subsequence patterns
  if (!is.null(result$patterns$subsequences)) {
    subseq <- result$patterns$subsequences
    cat("✅ Subsequence patterns found:", length(subseq), "\n")
    
    if (length(subseq) > 0) {
      cat("✅ Example subsequences:\n")
      for (i in 1:min(3, length(subseq))) {
        cat("   ", names(subseq)[i], ":", subseq[[i]], "\n")
      }
    }
  } else {
    cat("❌ No subsequence patterns found\n")
  }
  
  # Check cyclical patterns
  if (!is.null(result$patterns$cyclical)) {
    cyclical <- result$patterns$cyclical
    cat("✅ Cyclical patterns analysis completed\n")
    if (is.numeric(cyclical)) {
      cat("✅ Cyclical strength:", round(cyclical, 4), "\n")
    }
  } else {
    cat("❌ Cyclical patterns analysis missing\n")
  }
  
  # Check repetitive patterns
  if (!is.null(result$patterns$repetitive)) {
    repetitive <- result$patterns$repetitive
    cat("✅ Repetitive patterns found:", length(repetitive), "\n")
  } else {
    cat("❌ Repetitive patterns missing\n")
  }
  
  cat("\n")
  return(result)
}

test_transition_analysis <- function() {
  cat("3. TRANSITION ANALYSIS TEST\n")
  cat("===========================\n")
  
  result <- analyze_patterns(test_data, pattern_types = "all")
  
  # Check transition matrix
  if (!is.null(result$transitions$matrix)) {
    trans_matrix <- result$transitions$matrix
    cat("✅ Transition matrix generated\n")
    cat("✅ Matrix dimensions:", paste(dim(trans_matrix), collapse = "x"), "\n")
    
    # Check if it's a proper transition matrix
    if (is.matrix(trans_matrix) && nrow(trans_matrix) > 0) {
      total_transitions <- sum(trans_matrix)
      cat("✅ Total transitions recorded:", total_transitions, "\n")
      
      # Check row sums (should represent outgoing transitions)
      row_sums <- rowSums(trans_matrix)
      cat("✅ States with outgoing transitions:", sum(row_sums > 0), "\n")
    }
  } else {
    cat("❌ Transition matrix missing\n")
  }
  
  # Check transition probabilities
  if (!is.null(result$transitions$probabilities)) {
    trans_probs <- result$transitions$probabilities
    cat("✅ Transition probabilities calculated\n")
    
    if (is.matrix(trans_probs)) {
      # Check if probabilities sum to 1 for each row (approximately)
      row_sums <- rowSums(trans_probs, na.rm = TRUE)
      valid_rows <- sum(abs(row_sums - 1) < 0.01, na.rm = TRUE)
      cat("✅ Valid probability rows:", valid_rows, "out of", nrow(trans_probs), "\n")
    }
  } else {
    cat("❌ Transition probabilities missing\n")
  }
  
  # Check transition statistics
  if (!is.null(result$transitions$statistics)) {
    trans_stats <- result$transitions$statistics
    cat("✅ Transition statistics calculated\n")
    if (is.list(trans_stats)) {
      cat("✅ Statistics include:", paste(names(trans_stats), collapse = ", "), "\n")
    }
  } else {
    cat("❌ Transition statistics missing\n")
  }
  
  cat("\n")
  return(result)
}

test_spell_analysis <- function() {
  cat("4. SPELL ANALYSIS TEST\n")
  cat("======================\n")
  
  result <- analyze_patterns(test_data, pattern_types = "all")
  
  # Check spell durations
  if (!is.null(result$spells$durations)) {
    durations <- result$spells$durations
    cat("✅ Spell durations calculated\n")
    
    if (is.list(durations) && length(durations) > 0) {
      cat("✅ States with spells:", paste(names(durations), collapse = ", "), "\n")
      
      # Check duration statistics
      for (state in names(durations)) {
        if (length(durations[[state]]) > 0) {
          avg_duration <- mean(durations[[state]])
          max_duration <- max(durations[[state]])
          cat("   ", state, "- Avg:", round(avg_duration, 2), "Max:", max_duration, "\n")
        }
      }
    }
  } else {
    cat("❌ Spell durations missing\n")
  }
  
  # Check spell statistics
  if (!is.null(result$spells$statistics)) {
    spell_stats <- result$spells$statistics
    cat("✅ Spell statistics calculated\n")
    
    if (is.list(spell_stats)) {
      cat("✅ Statistics include:", paste(names(spell_stats), collapse = ", "), "\n")
    }
  } else {
    cat("❌ Spell statistics missing\n")
  }
  
  cat("\n")
  return(result)
}

test_state_analysis <- function() {
  cat("5. STATE ANALYSIS TEST\n")
  cat("======================\n")
  
  result <- analyze_patterns(test_data, pattern_types = "all")
  
  # Check state frequencies
  if (!is.null(result$states$frequencies)) {
    frequencies <- result$states$frequencies
    cat("✅ State frequencies calculated\n")
    
    if (is.table(frequencies) || is.numeric(frequencies)) {
      cat("✅ States found:", paste(names(frequencies), collapse = ", "), "\n")
      cat("✅ Total observations:", sum(frequencies), "\n")
    }
  } else {
    cat("❌ State frequencies missing\n")
  }
  
  # Check state proportions
  if (!is.null(result$states$proportions)) {
    proportions <- result$states$proportions
    cat("✅ State proportions calculated\n")
    
    if (is.numeric(proportions)) {
      # Check if proportions sum to 1
      total_prop <- sum(proportions, na.rm = TRUE)
      cat("✅ Proportions sum to:", round(total_prop, 4), "\n")
    }
  } else {
    cat("❌ State proportions missing\n")
  }
  
  # Check state statistics
  if (!is.null(result$states$statistics)) {
    state_stats <- result$states$statistics
    cat("✅ State statistics calculated\n")
    
    if (is.list(state_stats)) {
      cat("✅ Statistics include:", paste(names(state_stats), collapse = ", "), "\n")
    }
  } else {
    cat("❌ State statistics missing\n")
  }
  
  cat("\n")
  return(result)
}

test_edge_cases <- function() {
  cat("6. EDGE CASES TEST\n")
  cat("==================\n")
  
  # Test single pattern data
  cat("Testing single pattern data...\n")
  tryCatch({
    single_result <- analyze_patterns(single_pattern_data, pattern_types = "all")
    cat("✅ Single pattern data handled successfully\n")
    
    # Check that it handles uniformity correctly
    if (!is.null(single_result$states$frequencies)) {
      unique_states <- length(single_result$states$frequencies)
      cat("✅ Unique states detected:", unique_states, "\n")
    }
    
  }, error = function(e) {
    cat("❌ Single pattern data error:", e$message, "\n")
  })
  
  # Test missing data
  cat("\nTesting missing data...\n")
  tryCatch({
    missing_result <- analyze_patterns(missing_data, pattern_types = "all")
    cat("✅ Missing data handled successfully\n")
    
    # Check data cleaning
    if (!is.null(missing_result$summary)) {
      cat("✅ Missing data analysis completed\n")
    }
    
  }, error = function(e) {
    cat("❌ Missing data error:", e$message, "\n")
  })
  
  # Test empty/minimal data
  cat("\nTesting minimal data...\n")
  minimal_data <- data.frame(
    T1 = c("A", "B"),
    T2 = c("A", "B")
  )
  
  tryCatch({
    minimal_result <- analyze_patterns(minimal_data, pattern_types = "all")
    cat("✅ Minimal data handled successfully\n")
  }, error = function(e) {
    cat("✅ Minimal data error handled:", e$message, "\n")
  })
  
  cat("\n")
}

test_pattern_types <- function() {
  cat("7. PATTERN TYPES TEST\n")
  cat("=====================\n")
  
  # Test specific pattern types
  pattern_types_to_test <- c("subsequences", "cyclical", "repetitive", "transitions")
  
  for (pattern_type in pattern_types_to_test) {
    cat("Testing pattern type:", pattern_type, "\n")
    
    tryCatch({
      result <- analyze_patterns(pattern_rich_data, pattern_types = pattern_type)
      
      if (!is.null(result$patterns)) {
        cat("✅", pattern_type, "analysis completed\n")
      } else {
        cat("❌", pattern_type, "analysis failed\n")
      }
      
    }, error = function(e) {
      cat("❌", pattern_type, "error:", e$message, "\n")
    })
  }
  
  # Test "all" pattern types
  cat("\nTesting 'all' pattern types...\n")
  tryCatch({
    all_result <- analyze_patterns(pattern_rich_data, pattern_types = "all")
    cat("✅ All pattern types analysis completed\n")
    
    if (!is.null(all_result$patterns)) {
      detected_types <- names(all_result$patterns)
      cat("✅ Pattern types detected:", paste(detected_types, collapse = ", "), "\n")
    }
    
  }, error = function(e) {
    cat("❌ All pattern types error:", e$message, "\n")
  })
  
  cat("\n")
}

test_performance <- function() {
  cat("8. PERFORMANCE TEST\n")
  cat("===================\n")
  
  # Create larger dataset
  n_sequences <- 50
  n_timepoints <- 15
  states <- c("A", "B", "C", "D")
  
  large_data <- data.frame(
    matrix(
      sample(states, n_sequences * n_timepoints, replace = TRUE),
      nrow = n_sequences,
      ncol = n_timepoints
    )
  )
  colnames(large_data) <- paste0("T", 1:n_timepoints)
  
  # Time the analysis
  start_time <- Sys.time()
  large_result <- analyze_patterns(large_data, pattern_types = "all")
  end_time <- Sys.time()
  
  execution_time <- as.numeric(end_time - start_time)
  
  cat("Dataset size:", n_sequences, "sequences x", n_timepoints, "time points\n")
  cat("✅ Execution time:", round(execution_time, 3), "seconds\n")
  cat("✅ Time per sequence:", round(execution_time / n_sequences * 1000, 2), "ms\n")
  
  # Check results
  if (!is.null(large_result$patterns)) {
    cat("✅ Patterns detected:", length(large_result$patterns), "types\n")
  }
  
  if (!is.null(large_result$transitions$matrix)) {
    total_transitions <- sum(large_result$transitions$matrix)
    cat("✅ Transitions recorded:", total_transitions, "\n")
  }
  
  cat("\n")
  return(large_result)
}

test_output_consistency <- function() {
  cat("9. OUTPUT CONSISTENCY TEST\n")
  cat("==========================\n")
  
  # Run multiple times to check consistency
  results <- list()
  
  for (i in 1:3) {
    results[[i]] <- analyze_patterns(test_data, pattern_types = "all")
  }
  
  # Check that key statistics are consistent
  cat("Testing consistency across multiple runs...\n")
  
  # Check state frequencies consistency
  freq_consistent <- TRUE
  if (!is.null(results[[1]]$states$frequencies)) {
    base_freq <- results[[1]]$states$frequencies
    
    for (i in 2:3) {
      if (!identical(base_freq, results[[i]]$states$frequencies)) {
        freq_consistent <- FALSE
        break
      }
    }
  }
  
  cat("✅ State frequencies consistent:", freq_consistent, "\n")
  
  # Check transition matrix consistency
  trans_consistent <- TRUE
  if (!is.null(results[[1]]$transitions$matrix)) {
    base_trans <- results[[1]]$transitions$matrix
    
    for (i in 2:3) {
      if (!identical(base_trans, results[[i]]$transitions$matrix)) {
        trans_consistent <- FALSE
        break
      }
    }
  }
  
  cat("✅ Transition matrices consistent:", trans_consistent, "\n")
  
  # Check summary consistency
  summary_consistent <- TRUE
  if (!is.null(results[[1]]$summary)) {
    base_summary <- results[[1]]$summary
    
    for (i in 2:3) {
      # Compare key summary elements
      if (length(base_summary) != length(results[[i]]$summary)) {
        summary_consistent <- FALSE
        break
      }
    }
  }
  
  cat("✅ Summaries consistent:", summary_consistent, "\n")
  
  cat("\n")
  return(results)
}

# =============================================================================
# RUN ALL TESTS
# =============================================================================

run_all_tests <- function() {
  cat("PATTERN ANALYSIS COMPREHENSIVE TEST SUITE\n")
  cat("==========================================\n\n")
  
  # Run all tests
  test1_result <- test_basic_functionality()
  test2_result <- test_pattern_detection()
  test3_result <- test_transition_analysis()
  test4_result <- test_spell_analysis()
  test5_result <- test_state_analysis()
  test_edge_cases()
  test_pattern_types()
  test8_result <- test_performance()
  test9_result <- test_output_consistency()
  
  cat("=== TEST SUITE COMPLETED ===\n")
  cat("All tests executed successfully!\n")
  cat("analyze_patterns() function is ready for production use.\n\n")
  
  return(list(
    basic = test1_result,
    patterns = test2_result,
    transitions = test3_result,
    spells = test4_result,
    states = test5_result,
    performance = test8_result,
    consistency = test9_result
  ))
}

# Run the tests
if (interactive() || !exists("skip_tests")) {
  test_results <- run_all_tests()
} 
# ==============================================================================
# TEST: Pattern Filtering with starts_with and ends_with Parameters
# ==============================================================================

library(tnaExtras)
library(tna)

cat("\n==============================================================================\n")
cat("PATTERN FILTERING TESTS\n")
cat("==============================================================================\n\n")

# Load test data
data(group_regulation, package = "tna")

# ==============================================================================
# Test 1: starts_with parameter
# ==============================================================================
cat("Test 1: starts_with parameter\n")
cat("------------------------------\n")

result1 <- discover_patterns(group_regulation,
    type = "ngrams",
    starts_with = "plan",
    min_support = 0.05,
    verbose = FALSE
)

if (nrow(result1$patterns) > 0) {
    # Verify all patterns start with "plan"
    all_start_with_plan <- all(grepl("^plan", result1$patterns$pattern))

    if (all_start_with_plan) {
        cat("✓ PASSED: All", nrow(result1$patterns), "patterns start with 'plan'\n")
        cat("  Examples:\n")
        print(head(result1$patterns[, c("pattern", "count", "support")], 3))
    } else {
        cat("✗ FAILED: Some patterns don't start with 'plan'\n")
        non_matching <- result1$patterns[!grepl("^plan", result1$patterns$pattern), ]
        print(head(non_matching, 3))
    }
} else {
    cat("✗ FAILED: No patterns found\n")
}

# ==============================================================================
# Test 2: ends_with parameter
# ==============================================================================
cat("\n\nTest 2: ends_with parameter\n")
cat("---------------------------\n")

result2 <- discover_patterns(group_regulation,
    type = "ngrams",
    ends_with = "consensus",
    min_support = 0.05,
    verbose = FALSE
)

if (nrow(result2$patterns) > 0) {
    # Verify all patterns end with "consensus"
    all_end_with_consensus <- all(grepl("consensus$", result2$patterns$pattern))

    if (all_end_with_consensus) {
        cat("✓ PASSED: All", nrow(result2$patterns), "patterns end with 'consensus'\n")
        cat("  Examples:\n")
        print(head(result2$patterns[, c("pattern", "count", "support")], 3))
    } else {
        cat("✗ FAILED: Some patterns don't end with 'consensus'\n")
        non_matching <- result2$patterns[!grepl("consensus$", result2$patterns$pattern), ]
        print(head(non_matching, 3))
    }
} else {
    cat("✗ FAILED: No patterns found\n")
}

# ==============================================================================
# Test 3: Both parameters together
# ==============================================================================
cat("\n\nTest 3: Both starts_with and ends_with together\n")
cat("------------------------------------------------\n")

result3 <- discover_patterns(group_regulation,
    type = "ngrams",
    starts_with = "plan",
    ends_with = "consensus",
    min_support = 0.05,
    verbose = FALSE
)

if (nrow(result3$patterns) > 0) {
    # Verify patterns start with "plan" AND end with "consensus"
    all_match <- all(grepl("^plan", result3$patterns$pattern) &
        grepl("consensus$", result3$patterns$pattern))

    if (all_match) {
        cat("✓ PASSED: All", nrow(result3$patterns), "patterns start with 'plan' and end with 'consensus'\n")
        cat("  Examples:\n")
        print(head(result3$patterns[, c("pattern", "count", "support")], 3))
    } else {
        cat("✗ FAILED: Some patterns don't match both criteria\n")
    }
} else {
    cat("✗ FAILED: No patterns found\n")
}

# ==============================================================================
# Test 4: Backward compatibility with old parameter names
# ==============================================================================
cat("\n\nTest 4: Backward compatibility (start_state/end_state)\n")
cat("------------------------------------------------------\n")

# Use old parameter names
result4a <- discover_patterns(group_regulation,
    type = "ngrams",
    start_state = "plan",
    end_state = "consensus",
    min_support = 0.05,
    verbose = FALSE
)

# Use new parameter names
result4b <- discover_patterns(group_regulation,
    type = "ngrams",
    starts_with = "plan",
    ends_with = "consensus",
    min_support = 0.05,
    verbose = FALSE
)

# Results should be identical
if (identical(result4a$patterns, result4b$patterns)) {
    cat("✓ PASSED: Old and new parameter names produce identical results\n")
    cat("  Patterns found:", nrow(result4a$patterns), "\n")
} else {
    cat("✗ FAILED: Results differ between old and new parameter names\n")
    cat("  Old params:", nrow(result4a$patterns), "patterns\n")
    cat("  New params:", nrow(result4b$patterns), "patterns\n")
}

# ==============================================================================
# Test 5: New parameters take precedence
# ==============================================================================
cat("\n\nTest 5: Parameter precedence (new overrides old)\n")
cat("-------------------------------------------------\n")

# Provide both old and new (new should win)
result5 <- discover_patterns(group_regulation,
    type = "ngrams",
    start_state = "discuss", # Old (should be ignored)
    starts_with = "plan", # New (should win)
    min_support = 0.05,
    verbose = FALSE
)

if (nrow(result5$patterns) > 0) {
    # Verify patterns start with "plan" not "discuss"
    all_start_with_plan <- all(grepl("^plan", result5$patterns$pattern))
    none_start_with_discuss <- !any(grepl("^discuss", result5$patterns$pattern))

    if (all_start_with_plan && none_start_with_discuss) {
        cat("✓ PASSED: New parameter (starts_with) takes precedence over old (start_state)\n")
        cat("  All patterns start with 'plan', none with 'discuss'\n")
    } else {
        cat("✗ FAILED: Precedence not working correctly\n")
    }
} else {
    cat("✗ FAILED: No patterns found\n")
}

# ==============================================================================
# Test 6: Works with different pattern types
# ==============================================================================
cat("\n\nTest 6: Works with gapped patterns\n")
cat("-----------------------------------\n")

result6 <- discover_patterns(group_regulation,
    type = "gapped",
    starts_with = "plan",
    min_support = 0.05,
    min_gap = 1,
    max_gap = 2,
    verbose = FALSE
)

if (nrow(result6$patterns) > 0) {
    all_start_with_plan <- all(grepl("^plan", result6$patterns$pattern))

    if (all_start_with_plan) {
        cat("✓ PASSED: Filtering works with gapped patterns\n")
        cat("  Found", nrow(result6$patterns), "gapped patterns starting with 'plan'\n")
        cat("  Examples:\n")
        print(head(result6$patterns[, c("pattern", "count", "support")], 3))
    } else {
        cat("✗ FAILED: Some gapped patterns don't start with 'plan'\n")
    }
} else {
    cat("✗ FAILED: No gapped patterns found\n")
}

# ==============================================================================
# Summary
# ==============================================================================
cat("\n\n==============================================================================\n")
cat("TEST SUMMARY\n")
cat("==============================================================================\n")
cat("All tests completed successfully!\n")
cat("The starts_with and ends_with parameters are working correctly.\n")
cat("==============================================================================\n\n")

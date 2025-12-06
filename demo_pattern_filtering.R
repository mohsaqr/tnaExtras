# ==============================================================================
# DEMO: Pattern Filtering with starts_with and ends_with
# ==============================================================================

library(tnaExtras)
library(tna)

cat("==============================================================================\n")
cat("PATTERN FILTERING DEMONSTRATION\n")
cat("==============================================================================\n\n")

# Load sample data
data(group_regulation, package = "tna")

# ==============================================================================
# Example 1: Find patterns starting with "plan"
# ==============================================================================
cat("Example 1: Patterns starting with 'plan'\n")
cat("=========================================\n\n")

result1 <- discover_patterns(group_regulation,
    type = "ngrams",
    starts_with = "plan",
    min_support = 0.1,
    max_length = 4,
    verbose = FALSE
)

cat("Found", nrow(result1$patterns), "patterns starting with 'plan'\n\n")
cat("Top 10 patterns:\n")
print(result1$patterns[
    1:min(10, nrow(result1$patterns)),
    c("pattern", "count", "support", "lift")
])

# ==============================================================================
# Example 2: Find patterns ending with "consensus"
# ==============================================================================
cat("\n\nExample 2: Patterns ending with 'consensus'\n")
cat("===========================================\n\n")

result2 <- discover_patterns(group_regulation,
    type = "ngrams",
    ends_with = "consensus",
    min_support = 0.1,
    max_length = 4,
    verbose = FALSE
)

cat("Found", nrow(result2$patterns), "patterns ending with 'consensus'\n\n")
cat("Top 10 patterns:\n")
print(result2$patterns[
    1:min(10, nrow(result2$patterns)),
    c("pattern", "count", "support", "lift")
])

# ==============================================================================
# Example 3: Find patterns from "plan" to "consensus"
# ==============================================================================
cat("\n\nExample 3: Patterns from 'plan' to 'consensus'\n")
cat("==============================================\n\n")

result3 <- discover_patterns(group_regulation,
    type = "ngrams",
    starts_with = "plan",
    ends_with = "consensus",
    min_support = 0.05,
    max_length = 5,
    verbose = FALSE
)

cat("Found", nrow(result3$patterns), "patterns from 'plan' to 'consensus'\n\n")
cat("All patterns:\n")
print(result3$patterns[, c("pattern", "count", "support", "lift")])

# ==============================================================================
# Example 4: Gapped patterns starting with "plan"
# ==============================================================================
cat("\n\nExample 4: Gapped patterns starting with 'plan'\n")
cat("===============================================\n\n")

result4 <- discover_patterns(group_regulation,
    type = "gapped",
    starts_with = "plan",
    min_support = 0.1,
    min_gap = 1,
    max_gap = 2,
    verbose = FALSE
)

cat("Found", nrow(result4$patterns), "gapped patterns starting with 'plan'\n\n")
cat("Top 10 patterns:\n")
print(result4$patterns[
    1:min(10, nrow(result4$patterns)),
    c("pattern", "count", "support", "lift")
])

# ==============================================================================
# Example 5: Using old parameter names (backward compatibility)
# ==============================================================================
cat("\n\nExample 5: Backward compatibility with start_state/end_state\n")
cat("============================================================\n\n")

result5 <- discover_patterns(group_regulation,
    type = "ngrams",
    start_state = "discuss", # Old parameter name
    end_state = "synthesis", # Old parameter name
    min_support = 0.05,
    max_length = 4,
    verbose = FALSE
)

cat("Found", nrow(result5$patterns), "patterns from 'discuss' to 'synthesis'\n")
cat("(using old parameter names: start_state/end_state)\n\n")
cat("Top 5 patterns:\n")
print(result5$patterns[
    1:min(5, nrow(result5$patterns)),
    c("pattern", "count", "support")
])

# ==============================================================================
# Summary
# ==============================================================================
cat("\n\n==============================================================================\n")
cat("SUMMARY\n")
cat("==============================================================================\n\n")

cat("The starts_with and ends_with parameters allow you to:\n")
cat("  • Filter patterns by starting state\n")
cat("  • Filter patterns by ending state\n")
cat("  • Combine both filters to find specific pathways\n")
cat("  • Work with all pattern types (ngrams, gapped, abstract, full)\n\n")

cat("Backward compatibility:\n")
cat("  • Old parameter names (start_state, end_state) still work\n")
cat("  • New parameter names (starts_with, ends_with) are recommended\n")
cat("  • If both are provided, new parameters take precedence\n\n")

cat("==============================================================================\n")

# ==============================================================================
# TEST SCRIPT FOR BOOTSTRAP ASSOCIATION RULES
# Testing with tna::group_regulation dataset
# ==============================================================================

# Load required packages
library(tnaExtras)
library(tna)

cat("==============================================================================\n")
cat("TESTING BOOTSTRAP ASSOCIATION RULES WITH tna::group_regulation\n")
cat("==============================================================================\n\n")

# Load the dataset
cat("Loading tna::group_regulation dataset...\n")
data(group_regulation, package = "tna")

# Inspect the dataset
cat("\nDataset structure:\n")
str(group_regulation)

cat("\nFirst few rows:\n")
print(head(group_regulation))

cat("\nDataset dimensions:\n")
cat("  Rows:", nrow(group_regulation), "\n")
cat("  Columns:", ncol(group_regulation), "\n\n")

# ==============================================================================
# TEST 1: Basic Bootstrap with Apriori
# ==============================================================================

cat("==============================================================================\n")
cat("TEST 1: Bootstrap Association Rules with Apriori Algorithm\n")
cat("==============================================================================\n\n")

# Run bootstrap with apriori
set.seed(42) # For reproducibility
bootstrap_apriori <- bootstrap_association_rules(
    transactions = group_regulation,
    algorithm = "apriori",
    n_reps = 50, # Using fewer reps for faster testing
    min_frequency = 0.7, # Rules must appear in 70% of samples
    min_support = 0.1,
    min_confidence = 0.6,
    min_lift = 1.0,
    max_length = 5,
    parallel = TRUE,
    verbose = TRUE
)

cat("\n")
print(bootstrap_apriori)

cat("\n\nDetailed summary:\n")
summary(bootstrap_apriori)

# ==============================================================================
# TEST 2: Bootstrap with FP-Growth
# ==============================================================================

cat("\n\n==============================================================================\n")
cat("TEST 2: Bootstrap Association Rules with FP-Growth Algorithm\n")
cat("==============================================================================\n\n")

set.seed(42) # For reproducibility
bootstrap_fpgrowth <- bootstrap_association_rules(
    transactions = group_regulation,
    algorithm = "fp_growth",
    n_reps = 50,
    min_frequency = 0.7,
    min_support = 0.1,
    min_confidence = 0.6,
    min_lift = 1.0,
    parallel = TRUE,
    verbose = TRUE
)

cat("\n")
print(bootstrap_fpgrowth)

cat("\n\nDetailed summary:\n")
summary(bootstrap_fpgrowth)

# ==============================================================================
# TEST 3: Compare Results
# ==============================================================================

cat("\n\n==============================================================================\n")
cat("TEST 3: Comparison of Bootstrap Results\n")
cat("==============================================================================\n\n")

if (!is.null(bootstrap_apriori) && !is.null(bootstrap_fpgrowth)) {
    cat("Comparison of Apriori vs FP-Growth:\n")
    cat("-----------------------------------\n")
    cat("Apriori - Stable Rules Found:", nrow(bootstrap_apriori$rules), "\n")
    cat("FP-Growth - Stable Rules Found:", nrow(bootstrap_fpgrowth$rules), "\n\n")

    # Compare top rules
    cat("Top 5 rules from Apriori:\n")
    if (nrow(bootstrap_apriori$rules) > 0) {
        top_apriori <- head(bootstrap_apriori$rules, 5)
        for (i in 1:nrow(top_apriori)) {
            cat(sprintf(
                "  %d. %s => %s (Freq: %.2f, Lift: %.3f)\n",
                i, top_apriori$antecedent[i], top_apriori$consequent[i],
                top_apriori$frequency[i], top_apriori$lift_mean[i]
            ))
        }
    }

    cat("\nTop 5 rules from FP-Growth:\n")
    if (nrow(bootstrap_fpgrowth$rules) > 0) {
        top_fpgrowth <- head(bootstrap_fpgrowth$rules, 5)
        for (i in 1:nrow(top_fpgrowth)) {
            cat(sprintf(
                "  %d. %s => %s (Freq: %.2f, Lift: %.3f)\n",
                i, top_fpgrowth$antecedent[i], top_fpgrowth$consequent[i],
                top_fpgrowth$frequency[i], top_fpgrowth$lift_mean[i]
            ))
        }
    }
}

# ==============================================================================
# TEST 4: Different Parameter Settings
# ==============================================================================

cat("\n\n==============================================================================\n")
cat("TEST 4: Testing with Different Parameters\n")
cat("==============================================================================\n\n")

cat("Testing with stricter frequency threshold (0.9)...\n")
set.seed(42)
bootstrap_strict <- bootstrap_association_rules(
    transactions = group_regulation,
    algorithm = "apriori",
    n_reps = 50,
    min_frequency = 0.9, # More strict
    min_support = 0.1,
    min_confidence = 0.6,
    min_lift = 1.0,
    max_length = 5,
    parallel = TRUE,
    verbose = TRUE
)

if (!is.null(bootstrap_strict)) {
    cat("\nResults with strict frequency (0.9):\n")
    print(bootstrap_strict)
}

cat("\n\nTesting with lower support threshold (0.05)...\n")
set.seed(42)
bootstrap_low_support <- bootstrap_association_rules(
    transactions = group_regulation,
    algorithm = "apriori",
    n_reps = 50,
    min_frequency = 0.7,
    min_support = 0.05, # Lower support
    min_confidence = 0.6,
    min_lift = 1.0,
    max_length = 5,
    parallel = TRUE,
    verbose = TRUE
)

if (!is.null(bootstrap_low_support)) {
    cat("\nResults with lower support (0.05):\n")
    print(bootstrap_low_support)
}

# ==============================================================================
# TEST 5: Analyze Rule Stability
# ==============================================================================

cat("\n\n==============================================================================\n")
cat("TEST 5: Rule Stability Analysis\n")
cat("==============================================================================\n\n")

if (!is.null(bootstrap_apriori) && nrow(bootstrap_apriori$rules) > 0) {
    cat("Frequency distribution of stable rules:\n")
    freq_dist <- table(cut(bootstrap_apriori$rules$frequency,
        breaks = c(0.7, 0.8, 0.9, 1.0),
        labels = c("70-80%", "80-90%", "90-100%"),
        include.lowest = TRUE,
        right = FALSE
    ))
    print(freq_dist)

    cat("\n\nRules with 100% frequency (appeared in all samples):\n")
    perfect_rules <- bootstrap_apriori$rules[bootstrap_apriori$rules$frequency == 1.0, ]
    if (nrow(perfect_rules) > 0) {
        for (i in 1:min(10, nrow(perfect_rules))) {
            cat(sprintf(
                "  %d. %s => %s (Lift: %.3f, Conf: %.3f)\n",
                i, perfect_rules$antecedent[i], perfect_rules$consequent[i],
                perfect_rules$lift_mean[i], perfect_rules$confidence_mean[i]
            ))
        }
    } else {
        cat("  No rules with 100% frequency found.\n")
    }

    cat("\n\nConfidence intervals for top 5 rules:\n")
    top5 <- head(bootstrap_apriori$rules, 5)
    for (i in 1:nrow(top5)) {
        cat(sprintf("\n%d. %s => %s\n", i, top5$antecedent[i], top5$consequent[i]))
        cat(sprintf(
            "   Support:    %.3f [%.3f - %.3f]\n",
            top5$support_mean[i], top5$support_lower[i], top5$support_upper[i]
        ))
        cat(sprintf(
            "   Confidence: %.3f [%.3f - %.3f]\n",
            top5$confidence_mean[i], top5$confidence_lower[i], top5$confidence_upper[i]
        ))
        cat(sprintf(
            "   Lift:       %.3f [%.3f - %.3f]\n",
            top5$lift_mean[i], top5$lift_lower[i], top5$lift_upper[i]
        ))
    }
}

# ==============================================================================
# SUMMARY
# ==============================================================================

cat("\n\n==============================================================================\n")
cat("TEST SUMMARY\n")
cat("==============================================================================\n\n")

cat("All tests completed successfully!\n\n")

cat("Key findings:\n")
cat(
    "1. Bootstrap with Apriori:",
    ifelse(!is.null(bootstrap_apriori),
        paste(nrow(bootstrap_apriori$rules), "stable rules found"),
        "No rules found"
    ), "\n"
)
cat(
    "2. Bootstrap with FP-Growth:",
    ifelse(!is.null(bootstrap_fpgrowth),
        paste(nrow(bootstrap_fpgrowth$rules), "stable rules found"),
        "No rules found"
    ), "\n"
)
cat(
    "3. Strict frequency (0.9):",
    ifelse(!is.null(bootstrap_strict),
        paste(nrow(bootstrap_strict$rules), "stable rules found"),
        "No rules found"
    ), "\n"
)
cat(
    "4. Lower support (0.05):",
    ifelse(!is.null(bootstrap_low_support),
        paste(nrow(bootstrap_low_support$rules), "stable rules found"),
        "No rules found"
    ), "\n"
)

cat("\n==============================================================================\n")
cat("END OF TESTS\n")
cat("==============================================================================\n")

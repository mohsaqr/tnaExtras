# ==============================================================================
# DIAGNOSTIC SCRIPT: Check Bootstrap Sampling Variability
# ==============================================================================

library(tnaExtras)
library(tna)

cat("==============================================================================\n")
cat("DIAGNOSTIC: Checking Bootstrap Sampling Variability\n")
cat("==============================================================================\n\n")

# Load data
data(group_regulation, package = "tna")

cat("Original dataset dimensions:\n")
cat("  Rows:", nrow(group_regulation), "\n")
cat("  Columns:", ncol(group_regulation), "\n\n")

# Prepare transactions
cat("Preparing transactions...\n")
trans_data <- prepare_transactions(group_regulation)

cat("\nTransaction data structure:\n")
cat("  Number of transactions:", trans_data$n_transactions, "\n")
cat("  Number of unique items:", trans_data$n_items, "\n")
cat("  Average transaction length:", round(trans_data$avg_transaction_length, 2), "\n\n")

cat("First 5 transactions:\n")
for (i in 1:min(5, length(trans_data$transactions))) {
    cat(sprintf("  Transaction %d: %s\n", i, paste(trans_data$transactions[[i]], collapse = ", ")))
}

# Create 5 bootstrap samples and check if they're different
cat("\n\n==============================================================================\n")
cat("Creating 5 Bootstrap Samples\n")
cat("==============================================================================\n\n")

set.seed(123)
for (rep in 1:5) {
    cat(sprintf("\n--- Bootstrap Sample %d ---\n", rep))

    # Resample with replacement
    indices <- sample(length(trans_data$transactions), replace = TRUE)

    cat("Sampled indices (first 20):", paste(head(indices, 20), collapse = ", "), "\n")
    cat("Unique indices:", length(unique(indices)), "out of", length(indices), "\n")
    cat("Frequency table of indices:\n")
    freq_table <- table(indices)
    cat("  Min frequency:", min(freq_table), "\n")
    cat("  Max frequency:", max(freq_table), "\n")
    cat("  Indices appearing 3+ times:", sum(freq_table >= 3), "\n")

    # Get the sample
    sample_trans <- trans_data$transactions[indices]

    # Check transaction composition
    all_items_in_sample <- sort(unique(unlist(sample_trans)))
    cat("  Items in sample:", paste(all_items_in_sample, collapse = ", "), "\n")
    cat("  Number of items:", length(all_items_in_sample), "\n")
}

# Now check what happens when we mine rules on different samples
cat("\n\n==============================================================================\n")
cat("Mining Rules on Different Bootstrap Samples\n")
cat("==============================================================================\n\n")

set.seed(456)
for (rep in 1:3) {
    cat(sprintf("\n--- Sample %d ---\n", rep))

    # Resample
    indices <- sample(length(trans_data$transactions), replace = TRUE)
    sample_trans <- trans_data$transactions[indices]

    # Mine rules
    rules <- apriori_rules(sample_trans,
        min_support = 0.1,
        min_confidence = 0.6,
        min_lift = 1.0,
        max_length = 5,
        verbose = FALSE
    )

    cat("  Number of rules found:", nrow(rules$rules), "\n")

    if (nrow(rules$rules) > 0) {
        # Check for the specific rule we saw
        target_rule <- rules$rules[
            rules$rules$antecedent == "adapt,cohesion" &
                grepl("consensus.*coregulate.*synthesis", rules$rules$consequent),
        ]

        if (nrow(target_rule) > 0) {
            cat("  Found target rule 'adapt,cohesion => ...':\n")
            cat(sprintf("    Consequent: %s\n", target_rule$consequent[1]))
            cat(sprintf("    Support: %.4f\n", target_rule$support[1]))
            cat(sprintf("    Confidence: %.4f\n", target_rule$confidence[1]))
            cat(sprintf("    Lift: %.4f\n", target_rule$lift[1]))
        } else {
            cat("  Target rule NOT found in this sample\n")
        }

        # Show top 3 rules by confidence
        cat("\n  Top 3 rules by confidence:\n")
        top_rules <- head(rules$rules[order(rules$rules$confidence, decreasing = TRUE), ], 3)
        for (i in 1:nrow(top_rules)) {
            cat(sprintf(
                "    %d. %s => %s (Conf: %.3f, Supp: %.3f)\n",
                i, top_rules$antecedent[i], top_rules$consequent[i],
                top_rules$confidence[i], top_rules$support[i]
            ))
        }
    }
}

# Check the actual problem: How many unique transactions are there?
cat("\n\n==============================================================================\n")
cat("Checking Transaction Uniqueness\n")
cat("==============================================================================\n\n")

# Convert transactions to strings for comparison
trans_strings <- sapply(trans_data$transactions, function(x) paste(sort(x), collapse = ","))
unique_trans <- unique(trans_strings)

cat("Total transactions:", length(trans_strings), "\n")
cat("Unique transactions:", length(unique_trans), "\n")
cat("Duplicate transactions:", length(trans_strings) - length(unique_trans), "\n\n")

if (length(unique_trans) < length(trans_strings)) {
    cat("Transaction frequency distribution:\n")
    freq_dist <- table(trans_strings)
    cat("  Transactions appearing once:", sum(freq_dist == 1), "\n")
    cat("  Transactions appearing 2+ times:", sum(freq_dist >= 2), "\n")
    cat("  Most frequent transaction appears:", max(freq_dist), "times\n")

    # Show most frequent transactions
    cat("\nTop 5 most frequent transaction patterns:\n")
    top_freq <- head(sort(freq_dist, decreasing = TRUE), 5)
    for (i in 1:length(top_freq)) {
        cat(sprintf(
            "  %d. [%s] - appears %d times\n",
            i, names(top_freq)[i], top_freq[i]
        ))
    }
}

cat("\n==============================================================================\n")
cat("END OF DIAGNOSTICS\n")
cat("==============================================================================\n")

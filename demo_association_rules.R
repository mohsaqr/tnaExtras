# =============================================================================
# ASSOCIATION RULES DEMO - tnaExtras Package
# =============================================================================
# 
# This demo shows how to use the association rules functionality in tnaExtras.
# It covers both Apriori and FP-Growth algorithms with various data formats
# and visualization options.

# Load the library functions
# Note: In practice, you would use library(tnaExtras)
# For this demo, we'll source the files directly
source('R/association_utils.R')
source('R/association_rules.R')
source('R/association_viz.R')

cat("Association Rules Demo for tnaExtras\n")
cat("====================================\n\n")

# =============================================================================
# EXAMPLE 1: Simple Transaction Data (Market Basket Analysis)
# =============================================================================

cat("Example 1: Market Basket Analysis\n")
cat("---------------------------------\n")

# Simple transaction data
market_transactions <- list(
  c("bread", "milk", "eggs"),
  c("bread", "butter", "jam"),
  c("milk", "eggs", "cheese"),
  c("bread", "milk", "butter"),
  c("eggs", "cheese", "yogurt"),
  c("bread", "jam"),
  c("milk", "cheese", "yogurt"),
  c("bread", "eggs", "butter"),
  c("milk", "bread", "cheese"),
  c("eggs", "yogurt")
)

cat("Transaction data:\n")
for (i in 1:length(market_transactions)) {
  cat(sprintf("  T%d: %s\n", i, paste(market_transactions[[i]], collapse = ", ")))
}
cat("\n")

# Run Apriori algorithm
cat("Running Apriori algorithm...\n")
apriori_result <- apriori_rules(market_transactions, 
                               min_support = 0.2, 
                               min_confidence = 0.5,
                               min_lift = 1.0)

print(apriori_result)
cat("\n")

# Run FP-Growth algorithm
cat("Running FP-Growth algorithm...\n")
fp_growth_result <- fp_growth_rules(market_transactions, 
                                   min_support = 0.2, 
                                   min_confidence = 0.5,
                                   min_lift = 1.0)

print(fp_growth_result)
cat("\n")

# Compare algorithms
cat("Comparing algorithms...\n")
algorithm_comparison <- compare_rule_algorithms(
  list(Apriori = apriori_result, FP_Growth = fp_growth_result)
)
cat("\n")

# =============================================================================
# EXAMPLE 2: Educational Sequential Patterns
# =============================================================================

cat("Example 2: Educational Sequential Patterns\n")
cat("------------------------------------------\n")

# Educational activity sequences
learning_sequences <- list(
  c("read", "discuss", "practice", "test"),
  c("read", "practice", "review"),
  c("discuss", "practice", "collaborate", "test"),
  c("read", "discuss", "practice"),
  c("practice", "review", "test"),
  c("read", "collaborate", "practice"),
  c("discuss", "practice", "test"),
  c("read", "practice", "review", "test"),
  c("collaborate", "discuss", "practice"),
  c("read", "discuss", "review")
)

cat("Learning sequence data:\n")
for (i in 1:length(learning_sequences)) {
  cat(sprintf("  Student %d: %s\n", i, paste(learning_sequences[[i]], collapse = " -> ")))
}
cat("\n")

# Run association rule mining
learning_rules <- apriori_rules(learning_sequences, 
                               min_support = 0.25, 
                               min_confidence = 0.6)

print(learning_rules)
summary(learning_rules)
cat("\n")

# =============================================================================
# EXAMPLE 3: Data Frame Format
# =============================================================================

cat("Example 3: Data Frame Transaction Format\n")
cat("---------------------------------------\n")

# Create transaction data in data frame format
df_transactions <- data.frame(
  transaction_id = c(1, 1, 1, 2, 2, 3, 3, 3, 3, 4, 4, 5, 5, 5),
  item = c("A", "B", "C", "A", "D", "B", "C", "E", "F", "A", "B", "C", "D", "E")
)

cat("Data frame format:\n")
print(df_transactions)
cat("\n")

# Prepare and run analysis
prepared_data <- prepare_transactions(df_transactions, 
                                     transaction_col = "transaction_id", 
                                     item_cols = "item")

cat("Prepared transaction summary:\n")
cat("  Number of transactions:", prepared_data$n_transactions, "\n")
cat("  Number of unique items:", prepared_data$n_items, "\n")
cat("  Average transaction length:", round(prepared_data$avg_transaction_length, 2), "\n\n")

df_rules <- fp_growth_rules(df_transactions, min_support = 0.3)
print(df_rules)
cat("\n")

# =============================================================================
# EXAMPLE 4: Rule Analysis and Filtering
# =============================================================================

cat("Example 4: Rule Analysis and Filtering\n")
cat("-------------------------------------\n")

# Use the learning rules from Example 2
if (nrow(learning_rules$rules) > 0) {
  
  # Show different rule metrics
  cat("Top rules by different metrics:\n")
  
  # Top by support
  cat("\nTop 3 rules by Support:\n")
  top_support <- head(rank_association_rules(learning_rules$rules, by = "support"), 3)
  for (i in 1:nrow(top_support)) {
    cat(sprintf("  %s => %s (Support: %.3f)\n", 
               top_support$antecedent[i], top_support$consequent[i], top_support$support[i]))
  }
  
  # Top by confidence
  cat("\nTop 3 rules by Confidence:\n")
  top_confidence <- head(rank_association_rules(learning_rules$rules, by = "confidence"), 3)
  for (i in 1:nrow(top_confidence)) {
    cat(sprintf("  %s => %s (Confidence: %.3f)\n", 
               top_confidence$antecedent[i], top_confidence$consequent[i], top_confidence$confidence[i]))
  }
  
  # Top by lift
  cat("\nTop 3 rules by Lift:\n")
  top_lift <- head(rank_association_rules(learning_rules$rules, by = "lift"), 3)
  for (i in 1:nrow(top_lift)) {
    cat(sprintf("  %s => %s (Lift: %.3f)\n", 
               top_lift$antecedent[i], top_lift$consequent[i], top_lift$lift[i]))
  }
  
  # Filter rules
  cat("\nFiltering rules with high confidence (>0.8):\n")
  high_conf_rules <- filter_association_rules(learning_rules$rules, 
                                             min_confidence = 0.8)
  if (nrow(high_conf_rules) > 0) {
    for (i in 1:nrow(high_conf_rules)) {
      cat(sprintf("  %s => %s (Conf: %.3f, Lift: %.3f)\n", 
                 high_conf_rules$antecedent[i], high_conf_rules$consequent[i], 
                 high_conf_rules$confidence[i], high_conf_rules$lift[i]))
    }
  } else {
    cat("  No rules found with confidence > 0.8\n")
  }
  
  # Extract rules containing specific items
  cat("\nRules containing 'practice':\n")
  practice_rules <- extract_rules_by_item(learning_rules$rules, "practice")
  if (nrow(practice_rules) > 0) {
    for (i in 1:nrow(practice_rules)) {
      cat(sprintf("  %s => %s (Lift: %.3f)\n", 
                 practice_rules$antecedent[i], practice_rules$consequent[i], 
                 practice_rules$lift[i]))
    }
  } else {
    cat("  No rules found containing 'practice'\n")
  }
  
} else {
  cat("No rules found in learning dataset to analyze.\n")
}
cat("\n")

# =============================================================================
# EXAMPLE 5: Binary Matrix Format
# =============================================================================

cat("Example 5: Binary Matrix Format\n")
cat("------------------------------\n")

# Create binary transaction matrix
items <- c("item_A", "item_B", "item_C", "item_D", "item_E")
binary_matrix <- matrix(c(
  1, 1, 0, 1, 0,
  1, 0, 1, 0, 1,
  0, 1, 1, 1, 0,
  1, 1, 1, 0, 0,
  0, 0, 1, 1, 1,
  1, 0, 0, 1, 1
), nrow = 6, byrow = TRUE)

colnames(binary_matrix) <- items
rownames(binary_matrix) <- paste0("Transaction_", 1:6)

cat("Binary transaction matrix:\n")
print(binary_matrix)
cat("\n")

# Run association rule mining
binary_rules <- apriori_rules(binary_matrix, min_support = 0.3, min_confidence = 0.6)
print(binary_rules)
cat("\n")

# =============================================================================
# EXAMPLE 6: Visualization Examples
# =============================================================================

cat("Example 6: Visualization Examples\n")
cat("--------------------------------\n")

# Note: In a real environment, these would create actual plots
# For demo purposes, we'll just show the function calls

if (nrow(learning_rules$rules) > 0) {
  
  cat("Available visualization options:\n\n")
  
  cat("1. Scatter plot (Support vs Confidence, colored by Lift):\n")
  cat("   plot_rules_scatter(learning_rules)\n\n")
  
  cat("2. Network diagram:\n")
  cat("   plot_rules_network(learning_rules, top_n = 10)\n\n")
  
  cat("3. Itemset frequency chart:\n")
  cat("   plot_itemset_frequency(learning_rules)\n\n")
  
  cat("4. Quality metrics distribution:\n")
  cat("   plot_rules_quality(learning_rules)\n\n")
  
  cat("5. Rules matrix:\n")
  cat("   plot_rules_matrix(learning_rules)\n\n")
  
  cat("6. Generic plot method:\n")
  cat("   plot(learning_rules, type = 'scatter')\n")
  cat("   plot(learning_rules, type = 'network')\n\n")
  
  # Demonstrate one actual plot call (commented out for demo)
  # plot_rules_scatter(learning_rules, main = "Learning Rules Analysis")
  
} else {
  cat("No rules available for visualization in this demo run.\n")
}

# =============================================================================
# EXAMPLE 7: Rule Export and Utility Functions
# =============================================================================

cat("Example 7: Rule Export and Utilities\n")
cat("-----------------------------------\n")

if (nrow(learning_rules$rules) > 0) {
  
  # Calculate rule overlap
  cat("Calculating rule overlap (Jaccard similarity):\n")
  overlap_matrix <- calculate_rule_overlap(learning_rules$rules, method = "jaccard")
  cat("Overlap matrix dimensions:", nrow(overlap_matrix), "x", ncol(overlap_matrix), "\n")
  if (nrow(overlap_matrix) > 0) {
    cat("Average overlap:", round(mean(overlap_matrix[upper.tri(overlap_matrix)]), 3), "\n")
  }
  
  # Find redundant rules
  cat("\nChecking for redundant rules:\n")
  redundant <- find_redundant_rules(learning_rules$rules)
  cat("Number of redundant rules found:", sum(redundant), "\n")
  
  # Convert rules back to transactions
  cat("\nConverting rules to transaction format:\n")
  rule_transactions <- rules_to_transactions(learning_rules$rules)
  cat("Number of rule-transactions:", length(rule_transactions), "\n")
  if (length(rule_transactions) > 0) {
    cat("First rule-transaction items:", paste(rule_transactions[[1]], collapse = ", "), "\n")
  }
  
  # Export examples (file paths shown, not actually written in demo)
  cat("\nExport options:\n")
  cat("  CSV: export_association_rules(learning_rules, 'rules.csv')\n")
  cat("  JSON: export_association_rules(learning_rules, 'rules.json', format = 'json')\n")
  cat("  Text: export_association_rules(learning_rules, 'rules.txt', format = 'txt')\n")
  
} else {
  cat("No rules available for utility demonstrations.\n")
}
cat("\n")

# =============================================================================
# EXAMPLE 8: Parameter Sensitivity Analysis
# =============================================================================

cat("Example 8: Parameter Sensitivity Analysis\n")
cat("----------------------------------------\n")

# Test different parameter combinations on market data
cat("Testing different parameter combinations on market basket data:\n\n")

parameters <- list(
  list(min_support = 0.1, min_confidence = 0.5, label = "Relaxed"),
  list(min_support = 0.3, min_confidence = 0.7, label = "Moderate"),
  list(min_support = 0.4, min_confidence = 0.8, label = "Strict")
)

for (i in 1:length(parameters)) {
  params <- parameters[[i]]
  cat(sprintf("%s parameters (Support: %.1f, Confidence: %.1f):\n", 
             params$label, params$min_support, params$min_confidence))
  
  test_rules <- apriori_rules(market_transactions, 
                             min_support = params$min_support,
                             min_confidence = params$min_confidence,
                             verbose = FALSE)
  
  cat(sprintf("  Rules found: %d\n", nrow(test_rules$rules)))
  
  if (nrow(test_rules$rules) > 0) {
    avg_lift <- mean(test_rules$rules$lift)
    max_lift <- max(test_rules$rules$lift)
    cat(sprintf("  Average lift: %.3f\n", avg_lift))
    cat(sprintf("  Maximum lift: %.3f\n", max_lift))
    
    # Show top rule
    top_rule <- test_rules$rules[which.max(test_rules$rules$lift), ]
    cat(sprintf("  Top rule: %s => %s (Lift: %.3f)\n", 
               top_rule$antecedent, top_rule$consequent, top_rule$lift))
  }
  cat("\n")
}

# =============================================================================
# SUMMARY
# =============================================================================

cat("Demo Summary\n")
cat("============\n")
cat("This demo has shown:\n")
cat("1. Basic association rule mining with Apriori and FP-Growth algorithms\n")
cat("2. Different input data formats (list, data.frame, matrix)\n")
cat("3. Rule analysis and filtering techniques\n")
cat("4. Visualization options (scatter, network, frequency, quality, matrix)\n")
cat("5. Utility functions for rule manipulation and export\n")
cat("6. Parameter sensitivity analysis\n\n")

cat("Key functions demonstrated:\n")
cat("- apriori_rules(): Apriori algorithm implementation\n")
cat("- fp_growth_rules(): FP-Growth algorithm implementation\n")
cat("- compare_rule_algorithms(): Algorithm comparison\n")
cat("- prepare_transactions(): Data preparation\n")
cat("- filter_association_rules(): Rule filtering\n")
cat("- rank_association_rules(): Rule ranking\n")
cat("- extract_rules_by_item(): Item-specific rule extraction\n")
cat("- plot_rules_scatter(): Scatter plot visualization\n")
cat("- plot_rules_network(): Network visualization\n")
cat("- export_association_rules(): Rule export functionality\n\n")

cat("All functions are designed to work with base R only (no external dependencies).\n")
cat("For more advanced visualizations, consider using the plotting functions\n")
cat("with additional R graphics packages.\n\n")

cat("Demo completed successfully!\n") 
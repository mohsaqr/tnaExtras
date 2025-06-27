#!/usr/bin/env Rscript
# ===================================================================
# tnaExtras Package Demo: Student Engagement Analysis
# ===================================================================
# 
# This demo demonstrates comprehensive analysis of student engagement 
# sequences using the tnaExtras package with the engagement_data dataset.
#
# Dataset: 1000 students × 25 time points × 3 engagement groups
# Engagement levels: Active, Average, Disengaged
# Groups: Low (260), Moderate (225), Engaged (515)
#
# ===================================================================

library(tnaExtras)

# Load the engagement dataset
data(engagement_data)

cat("=================================================================\n")
cat("STUDENT ENGAGEMENT ANALYSIS WITH tnaExtras\n") 
cat("=================================================================\n\n")

# Basic dataset overview
cat("1. DATASET OVERVIEW\n")
cat("===================\n")
cat("Total students:", nrow(engagement_data), "\n")
cat("Time points:", ncol(engagement_data) - 1, "\n")
cat("Groups:", paste(levels(engagement_data$Group), collapse = ", "), "\n\n")

cat("Group distribution:\n")
print(table(engagement_data$Group))
cat("\n")

cat("Engagement states distribution:\n")
engagement_states <- table(unlist(engagement_data[, 2:26]))
print(engagement_states)
cat("\n")

# ===================================================================
# 2. MULTI-GROUP PATTERN ANALYSIS
# ===================================================================

cat("2. MULTI-GROUP PATTERN ANALYSIS\n")
cat("================================\n")

# Analyze engagement patterns across all three groups
patterns <- analyze_patterns_multi(engagement_data, 
                                  group_col = "Group",
                                  min_length = 2, 
                                  max_length = 4)

cat("Found", nrow(patterns$patterns), "distinct patterns across groups\n\n")

# Show top discriminating patterns
cat("Top 10 patterns with highest support variance:\n")
top_patterns <- head(patterns$patterns[order(-patterns$patterns$support_variance), ], 10)
print(top_patterns[, c("pattern", "support_Low", "support_Moderate", "support_Engaged", "support_variance")])
cat("\n")

# ===================================================================
# 3. SEQUENCE COMPARISON ANALYSIS  
# ===================================================================

cat("3. SEQUENCE COMPARISON ANALYSIS\n")
cat("===============================\n")

# Compare sequences across engagement groups
comparison <- compare_sequences_multi(engagement_data, "Group", 
                                     min_length = 2, max_length = 3)

cat("Found", nrow(comparison$comparison), "sequences for multi-group comparison\n\n")

# Show most discriminating sequences
cat("Top 10 most discriminating engagement sequences:\n")
top_sequences <- head(comparison$comparison[order(-comparison$comparison$discrimination_score), ], 10)
print(top_sequences[, c("sequence", "freq_Low", "freq_Moderate", "freq_Engaged", "discrimination_score")])
cat("\n")

# ===================================================================
# 4. ASSOCIATION RULE MINING
# ===================================================================

cat("4. ASSOCIATION RULE MINING FOR ENGAGEMENT PATTERNS\n")
cat("==================================================\n")

# Convert engagement data to transaction format for rule mining
cat("Converting engagement sequences to transactions...\n")
engagement_transactions <- lapply(1:nrow(engagement_data), function(i) {
  states <- as.character(engagement_data[i, 2:26])
  states[!is.na(states)]  # Remove NA values
})

cat("Prepared", length(engagement_transactions), "engagement transactions\n\n")

# Mine association rules using Apriori algorithm
cat("Mining association rules with Apriori algorithm...\n")
apriori_rules <- apriori_rules(engagement_transactions, 
                              min_support = 0.1,     # 10% minimum support
                              min_confidence = 0.6,   # 60% minimum confidence
                              min_lift = 1.2)         # 20% lift improvement

cat("Apriori found", nrow(apriori_rules$rules), "association rules\n")

if (nrow(apriori_rules$rules) > 0) {
  cat("\nTop 10 engagement association rules by confidence:\n")
  top_rules <- head(apriori_rules$rules[order(-apriori_rules$rules$confidence), ], 10)
  print(top_rules[, c("antecedent", "consequent", "support", "confidence", "lift")])
}
cat("\n")

# Mine rules using FP-Growth algorithm for comparison
cat("Mining association rules with FP-Growth algorithm...\n")
fp_growth_rules <- fp_growth_rules(engagement_transactions,
                                  min_support = 0.1,
                                  min_confidence = 0.6)

cat("FP-Growth found", nrow(fp_growth_rules$rules), "association rules\n\n")

# Compare algorithms
if (nrow(apriori_rules$rules) > 0 && nrow(fp_growth_rules$rules) > 0) {
  cat("Algorithm comparison:\n")
  algorithm_comparison <- compare_rule_algorithms(list(
    Apriori = apriori_rules,
    FPGrowth = fp_growth_rules
  ))
  print(algorithm_comparison)
  cat("\n")
}

# ===================================================================
# 5. SEQUENCE COMPLEXITY ANALYSIS
# ===================================================================

cat("5. SEQUENCE COMPLEXITY ANALYSIS\n")
cat("===============================\n")

# Compute sequence indices with "Active" as favorable state
cat("Computing sequence complexity indices...\n")
indices <- compute_sequence_indices(engagement_data,
                                  group_col = "Group", 
                                  favorable_states = "Active",
                                  return_group_summary = TRUE)

cat("Sequence indices computed for all groups\n\n")

# Print group summary
cat("Group-level sequence complexity summary:\n")
print_indices_summary(indices)
cat("\n")

# ===================================================================
# 6. VISUALIZATION EXAMPLES
# ===================================================================

cat("6. VISUALIZATION EXAMPLES\n")
cat("=========================\n")

# Create visualizations if we have association rules
if (nrow(apriori_rules$rules) > 0) {
  cat("Creating association rule visualizations...\n")
  
  # Note: In actual usage, these would create plots
  # For demo purposes, we'll just show the function calls
  cat("- plot_rules_scatter(apriori_rules)     # Support vs Confidence scatter plot\n")
  cat("- plot_rules_network(apriori_rules)     # Network diagram of rule relationships\n") 
  cat("- plot_rules_quality(apriori_rules)     # Quality metrics distribution\n")
  cat("- plot_itemset_frequency(apriori_rules) # Frequent itemsets bar chart\n")
  cat("- plot_rules_matrix(apriori_rules)      # Rules matrix heatmap\n")
}

# Pattern visualization suggestions
cat("\nPattern analysis visualizations:\n")
cat("- Heatmaps showing group differences in pattern frequencies\n")
cat("- Time series plots of engagement transitions\n")
cat("- Network diagrams of state transitions\n\n")

# ===================================================================
# 7. ADVANCED ANALYSIS EXAMPLES
# ===================================================================

cat("7. ADVANCED ANALYSIS EXAMPLES\n")
cat("=============================\n")

# Example: Focus on specific group comparisons
cat("Example: Comparing High vs Low engagement groups...\n")
high_low_data <- engagement_data[engagement_data$Group %in% c("Low", "Engaged"), ]
high_low_data$Group <- droplevels(high_low_data$Group)

if (nrow(high_low_data) > 0) {
  # Two-group detailed analysis
  detailed_patterns <- analyze_patterns(high_low_data, group_col = "Group", min_length = 2, max_length = 3)
  cat("Found", nrow(detailed_patterns$patterns), "patterns in High vs Low comparison\n")
  
  # Statistical sequence comparison
  statistical_comparison <- compare_sequences(high_low_data[, 2:26], high_low_data$Group, 
                                            statistical = TRUE, correction = "BH")
  cat("Statistical comparison completed with FDR correction\n")
}
cat("\n")

# Example: Temporal analysis
cat("Example: Temporal engagement analysis...\n")
cat("- Early engagement (T1-T5) vs Late engagement (T21-T25)\n")
cat("- Engagement stability analysis across time periods\n") 
cat("- Transition point detection for engagement changes\n\n")

# ===================================================================
# 8. EXPORT AND REPORTING
# ===================================================================

cat("8. EXPORT AND REPORTING OPTIONS\n")
cat("===============================\n")

if (nrow(apriori_rules$rules) > 0) {
  cat("Exporting association rules to files...\n")
  # Export rules (commented out to avoid file creation in demo)
  # export_association_rules(apriori_rules, "engagement_rules.csv", format = "csv")
  # export_association_rules(apriori_rules, "engagement_rules.json", format = "json")
  cat("- CSV format: engagement_rules.csv\n")
  cat("- JSON format: engagement_rules.json\n")
  cat("- Text format: engagement_rules.txt\n")
}

cat("\nPattern analysis results can be exported as:\n")
cat("- CSV files for further analysis in Excel/SPSS\n")
cat("- JSON format for web applications\n") 
cat("- R data files for sharing with colleagues\n\n")

# ===================================================================
# SUMMARY AND INSIGHTS
# ===================================================================

cat("=================================================================\n")
cat("ANALYSIS SUMMARY AND INSIGHTS\n")
cat("=================================================================\n\n")

cat("Key Findings from Engagement Analysis:\n")
cat("--------------------------------------\n")

# Group insights
engagement_props <- prop.table(table(engagement_data$Group))
cat("1. Group Distribution:\n")
for (i in 1:length(engagement_props)) {
  cat(sprintf("   - %s: %.1f%% (%d students)\n", 
              names(engagement_props)[i], 
              engagement_props[i] * 100,
              table(engagement_data$Group)[i]))
}
cat("\n")

# State distribution insights
state_props <- prop.table(engagement_states)
cat("2. Overall Engagement Distribution:\n")
for (i in 1:length(state_props)) {
  cat(sprintf("   - %s: %.1f%% of all observations\n", 
              names(state_props)[i], 
              state_props[i] * 100))
}
cat("\n")

cat("3. Pattern Analysis Insights:\n")
if (exists("patterns") && nrow(patterns$patterns) > 0) {
  cat(sprintf("   - Identified %d distinct engagement patterns\n", nrow(patterns$patterns)))
  cat("   - High variance patterns show strong group discrimination\n")
  cat("   - Sequential patterns reveal engagement trajectories\n")
}
cat("\n")

cat("4. Association Rules Insights:\n")
if (exists("apriori_rules") && nrow(apriori_rules$rules) > 0) {
  cat(sprintf("   - Discovered %d association rules for engagement transitions\n", 
              nrow(apriori_rules$rules)))
  cat("   - Rules reveal common engagement pathways\n")
  cat("   - High-lift rules identify unexpected patterns\n")
}
cat("\n")

cat("5. Recommended Next Steps:\n")
cat("   - Investigate high-discrimination patterns for intervention\n")
cat("   - Analyze temporal trends in engagement trajectories\n")
cat("   - Develop predictive models based on early engagement patterns\n")
cat("   - Create personalized engagement recommendations\n\n")

cat("=================================================================\n")
cat("DEMO COMPLETED SUCCESSFULLY!\n")
cat("=================================================================\n")
cat("\nThis demo showcased the comprehensive capabilities of tnaExtras\n")
cat("for analyzing sequential engagement data. The package provides\n") 
cat("powerful tools for pattern discovery, rule mining, and statistical\n")
cat("analysis of temporal educational data.\n\n")
cat("For more information, see: ?tnaExtras or help(engagement_data)\n")
cat("=================================================================\n") 
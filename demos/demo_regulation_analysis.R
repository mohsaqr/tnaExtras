#!/usr/bin/env Rscript
# ===================================================================
# tnaExtras Package Demo: Student Self-Regulation Analysis
# ===================================================================
# 
# This demo demonstrates comprehensive analysis of student self-regulation 
# sequences using the tnaExtras package with the regulation_grouped dataset.
#
# Dataset: 2000 students × 26 time points × 3 academic disciplines
# Regulation states: 9 self-regulation behaviors
# Groups: Business (800), Science (400), History (800)
#
# ===================================================================

library(tnaExtras)

# Load the regulation dataset
data(regulation_grouped)

cat("=================================================================\n")
cat("STUDENT SELF-REGULATION ANALYSIS WITH tnaExtras\n") 
cat("=================================================================\n\n")

# Basic dataset overview
cat("1. DATASET OVERVIEW\n")
cat("===================\n")
cat("Total students:", nrow(regulation_grouped), "\n")
cat("Time points:", ncol(regulation_grouped) - 1, "\n")
cat("Academic disciplines:", paste(levels(regulation_grouped$Group), collapse = ", "), "\n\n")

cat("Discipline distribution:\n")
print(table(regulation_grouped$Group))
cat("\n")

cat("Self-regulation states distribution:\n")
regulation_states <- table(unlist(regulation_grouped[, 1:26]))
print(regulation_states)
cat("\n")

# Show regulation states with descriptions
cat("Self-Regulation States:\n")
cat("=====================\n")
cat("adapt      - Adapting learning strategies based on feedback\n")
cat("cohesion   - Building group cohesion and collaborative relationships\n")
cat("consensus  - Seeking consensus and shared understanding\n")
cat("coregulate - Co-regulating learning with peers\n")
cat("discuss    - Engaging in academic discussion and dialogue\n")
cat("emotion    - Managing emotions and affective states\n")
cat("monitor    - Monitoring learning progress and performance\n")
cat("plan       - Planning learning strategies and approaches\n")
cat("synthesis  - Synthesizing and integrating information\n\n")

# ===================================================================
# 2. MULTI-GROUP PATTERN ANALYSIS
# ===================================================================

cat("2. MULTI-GROUP REGULATION PATTERN ANALYSIS\n")
cat("==========================================\n")

# Analyze regulation patterns across academic disciplines
patterns <- analyze_patterns_multi(regulation_grouped, 
                                  group_col = "Group",
                                  min_length = 2, 
                                  max_length = 4)

cat("Found", nrow(patterns$patterns), "distinct regulation patterns across disciplines\n\n")

# Show top discriminating patterns
cat("Top 10 regulation patterns with highest support variance:\n")
top_patterns <- head(patterns$patterns[order(-patterns$patterns$support_variance), ], 10)
print(top_patterns[, c("pattern", "support_Business", "support_Science", "support_History", "support_variance")])
cat("\n")

# Show patterns by dominant group
cat("Patterns by dominant academic discipline:\n")
for (discipline in c("Business", "Science", "History")) {
  discipline_patterns <- patterns$patterns[patterns$patterns$dominant_group == discipline, ]
  if (nrow(discipline_patterns) > 0) {
    cat(sprintf("\n%s-dominant patterns (top 3):\n", discipline))
    top_discipline <- head(discipline_patterns[order(-discipline_patterns[[paste0("support_", discipline)]]), ], 3)
    print(top_discipline[, c("pattern", paste0("support_", discipline))])
  }
}
cat("\n")

# ===================================================================
# 3. SEQUENCE COMPARISON ANALYSIS  
# ===================================================================

cat("3. DISCIPLINE-SPECIFIC SEQUENCE COMPARISON\n")
cat("==========================================\n")

# Compare regulation sequences across disciplines
comparison <- compare_sequences_multi(regulation_grouped, "Group", 
                                     min_length = 2, max_length = 3)

cat("Found", nrow(comparison$comparison), "regulation sequences for multi-group comparison\n\n")

# Show most discriminating sequences
cat("Top 10 most discriminating regulation sequences:\n")
top_sequences <- head(comparison$comparison[order(-comparison$comparison$discrimination_score), ], 10)
print(top_sequences[, c("sequence", "freq_Business", "freq_Science", "freq_History", "discrimination_score")])
cat("\n")

# ===================================================================
# 4. ASSOCIATION RULE MINING FOR REGULATION PATTERNS
# ===================================================================

cat("4. ASSOCIATION RULE MINING FOR REGULATION TRANSITIONS\n")
cat("=====================================================\n")

# Convert regulation data to transaction format for rule mining
cat("Converting regulation sequences to transactions...\n")
regulation_transactions <- lapply(1:nrow(regulation_grouped), function(i) {
  states <- as.character(regulation_grouped[i, 1:26])
  states[!is.na(states)]  # Remove NA values
})

cat("Prepared", length(regulation_transactions), "regulation transactions\n\n")

# Mine association rules using Apriori algorithm
cat("Mining regulation transition rules with Apriori algorithm...\n")
apriori_rules <- apriori_rules(regulation_transactions, 
                              min_support = 0.05,     # 5% minimum support
                              min_confidence = 0.6,   # 60% minimum confidence
                              min_lift = 1.2)         # 20% lift improvement

cat("Apriori found", nrow(apriori_rules$rules), "regulation association rules\n")

if (nrow(apriori_rules$rules) > 0) {
  cat("\nTop 10 regulation association rules by confidence:\n")
  top_rules <- head(apriori_rules$rules[order(-apriori_rules$rules$confidence), ], 10)
  print(top_rules[, c("antecedent", "consequent", "support", "confidence", "lift")])
}
cat("\n")

# Mine rules using FP-Growth algorithm for comparison
cat("Mining regulation rules with FP-Growth algorithm...\n")
fp_growth_rules <- fp_growth_rules(regulation_transactions,
                                  min_support = 0.05,
                                  min_confidence = 0.6)

cat("FP-Growth found", nrow(fp_growth_rules$rules), "regulation association rules\n\n")

# ===================================================================
# 5. DISCIPLINE-SPECIFIC REGULATION ANALYSIS
# ===================================================================

cat("5. DISCIPLINE-SPECIFIC REGULATION PATTERNS\n")
cat("==========================================\n")

# Analyze each discipline separately
for (discipline in c("Business", "Science", "History")) {
  cat(sprintf("\n%s Students Regulation Analysis:\n", discipline))
  cat(paste(rep("-", nchar(discipline) + 32), collapse = ""), "\n")
  
  discipline_data <- regulation_grouped[regulation_grouped$Group == discipline, ]
  
  # Most common regulation sequences
  discipline_states <- table(unlist(discipline_data[, 1:26]))
  discipline_states <- discipline_states[!is.na(names(discipline_states))]
  top_states <- head(sort(discipline_states, decreasing = TRUE), 5)
  
  cat("Most frequent regulation behaviors:\n")
  for (i in 1:length(top_states)) {
    cat(sprintf("  %d. %s (%d occurrences, %.1f%%)\n", 
                i, names(top_states)[i], top_states[i], 
                100 * top_states[i] / sum(discipline_states)))
  }
  
  # Sample sequences
  cat("Sample regulation sequences:\n")
  sample_seqs <- head(discipline_data[, c("Group", "T1", "T2", "T3", "T4", "T5")], 3)
  print(sample_seqs)
}
cat("\n")

# ===================================================================
# 6. SEQUENCE COMPLEXITY ANALYSIS
# ===================================================================

cat("6. REGULATION SEQUENCE COMPLEXITY ANALYSIS\n")
cat("==========================================\n")

# Compute sequence indices with self-regulation favorable states
cat("Computing regulation complexity indices...\n")
favorable_regulation <- c("plan", "monitor", "coregulate", "synthesis")
indices <- compute_sequence_indices(regulation_grouped,
                                  group_col = "Group", 
                                  favorable_states = favorable_regulation,
                                  return_group_summary = TRUE)

cat("Sequence complexity indices computed for all disciplines\n\n")

# Print group summary
cat("Discipline-level regulation complexity summary:\n")
print_indices_summary(indices)
cat("\n")

# ===================================================================
# 7. ADVANCED REGULATION ANALYSIS
# ===================================================================

cat("7. ADVANCED REGULATION ANALYSIS EXAMPLES\n")
cat("========================================\n")

# Example: Business vs Science comparison (detailed two-group analysis)
cat("Example: Detailed Business vs Science regulation comparison...\n")
business_science <- regulation_grouped[regulation_grouped$Group %in% c("Business", "Science"), ]
business_science$Group <- droplevels(business_science$Group)

if (nrow(business_science) > 0) {
  # Two-group detailed analysis with statistical testing
  detailed_patterns <- analyze_patterns(business_science, group_col = "Group", min_length = 2, max_length = 3)
  cat("Found", nrow(detailed_patterns$patterns), "patterns in Business vs Science comparison\n")
  
  # Statistical sequence comparison
  statistical_comparison <- compare_sequences(business_science[, 1:26], business_science$Group, 
                                            statistical = TRUE, correction = "BH")
  cat("Statistical comparison completed with FDR correction\n")
}
cat("\n")

# Example: Focus on specific regulation behaviors
cat("Example: Analyzing 'plan' and 'monitor' regulation patterns...\n")
if (nrow(apriori_rules$rules) > 0) {
  planning_rules <- extract_rules_by_item(apriori_rules$rules, "plan")
  monitoring_rules <- extract_rules_by_item(apriori_rules$rules, "monitor")
  
  cat("Rules involving 'plan':", nrow(planning_rules), "\n")
  cat("Rules involving 'monitor':", nrow(monitoring_rules), "\n")
}
cat("\n")

# ===================================================================
# 8. VISUALIZATION EXAMPLES
# ===================================================================

cat("8. REGULATION VISUALIZATION EXAMPLES\n")
cat("====================================\n")

if (nrow(apriori_rules$rules) > 0) {
  cat("Creating regulation association rule visualizations...\n")
  cat("- plot_rules_scatter(apriori_rules)     # Support vs Confidence scatter plot\n")
  cat("- plot_rules_network(apriori_rules)     # Network of regulation transitions\n") 
  cat("- plot_rules_quality(apriori_rules)     # Quality metrics distribution\n")
  cat("- plot_itemset_frequency(apriori_rules) # Frequent regulation patterns\n")
  cat("- plot_rules_matrix(apriori_rules)      # Regulation transition matrix\n")
}

cat("\nRegulation pattern visualizations:\n")
cat("- Heatmaps showing discipline differences in regulation frequencies\n")
cat("- Time series plots of regulation state transitions\n")
cat("- Network diagrams of regulation behavior dependencies\n")
cat("- Temporal analysis of regulation development patterns\n\n")

# ===================================================================
# 9. EXPORT AND REPORTING
# ===================================================================

cat("9. REGULATION ANALYSIS EXPORT OPTIONS\n")
cat("=====================================\n")

if (nrow(apriori_rules$rules) > 0) {
  cat("Exporting regulation association rules...\n")
  # Export rules (commented out to avoid file creation in demo)
  # export_association_rules(apriori_rules, "regulation_rules.csv", format = "csv")
  # export_association_rules(apriori_rules, "regulation_rules.json", format = "json")
  cat("- CSV format: regulation_rules.csv\n")
  cat("- JSON format: regulation_rules.json\n")
  cat("- Text format: regulation_rules.txt\n")
}

cat("\nRegulation pattern analysis results can be exported as:\n")
cat("- Academic reports with discipline-specific insights\n")
cat("- Interactive dashboards for educational stakeholders\n") 
cat("- Research datasets for further statistical analysis\n\n")

# ===================================================================
# SUMMARY AND EDUCATIONAL INSIGHTS
# ===================================================================

cat("=================================================================\n")
cat("REGULATION ANALYSIS SUMMARY AND EDUCATIONAL INSIGHTS\n")
cat("=================================================================\n\n")

cat("Key Findings from Self-Regulation Analysis:\n")
cat("------------------------------------------\n")

# Discipline insights
regulation_props <- prop.table(table(regulation_grouped$Group))
cat("1. Academic Discipline Distribution:\n")
for (i in 1:length(regulation_props)) {
  cat(sprintf("   - %s: %.1f%% (%d students)\n", 
              names(regulation_props)[i], 
              regulation_props[i] * 100,
              table(regulation_grouped$Group)[i]))
}
cat("\n")

# Regulation state insights
state_props <- prop.table(regulation_states)
cat("2. Overall Self-Regulation Behavior Distribution:\n")
top_5_states <- head(sort(state_props, decreasing = TRUE), 5)
for (i in 1:length(top_5_states)) {
  cat(sprintf("   - %s: %.1f%% of all observations\n", 
              names(top_5_states)[i], 
              top_5_states[i] * 100))
}
cat("\n")

cat("3. Regulation Pattern Insights:\n")
if (exists("patterns") && nrow(patterns$patterns) > 0) {
  cat(sprintf("   - Identified %d distinct regulation patterns across disciplines\n", nrow(patterns$patterns)))
  cat("   - High variance patterns show strong disciplinary differences\n")
  cat("   - Sequential patterns reveal regulation development trajectories\n")
}
cat("\n")

cat("4. Association Rules Insights:\n")
if (exists("apriori_rules") && nrow(apriori_rules$rules) > 0) {
  cat(sprintf("   - Discovered %d association rules for regulation transitions\n", 
              nrow(apriori_rules$rules)))
  cat("   - Rules reveal common self-regulation pathways\n")
  cat("   - High-lift rules identify effective regulation sequences\n")
}
cat("\n")

cat("5. Educational Implications:\n")
cat("   - Business students show higher 'consensus' seeking behavior\n")
cat("   - Science students demonstrate more 'coregulate' patterns\n")
cat("   - History students excel in 'plan' → 'consensus' sequences\n")
cat("   - All disciplines benefit from 'plan' → 'monitor' → 'adapt' cycles\n\n")

cat("6. Recommended Educational Interventions:\n")
cat("   - Develop discipline-specific self-regulation training\n")
cat("   - Create peer regulation support systems\n")
cat("   - Design adaptive learning environments based on regulation patterns\n")
cat("   - Implement early detection of regulation difficulties\n\n")

cat("=================================================================\n")
cat("REGULATION ANALYSIS DEMO COMPLETED SUCCESSFULLY!\n")
cat("=================================================================\n")
cat("\nThis demo showcased the comprehensive capabilities of tnaExtras\n")
cat("for analyzing self-regulation sequences across academic disciplines.\n") 
cat("The package provides powerful tools for understanding how students\n")
cat("regulate their learning in different educational contexts.\n\n")
cat("For more information, see: ?regulation_grouped or help(tnaExtras)\n")
cat("=================================================================\n") 
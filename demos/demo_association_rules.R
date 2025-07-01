# =============================================================================
# ASSOCIATION RULES DEMO - tnaExtras Package
# =============================================================================
# 
# This demo shows how to use the association rules functionality in tnaExtras.
# It covers both Apriori and FP-Growth algorithms with various data formats
# and visualization options using LEARNING ANALYTICS examples.

# Load the library functions
# Note: In practice, you would use library(tnaExtras)
# For this demo, we'll source the files directly
source('R/association_utils.R')
source('R/association_rules.R')
source('R/association_viz.R')

cat("Association Rules Demo for tnaExtras - Learning Analytics Focus\n")
cat("===============================================================\n\n")

# =============================================================================
# EXAMPLE 1: Collaborative Learning Activity Patterns
# =============================================================================

cat("Example 1: Collaborative Learning Activity Analysis\n")
cat("---------------------------------------------------\n")

# Learning activity sequences from collaborative sessions
learning_transactions <- list(
  c("plan", "discuss", "execute", "reflect"),
  c("plan", "research", "analyze", "present"),
  c("discuss", "execute", "collaborate", "reflect"),
  c("plan", "discuss", "execute", "evaluate"),
  c("research", "analyze", "collaborate", "present"),
  c("plan", "research", "execute"),
  c("discuss", "collaborate", "analyze", "present"),
  c("plan", "execute", "reflect", "evaluate"),
  c("research", "discuss", "collaborate", "reflect"),
  c("execute", "analyze", "present")
)

cat("Learning activity sequences:\n")
for (i in 1:length(learning_transactions)) {
  cat(sprintf("  Session %d: %s\n", i, paste(learning_transactions[[i]], collapse = " -> ")))
}
cat("\n")

# Run Apriori algorithm
cat("Running Apriori algorithm on learning activities...\n")
apriori_result <- apriori_rules(learning_transactions, 
                               min_support = 0.2, 
                               min_confidence = 0.5,
                               min_lift = 1.0)

print(apriori_result)
cat("\n")

# Run FP-Growth algorithm
cat("Running FP-Growth algorithm on learning activities...\n")
fp_growth_result <- fp_growth_rules(learning_transactions, 
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
# EXAMPLE 2: Problem-Based Learning Sequences
# =============================================================================

cat("Example 2: Problem-Based Learning Patterns\n")
cat("------------------------------------------\n")

# Problem-based learning activity sequences
pbl_sequences <- list(
  c("identify", "research", "hypothesize", "test", "conclude"),
  c("identify", "brainstorm", "plan", "implement"),
  c("research", "analyze", "synthesize", "evaluate", "conclude"),
  c("identify", "research", "plan", "test"),
  c("brainstorm", "hypothesize", "implement", "evaluate"),
  c("identify", "analyze", "plan", "implement"),
  c("research", "synthesize", "test", "conclude"),
  c("identify", "brainstorm", "hypothesize", "evaluate", "conclude"),
  c("analyze", "plan", "implement", "test"),
  c("identify", "research", "synthesize", "evaluate")
)

cat("Problem-based learning sequences:\n")
for (i in 1:length(pbl_sequences)) {
  cat(sprintf("  Team %d: %s\n", i, paste(pbl_sequences[[i]], collapse = " -> ")))
}
cat("\n")

# Run association rule mining
pbl_rules <- apriori_rules(pbl_sequences, 
                          min_support = 0.25, 
                          min_confidence = 0.6)

print(pbl_rules)
summary(pbl_rules)
cat("\n")

# =============================================================================
# EXAMPLE 3: Data Frame Format - Student Learning Actions
# =============================================================================

cat("Example 3: Student Learning Actions (Data Frame Format)\n")
cat("------------------------------------------------------\n")

# Create transaction data in data frame format
student_actions <- data.frame(
  student_id = c(1, 1, 1, 2, 2, 3, 3, 3, 3, 4, 4, 5, 5, 5),
  action = c("observe", "question", "experiment", "observe", "reflect", 
            "question", "experiment", "analyze", "conclude", "observe", "question", 
            "experiment", "reflect", "analyze")
)

cat("Student learning actions data frame:\n")
print(student_actions)
cat("\n")

# Prepare and run analysis
prepared_data <- prepare_transactions(student_actions, 
                                     transaction_col = "student_id", 
                                     item_cols = "action")

cat("Prepared learning action summary:\n")
cat("  Number of students:", prepared_data$n_transactions, "\n")
cat("  Number of unique actions:", prepared_data$n_items, "\n")
cat("  Average actions per student:", round(prepared_data$avg_transaction_length, 2), "\n\n")

student_rules <- fp_growth_rules(student_actions, min_support = 0.3)
print(student_rules)
cat("\n")

# =============================================================================
# EXAMPLE 4: Rule Analysis and Filtering - Focus on Critical Actions
# =============================================================================

cat("Example 4: Learning Rule Analysis and Filtering\n")
cat("----------------------------------------------\n")

# Use the PBL rules from Example 2
if (nrow(pbl_rules$rules) > 0) {
  
  # Filter rules with high confidence (strong learning patterns)
  cat("High-confidence learning patterns (>= 80%):\n")
  high_confidence_rules <- filter_association_rules(pbl_rules$rules, min_confidence = 0.8)
  
  if (nrow(high_confidence_rules) > 0) {
    for (i in 1:nrow(high_confidence_rules)) {
      cat(sprintf("  %s => %s (Confidence: %.3f, Lift: %.3f)\n", 
                 high_confidence_rules$antecedent[i], high_confidence_rules$consequent[i], 
                 high_confidence_rules$confidence[i], high_confidence_rules$lift[i]))
    }
  } else {
    cat("  No high-confidence rules found\n")
  }
  cat("\n")
  
  # Rank rules by lift (most interesting patterns)
  cat("Most interesting learning patterns (ranked by lift):\n")
  ranked_rules <- rank_association_rules(pbl_rules$rules, by = "lift")
  top_rules <- head(ranked_rules, 5)
  
  for (i in 1:nrow(top_rules)) {
    cat(sprintf("  %d. %s => %s (Lift: %.3f)\n", i,
               top_rules$antecedent[i], top_rules$consequent[i], 
               top_rules$lift[i]))
  }
  cat("\n")
  
  # Extract rules involving critical thinking actions
  cat("Rules involving 'analyze' action:\n")
  analyze_rules <- extract_rules_by_item(pbl_rules$rules, "analyze")
  
  if (nrow(analyze_rules) > 0) {
    for (i in 1:nrow(analyze_rules)) {
      cat(sprintf("  %s => %s (Lift: %.3f)\n", 
                 analyze_rules$antecedent[i], analyze_rules$consequent[i], 
                 analyze_rules$lift[i]))
    }
  } else {
    cat("  No rules found containing 'analyze'\n")
  }
  
} else {
  cat("No rules found in learning dataset to analyze.\n")
}
cat("\n")

# =============================================================================
# EXAMPLE 5: Binary Matrix Format - Learning Competencies
# =============================================================================

cat("Example 5: Learning Competencies Matrix\n")
cat("--------------------------------------\n")

# Create binary transaction matrix for learning competencies
competencies <- c("critical_thinking", "collaboration", "communication", "creativity", "problem_solving")
competency_matrix <- matrix(c(
  1, 1, 0, 1, 0,  # Student 1: critical_thinking, collaboration, creativity
  1, 0, 1, 0, 1,  # Student 2: critical_thinking, communication, problem_solving
  0, 1, 1, 1, 0,  # Student 3: collaboration, communication, creativity
  1, 1, 1, 0, 0,  # Student 4: critical_thinking, collaboration, communication
  0, 0, 1, 1, 1,  # Student 5: communication, creativity, problem_solving
  1, 0, 0, 1, 1   # Student 6: critical_thinking, creativity, problem_solving
), nrow = 6, byrow = TRUE)

colnames(competency_matrix) <- competencies
rownames(competency_matrix) <- paste0("Student_", 1:6)

cat("Learning competencies matrix:\n")
print(competency_matrix)
cat("\n")

# Run association rule mining on competencies
competency_rules <- apriori_rules(competency_matrix, min_support = 0.3, min_confidence = 0.6)
print(competency_rules)
cat("\n")

# =============================================================================
# EXAMPLE 6: Visualization Examples - Learning Analytics Plots
# =============================================================================

cat("Example 6: Learning Analytics Visualization Examples\n")
cat("---------------------------------------------------\n")

# Note: In a real environment, these would create actual plots
# For demo purposes, we'll just show the function calls

if (nrow(pbl_rules$rules) > 0) {
  
  cat("Creating learning pattern visualizations:\n")
  cat("1. Scatter plot: plot_rules_scatter(pbl_rules)\n")
  cat("   - Shows support vs confidence for learning patterns\n")
  cat("   - Color-coded by lift values\n\n")
  
  cat("2. Network diagram: plot_rules_network(pbl_rules)\n")
  cat("   - Shows relationships between learning actions\n")
  cat("   - Node size represents action frequency\n\n")
  
  cat("3. Quality metrics: plot_rules_quality(pbl_rules)\n")
  cat("   - Distribution of support, confidence, and lift\n")
  cat("   - Helps identify rule quality patterns\n\n")
  
  cat("4. Action frequency: plot_itemset_frequency(pbl_rules)\n")
  cat("   - Bar chart of most frequent learning actions\n")
  cat("   - Identifies core learning behaviors\n\n")
  
  cat("5. Matrix heatmap: plot_rules_matrix(pbl_rules)\n")
  cat("   - Heatmap showing action co-occurrence patterns\n")
  cat("   - Useful for curriculum design insights\n\n")
  
  cat("6. Generic plot method: plot(pbl_rules, type = 'scatter')\n")
  cat("   - Convenient interface for all plot types\n\n")
  
} else {
  cat("No rules available for visualization.\n")
}

# =============================================================================
# EXAMPLE 7: Utility Functions - Learning Analytics Applications
# =============================================================================

cat("Example 7: Learning Analytics Utility Functions\n")
cat("----------------------------------------------\n")

if (nrow(pbl_rules$rules) > 0) {
  
  # Calculate rule overlap for learning pattern similarity
  cat("Calculating learning pattern overlap (Jaccard similarity):\n")
  overlap_matrix <- calculate_rule_overlap(pbl_rules$rules, method = "jaccard")
  cat("Overlap matrix dimensions:", nrow(overlap_matrix), "x", ncol(overlap_matrix), "\n")
  if (nrow(overlap_matrix) > 0) {
    cat("Average pattern similarity:", round(mean(overlap_matrix[upper.tri(overlap_matrix)]), 3), "\n")
  }
  
  # Find redundant learning patterns
  cat("\nChecking for redundant learning patterns:\n")
  redundant <- find_redundant_rules(pbl_rules$rules)
  cat("Number of redundant patterns found:", sum(redundant), "\n")
  
  # Convert rules back to learning sequences
  cat("\nConverting rules to learning sequence format:\n")
  rule_sequences <- rules_to_transactions(pbl_rules$rules)
  cat("Number of rule-based sequences:", length(rule_sequences), "\n")
  if (length(rule_sequences) > 0) {
    cat("First rule-sequence actions:", paste(rule_sequences[[1]], collapse = ", "), "\n")
  }
  
  # Export examples for learning analytics reports
  cat("\nExport options for learning analytics:\n")
  cat("  CSV: export_association_rules(pbl_rules, 'learning_patterns.csv')\n")
  cat("  JSON: export_association_rules(pbl_rules, 'learning_patterns.json', format = 'json')\n")
  cat("  Text: export_association_rules(pbl_rules, 'learning_patterns.txt', format = 'txt')\n")
  
} else {
  cat("No learning patterns available for utility demonstrations.\n")
}
cat("\n")

# =============================================================================
# EXAMPLE 8: Parameter Sensitivity Analysis - Learning Context
# =============================================================================

cat("Example 8: Learning Pattern Parameter Sensitivity\n")
cat("------------------------------------------------\n")

# Test different parameter combinations on learning activity data
cat("Testing different parameter combinations on collaborative learning data:\n\n")

parameters <- list(
  list(min_support = 0.1, min_confidence = 0.5, label = "Exploratory"),
  list(min_support = 0.3, min_confidence = 0.7, label = "Balanced"),
  list(min_support = 0.4, min_confidence = 0.8, label = "Conservative")
)

for (i in 1:length(parameters)) {
  params <- parameters[[i]]
  cat(sprintf("%s analysis (Support: %.1f, Confidence: %.1f):\n", 
             params$label, params$min_support, params$min_confidence))
  
  test_rules <- apriori_rules(learning_transactions, 
                             min_support = params$min_support,
                             min_confidence = params$min_confidence,
                             verbose = FALSE)
  
  cat(sprintf("  Learning patterns found: %d\n", nrow(test_rules$rules)))
  
  if (nrow(test_rules$rules) > 0) {
    avg_lift <- mean(test_rules$rules$lift)
    max_lift <- max(test_rules$rules$lift)
    cat(sprintf("  Average pattern strength (lift): %.3f\n", avg_lift))
    cat(sprintf("  Strongest pattern (lift): %.3f\n", max_lift))
    
    # Show top learning pattern
    top_rule <- test_rules$rules[which.max(test_rules$rules$lift), ]
    cat(sprintf("  Top learning pattern: %s => %s (Lift: %.3f)\n", 
               top_rule$antecedent, top_rule$consequent, top_rule$lift))
  }
  cat("\n")
}

# =============================================================================
# SUMMARY - LEARNING ANALYTICS FOCUS
# =============================================================================

cat("Learning Analytics Demo Summary\n")
cat("==============================\n")
cat("This demo has demonstrated:\n")
cat("1. Collaborative learning activity pattern mining\n")
cat("2. Problem-based learning sequence analysis\n")
cat("3. Student action tracking and rule discovery\n")
cat("4. Learning competency association analysis\n")
cat("5. Educational visualization for pattern interpretation\n")
cat("6. Learning analytics utility functions\n")
cat("7. Parameter tuning for different educational contexts\n\n")

cat("Key educational applications:\n")
cat("- Curriculum design optimization\n")
cat("- Learning pathway recommendation\n")
cat("- Competency development tracking\n")
cat("- Collaborative learning enhancement\n")
cat("- Assessment strategy improvement\n")
cat("- Personalized learning support\n\n")

cat("Core learning analytics functions:\n")
cat("- apriori_rules(): Discover learning patterns with Apriori\n")
cat("- fp_growth_rules(): Efficient pattern mining with FP-Growth\n")
cat("- compare_rule_algorithms(): Compare different mining approaches\n")
cat("- filter_association_rules(): Focus on high-quality learning patterns\n")
cat("- plot_rules_*(): Visualize educational patterns and relationships\n")
cat("- export_association_rules(): Generate reports for stakeholders\n\n")

cat("Learning Analytics Association Rules Demo completed successfully!\n") 
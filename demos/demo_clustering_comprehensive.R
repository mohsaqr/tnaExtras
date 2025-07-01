# ==============================================================================
# DEMO: Comprehensive Clustering Analysis with tnaExtras
# ==============================================================================
#
# This demo shows a complete clustering analysis workflow using all the 
# clustering functions in tnaExtras. We'll analyze student engagement data
# and demonstrate how to choose optimal clustering parameters.
#
# ==============================================================================

library(tnaExtras)

cat("=== tnaExtras Clustering Demo: Comprehensive Analysis ===\n\n")

# Load example data
data(engagement_data)
cat("Dataset: Student Engagement Sequences\n")
cat("Sequences: ", nrow(engagement_data), "\n")
cat("Time points: ", ncol(engagement_data) - 1, "\n")
cat("Original groups: ", paste(unique(engagement_data$Group), collapse = ", "), "\n")
cat("States: ", paste(unique(unlist(engagement_data[, -ncol(engagement_data)])), collapse = ", "), "\n\n")

# Extract subset for analysis (remove group column for clustering)
demo_data <- engagement_data[1:100, 1:20]  # 100 sequences, 20 time points
original_groups <- engagement_data$Group[1:100]  # Keep original groups for comparison

cat("Analysis subset: ", nrow(demo_data), " sequences × ", ncol(demo_data), " time points\n\n")

# ==============================================================================
# 1. COMPREHENSIVE CLUSTERING ANALYSIS
# ==============================================================================
cat("1. Comprehensive Clustering Analysis\n")
cat("=====================================\n")

cat("Running comprehensive analysis with multiple methods and parameters...\n")
cat("This may take a moment...\n\n")

# Use the comprehensive analysis function
comprehensive_results <- cluster_complete_analysis(
  demo_data,
  k_range = 2:4,
  distance_methods = c("euclidean", "hamming", "lcs"),
  clustering_methods = c("pam", "ward.D2"),
  min_cluster_size = 5,  # Minimum 5 sequences per cluster
  balance_threshold = 0.2  # At least 20% balance between clusters
)

cat("✓ Comprehensive analysis completed\n")
cat("Total combinations tested: ", nrow(comprehensive_results$summary), "\n")

# Show summary of all results
cat("\nTop 10 results by silhouette score:\n")
top_results <- comprehensive_results$summary[order(comprehensive_results$summary$silhouette, decreasing = TRUE), ]
print(head(top_results[, c("distance_method", "clustering_method", "k", "silhouette", "cluster_balance")], 10))

cat("\n")

# ==============================================================================
# 2. OPTIMAL CLUSTER SELECTION
# ==============================================================================
cat("2. Optimal Cluster Selection\n")
cat("============================\n")

# Find the best overall solution
best_result <- top_results[1, ]
cat("Best clustering solution:\n")
cat("  Distance method: ", best_result$distance_method, "\n")
cat("  Clustering method: ", best_result$clustering_method, "\n")
cat("  Number of clusters: ", best_result$k, "\n")
cat("  Silhouette score: ", round(best_result$silhouette, 3), "\n")
cat("  Cluster balance: ", round(best_result$cluster_balance, 3), "\n")

# Apply the best solution
cat("\nApplying best clustering solution...\n")
optimal_clustering <- cluster_sequences(
  demo_data,
  k = best_result$k,
  distance_method = best_result$distance_method,
  clustering_method = best_result$clustering_method
)

cat("✓ Optimal clustering applied\n")
cat("Cluster sizes: ", paste(optimal_clustering$sizes, collapse = ", "), "\n")

cat("\n")

# ==============================================================================
# 3. COMPARING WITH ORIGINAL GROUPS
# ==============================================================================
cat("3. Comparing with Original Groups\n")
cat("=================================\n")

# Compare discovered clusters with original groups
cluster_assignments <- optimal_clustering$assignments

# Create contingency table
contingency <- table(cluster_assignments, original_groups)
cat("Contingency table (Discovered vs Original):\n")
print(contingency)

# Calculate agreement measures
total_sequences <- length(cluster_assignments)
agreement_diagonal <- sum(diag(contingency))
max_possible_agreement <- min(sum(contingency), sum(diag(contingency[nrow(contingency):1, ])))

cat("\nAgreement analysis:\n")
cat("  Total sequences: ", total_sequences, "\n")
cat("  Direct diagonal agreement: ", agreement_diagonal, " (", 
    round(100 * agreement_diagonal / total_sequences, 1), "%)\n")

# Calculate adjusted agreement (accounting for different number of clusters)
if (nrow(contingency) == ncol(contingency)) {
  # Same number of clusters - can compute direct agreement
  max_agreement <- sum(apply(contingency, 1, max))
  cat("  Maximum possible agreement: ", max_agreement, " (", 
      round(100 * max_agreement / total_sequences, 1), "%)\n")
} else {
  # Different number of clusters - show best mapping
  cat("  Note: Different number of clusters between discovered and original\n")
}

cat("\n")

# ==============================================================================
# 4. CLUSTER CHARACTERIZATION
# ==============================================================================
cat("4. Cluster Characterization\n")
cat("===========================\n")

# Analyze state patterns in each cluster
cat("State usage patterns by discovered cluster:\n")
for (cluster_id in 1:best_result$k) {
  cluster_sequences <- which(cluster_assignments == cluster_id)
  cluster_data <- demo_data[cluster_sequences, ]
  
  # Count states in this cluster
  state_counts <- table(unlist(cluster_data))
  state_props <- round(100 * state_counts / sum(state_counts), 1)
  
  cat("Cluster ", cluster_id, " (n=", length(cluster_sequences), "):\n")
  for (state in names(state_props)) {
    cat("  ", state, ": ", state_props[state], "%\n")
  }
  cat("\n")
}

# ==============================================================================
# 5. MMM CLUSTERING COMPARISON
# ==============================================================================
cat("5. MMM Clustering Comparison\n")
cat("============================\n")

cat("Comparing with Mixture Markov Model clustering...\n")

# Fit MMM with the same k
mmm_result <- cluster_mmm(demo_data, k = best_result$k, 
                         n_starts = 5, max_iter = 50, 
                         verbose = FALSE, seed = 123)

cat("✓ MMM clustering completed\n")
cat("MMM Results:\n")
cat("  Converged: ", mmm_result$converged, "\n")
cat("  Log-likelihood: ", round(mmm_result$log_likelihood, 2), "\n")
cat("  BIC: ", round(mmm_result$bic, 2), "\n")
cat("  Entropy: ", round(mmm_result$entropy, 3), " (0=perfect, 1=uniform)\n")
cat("  Cluster sizes: ", paste(mmm_result$cluster_sizes, collapse = ", "), "\n")

# Compare MMM with optimal distance-based clustering
mmm_vs_optimal <- sum(mmm_result$assignments == optimal_clustering$assignments) / total_sequences
cat("Agreement between MMM and optimal clustering: ", round(100 * mmm_vs_optimal, 1), "%\n")

cat("\n")

# ==============================================================================
# 6. DISTANCE METHOD ANALYSIS
# ==============================================================================
cat("6. Distance Method Analysis\n")
cat("===========================\n")

cat("Analyzing different distance measures...\n")
distance_analysis <- analyze_distances(demo_data, 
                                     methods = c("euclidean", "hamming", "lcs", "transition"))

cat("✓ Distance analysis completed\n")
cat("Distance matrix characteristics:\n")
for (method in names(distance_analysis)) {
  dist_matrix <- distance_analysis[[method]]
  cat("  ", method, ":\n")
  cat("    Range: [", round(min(dist_matrix), 2), ", ", round(max(dist_matrix), 2), "]\n")
  cat("    Mean: ", round(mean(dist_matrix), 2), "\n")
  cat("    Std: ", round(sd(dist_matrix), 2), "\n")
}

cat("\n")

# ==============================================================================
# 7. PRACTICAL RECOMMENDATIONS
# ==============================================================================
cat("7. Practical Recommendations\n")
cat("=============================\n")

cat("Based on this analysis:\n\n")

cat("Best distance method: ", best_result$distance_method, "\n")
cat("Best clustering method: ", best_result$clustering_method, "\n")
cat("Optimal number of clusters: ", best_result$k, "\n")

cat("\nMethod recommendations:\n")
if (best_result$distance_method == "euclidean") {
  cat("• Euclidean distance worked best - suggests numeric encoding captures important patterns\n")
} else if (best_result$distance_method == "hamming") {
  cat("• Hamming distance worked best - position-wise differences are important\n")
} else if (best_result$distance_method == "lcs") {
  cat("• LCS distance worked best - common subsequences are important\n")
}

if (best_result$clustering_method == "pam") {
  cat("• PAM clustering worked best - robust method, works well with any distance\n")
} else if (best_result$clustering_method == "ward.D2") {
  cat("• Ward clustering worked best - suggests compact, spherical clusters\n")
}

cat("\nCluster quality assessment:\n")
if (best_result$silhouette > 0.5) {
  cat("• Excellent cluster separation (silhouette > 0.5)\n")
} else if (best_result$silhouette > 0.3) {
  cat("• Good cluster separation (silhouette > 0.3)\n")
} else {
  cat("• Moderate cluster separation - consider different approaches\n")
}

if (best_result$cluster_balance > 0.3) {
  cat("• Well-balanced clusters\n")
} else {
  cat("• Unbalanced clusters - some clusters much larger than others\n")
}

cat("\nNext steps:\n")
cat("• Use the optimal parameters for your full dataset\n")
cat("• Consider MMM if temporal dynamics are important\n")
cat("• Validate clusters with domain knowledge\n")
cat("• Use cluster assignments for further analysis\n")

cat("\n=== Comprehensive clustering analysis completed! ===\n")
cat("Access detailed results with: comprehensive_results$all_results\n")
cat("Best clustering assignments: optimal_clustering$assignments\n") 
# ==============================================================================
# DEMO: Mixture Markov Model Clustering with tnaExtras
# ==============================================================================
#
# This demo shows how to use the advanced Mixture Markov Model (MMM) clustering
# in tnaExtras. MMM clustering uses EM algorithm to fit probabilistic models
# that capture both initial state probabilities and transition dynamics.
#
# ==============================================================================

library(tnaExtras)

cat("=== tnaExtras Clustering Demo: Mixture Markov Models ===\n\n")

# Load example data
data(regulation_grouped)
cat("Loaded regulation_grouped: ", nrow(regulation_grouped), " sequences\n")
cat("Time points: ", ncol(regulation_grouped) - 1, "\n")
cat("Groups: ", unique(regulation_grouped$Group), "\n")
cat("States: ", length(unique(unlist(regulation_grouped[, -ncol(regulation_grouped)]))), "\n\n")

# Extract a subset for demonstration (MMM can be computationally intensive)
demo_data <- regulation_grouped[1:50, 1:15]  # First 50 sequences, 15 time points
cat("Demo subset: ", nrow(demo_data), " sequences × ", ncol(demo_data), " time points\n")
cat("States in demo: ", paste(sort(unique(unlist(demo_data))), collapse = ", "), "\n\n")

# ==============================================================================
# 1. BASIC MMM CLUSTERING
# ==============================================================================
cat("1. Basic MMM Clustering\n")
cat("=======================\n")

# Fit MMM with 2 clusters
cat("Fitting MMM with k=2 clusters...\n")
mmm_2 <- cluster_mmm(demo_data, k = 2, n_starts = 5, max_iter = 50, 
                     verbose = FALSE, seed = 123)

cat("✓ MMM clustering completed\n")
cat("  Converged: ", mmm_2$converged, "\n")
cat("  Log-likelihood: ", round(mmm_2$log_likelihood, 2), "\n")
cat("  AIC: ", round(mmm_2$aic, 2), "\n")
cat("  BIC: ", round(mmm_2$bic, 2), "\n")
cat("  Cluster sizes: ", paste(mmm_2$cluster_sizes, collapse = ", "), "\n")
cat("  Mixture weights: ", paste(round(mmm_2$mixture_weights, 3), collapse = ", "), "\n")
cat("  Entropy: ", mmm_2$entropy, " (0=perfect separation, 1=uniform)\n")

cat("\n")

# ==============================================================================
# 2. EXAMINING MMM RESULTS
# ==============================================================================
cat("2. Examining MMM Results\n")
cat("========================\n")

# Look at cluster assignments
cat("Cluster assignments (first 10 sequences):\n")
print(mmm_2$assignments[1:10])

# Look at posterior probabilities
cat("\nPosterior probabilities (first 5 sequences):\n")
print(round(mmm_2$responsibilities[1:5, ], 3))

# Examine Markov models for each cluster
cat("\nMarkov models for each cluster:\n")
for (i in 1:mmm_2$k) {
  cat("--- Cluster ", i, " ---\n")
  cat("Initial state probabilities:\n")
  print(round(mmm_2$models[[i]]$initial, 3))
  cat("Transition matrix (top 3x3):\n")
  transition_matrix <- mmm_2$models[[i]]$transition
  print(round(transition_matrix[1:min(3, nrow(transition_matrix)), 
                                1:min(3, ncol(transition_matrix))], 3))
  cat("\n")
}

# ==============================================================================
# 3. MODEL SELECTION WITH MULTIPLE K VALUES
# ==============================================================================
cat("3. Model Selection with Multiple K Values\n")
cat("==========================================\n")

# Test different numbers of clusters using BIC
cat("Testing k = 2, 3, 4 with BIC criterion...\n")
mmm_selection <- find_optimal_mmm(demo_data, 
                                 k_range = 2:4, 
                                 criterion = "bic",
                                 n_starts = 3,  # Fewer starts for speed
                                 max_iter = 30,
                                 verbose = FALSE,
                                 seed = 456)

cat("✓ Model selection completed\n")
cat("Optimal k: ", mmm_selection$optimal_k, "\n")
print(mmm_selection$comparison[, c("k", "bic", "aic", "avg_posterior_overall", "converged")])

cat("\n")

# ==============================================================================
# 4. COMPARING MMM WITH STANDARD CLUSTERING
# ==============================================================================
cat("4. Comparing MMM with Standard Clustering\n")
cat("==========================================\n")

# Compare MMM with PAM clustering
cat("Comparing MMM vs PAM clustering (k=3)...\n")

# MMM clustering
mmm_3 <- cluster_mmm(demo_data, k = 3, n_starts = 3, max_iter = 30, 
                     verbose = FALSE, seed = 789)

# PAM clustering
pam_3 <- cluster_sequences(demo_data, k = 3, distance_method = "euclidean")

cat("✓ Comparison completed\n")
cat("MMM results:\n")
cat("  Log-likelihood: ", round(mmm_3$log_likelihood, 2), "\n")
cat("  BIC: ", round(mmm_3$bic, 2), "\n")
cat("  Cluster sizes: ", paste(mmm_3$cluster_sizes, collapse = ", "), "\n")
cat("  Entropy: ", mmm_3$entropy, "\n")

cat("PAM results:\n")
cat("  Silhouette: ", round(pam_3$silhouette, 3), "\n")
cat("  Cluster sizes: ", paste(pam_3$sizes, collapse = ", "), "\n")

# Compare cluster assignments
agreement <- sum(mmm_3$assignments == pam_3$assignments) / length(mmm_3$assignments)
cat("Cluster assignment agreement: ", round(agreement * 100, 1), "%\n")

cat("\n")

# ==============================================================================
# 5. ADVANCED MMM ANALYSIS
# ==============================================================================
cat("5. Advanced MMM Analysis\n")
cat("========================\n")

# Analyze cluster characteristics
optimal_mmm <- mmm_selection$best_result

cat("Analyzing optimal MMM model (k=", optimal_mmm$k, "):\n")

# Cluster quality metrics
cat("Quality metrics:\n")
cat("  Average posterior per cluster: ", 
    paste(round(optimal_mmm$avg_posterior_per_cluster, 3), collapse = ", "), "\n")
cat("  Overall average posterior: ", round(optimal_mmm$avg_posterior_overall, 3), "\n")
cat("  Normalized entropy: ", optimal_mmm$entropy, "\n")

# State preferences by cluster
cat("\nState usage patterns by cluster:\n")
for (i in 1:optimal_mmm$k) {
  cluster_seqs <- which(optimal_mmm$assignments == i)
  if (length(cluster_seqs) > 0) {
    cluster_data <- demo_data[cluster_seqs, ]
    state_counts <- table(unlist(cluster_data))
    top_states <- head(sort(state_counts, decreasing = TRUE), 3)
    cat("  Cluster ", i, " (n=", length(cluster_seqs), "): ", 
        paste(names(top_states), "(", top_states, ")", collapse = ", "), "\n")
  }
}

cat("\n")

# ==============================================================================
# 6. PRACTICAL RECOMMENDATIONS
# ==============================================================================
cat("6. Practical Recommendations\n")
cat("=============================\n")

cat("MMM Clustering best practices:\n")
cat("• Use multiple random starts (n_starts >= 10) for robust results\n")
cat("• Start with BIC for model selection, use AIC for comparison\n")
cat("• Check convergence status - rerun with more iterations if needed\n")
cat("• Lower entropy values indicate better cluster separation\n")
cat("• Higher average posterior probabilities indicate confident assignments\n")

cat("\nWhen to use MMM vs. standard clustering:\n")
cat("• MMM: When temporal dynamics matter, probabilistic assignments needed\n")
cat("• Standard: When simpler distance-based clustering is sufficient\n")
cat("• MMM: Better for sequential/temporal data with meaningful state transitions\n")
cat("• Standard: Faster, more straightforward interpretation\n")

cat("\nMMM parameters guidance:\n")
cat("• n_starts: 10-20 for production, 3-5 for exploration\n")
cat("• max_iter: 100-500 for convergence, 30-50 for quick tests\n")
cat("• k_range: Test 2 to sqrt(n_sequences), use BIC/AIC for selection\n")

cat("\n=== MMM Demo completed successfully! ===\n")
cat("Explore the fitted models with: mmm_result$models[[cluster_number]]\n")
cat("Use ?cluster_mmm and ?find_optimal_mmm for detailed documentation.\n") 
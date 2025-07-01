# ==============================================================================
# DEMO: Basic Sequence Clustering with tnaExtras
# ==============================================================================
#
# This demo shows how to use the basic clustering functions in tnaExtras
# for sequence analysis and clustering. We demonstrate multiple distance 
# methods and clustering algorithms.
#
# ==============================================================================

library(tnaExtras)

cat("=== tnaExtras Clustering Demo: Basic Functionality ===\n\n")

# Load example data
data(engagement_data)
cat("Loaded engagement_data: ", nrow(engagement_data), " sequences\n")
cat("Time points: ", ncol(engagement_data) - 1, "\n")
cat("Groups: ", unique(engagement_data$Group), "\n\n")

# Extract a subset for demonstration
demo_data <- engagement_data[1:20, 1:10]  # First 20 sequences, 10 time points
cat("Demo subset: ", nrow(demo_data), " sequences × ", ncol(demo_data), " time points\n\n")

# ==============================================================================
# 1. DISTANCE MATRIX COMPUTATION
# ==============================================================================
cat("1. Computing Distance Matrices\n")
cat("===============================\n")

# Euclidean distance (fastest, good for exploration)
cat("Computing Euclidean distance...\n")
dist_euclidean <- compute_distance_matrix(demo_data, "euclidean")
cat("✓ Euclidean distance: ", length(dist_euclidean), " pairwise distances\n")

# Hamming distance (position-wise differences)
cat("Computing Hamming distance...\n")
dist_hamming <- compute_distance_matrix(demo_data, "hamming")
cat("✓ Hamming distance: ", length(dist_hamming), " pairwise distances\n")

# LCS distance (longest common subsequence)
cat("Computing LCS distance...\n")
dist_lcs <- compute_distance_matrix(demo_data, "lcs")
cat("✓ LCS distance: ", length(dist_lcs), " pairwise distances\n")

# Advanced transition distance (Markov-based)
cat("Computing transition distance...\n")
dist_transition <- compute_distance_matrix(demo_data, "transition")
cat("✓ Transition distance: ", length(dist_transition), " pairwise distances\n")

cat("\n")

# ==============================================================================
# 2. CLUSTERING WITH DIFFERENT ALGORITHMS
# ==============================================================================
cat("2. Clustering with Different Algorithms\n")
cat("========================================\n")

# PAM clustering (Partitioning Around Medoids)
cat("PAM clustering with Euclidean distance...\n")
clust_pam <- cluster_sequences(demo_data, k = 3, 
                              distance_method = "euclidean", 
                              clustering_method = "pam")
cat("✓ PAM clustering completed\n")
cat("  Silhouette score: ", round(clust_pam$silhouette, 3), "\n")
cat("  Cluster sizes: ", paste(clust_pam$sizes, collapse = ", "), "\n")

# Hierarchical clustering (Ward's method)
cat("\nHierarchical clustering with Hamming distance...\n")
clust_ward <- cluster_sequences(demo_data, k = 3, 
                               distance_method = "hamming", 
                               clustering_method = "ward.D2")
cat("✓ Ward clustering completed\n")
cat("  Silhouette score: ", round(clust_ward$silhouette, 3), "\n")
cat("  Cluster sizes: ", paste(clust_ward$sizes, collapse = ", "), "\n")

# Complete linkage clustering
cat("\nComplete linkage clustering with LCS distance...\n")
clust_complete <- cluster_sequences(demo_data, k = 3, 
                                   distance_method = "lcs", 
                                   clustering_method = "complete")
cat("✓ Complete linkage clustering completed\n")
cat("  Silhouette score: ", round(clust_complete$silhouette, 3), "\n")
cat("  Cluster sizes: ", paste(clust_complete$sizes, collapse = ", "), "\n")

cat("\n")

# ==============================================================================
# 3. ANALYZING MULTIPLE DISTANCE METHODS
# ==============================================================================
cat("3. Analyzing Multiple Distance Methods\n")
cat("======================================\n")

# Compare different distance measures
cat("Analyzing distances with multiple methods...\n")
distance_analysis <- analyze_distances(demo_data, 
                                     methods = c("euclidean", "hamming", "lcs", "transition"))

cat("✓ Distance analysis completed\n")
for (method in names(distance_analysis)) {
  dist_matrix <- distance_analysis[[method]]
  cat("  ", method, ": ", length(dist_matrix), " distances, range: [", 
      round(min(dist_matrix), 2), ", ", round(max(dist_matrix), 2), "]\n")
}

cat("\n")

# ==============================================================================
# 4. COMPARING CLUSTERING METHODS
# ==============================================================================
cat("4. Comparing Clustering Methods\n")
cat("===============================\n")

# Compare different clustering algorithms
cat("Comparing PAM vs Ward clustering...\n")
method_comparison <- compare_clustering_methods(demo_data, 
                                              k = 3,
                                              distance_method = "euclidean",
                                              clustering_methods = c("pam", "ward.D2", "complete"))

cat("✓ Method comparison completed\n")
print(method_comparison[, c("method", "silhouette", "cluster_balance")])

cat("\n")

# ==============================================================================
# 5. TESTING MULTIPLE K VALUES
# ==============================================================================
cat("5. Testing Multiple K Values\n")
cat("=============================\n")

# Test different numbers of clusters
cat("Testing k = 2, 3, 4 with PAM clustering...\n")
k_analysis <- find_clusters_range(demo_data, 
                                 distance_method = "euclidean",
                                 k_range = 2:4,
                                 clustering_method = "pam")

cat("✓ K-range analysis completed\n")
for (k_name in names(k_analysis$all_results)) {
  result <- k_analysis$all_results[[k_name]]
  cat("  k=", result$k, ": silhouette=", round(result$silhouette, 3), 
      ", sizes=(", paste(result$sizes, collapse = ","), ")\n")
}

cat("\n")

# ==============================================================================
# 6. PRACTICAL RECOMMENDATIONS
# ==============================================================================
cat("6. Practical Recommendations\n")
cat("=============================\n")

# Find the best silhouette score
best_k <- NULL
best_silhouette <- -1
for (k_name in names(k_analysis$all_results)) {
  result <- k_analysis$all_results[[k_name]]
  if (result$silhouette > best_silhouette) {
    best_silhouette <- result$silhouette
    best_k <- result$k
  }
}

cat("Best clustering solution:\n")
cat("✓ Optimal k: ", best_k, "\n")
cat("✓ Silhouette score: ", round(best_silhouette, 3), "\n")

# Performance comparison
cat("\nDistance method characteristics:\n")
cat("• Euclidean: Fastest, good for exploration\n")
cat("• Hamming: Simple position-wise differences\n")
cat("• LCS: Considers common subsequences\n")
cat("• Transition: Advanced Markov-based analysis\n")

cat("\nClustering method characteristics:\n")
cat("• PAM: Robust, works well with any distance\n")
cat("• Ward: Good for compact, spherical clusters\n")
cat("• Complete: Conservative, tends to create tight clusters\n")

cat("\n=== Demo completed successfully! ===\n")
cat("Try modifying the parameters above to explore different options.\n")
cat("Use ?cluster_sequences for detailed documentation.\n") 
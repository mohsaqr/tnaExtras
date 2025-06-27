# ==============================================================================
# MAIN USER-FACING FUNCTIONS FOR TNA CLUSTERING
# ==============================================================================
# 
# This module contains the main functions for sequence analysis and clustering.
# All functions are fully documented and optimized for performance.
#
# ==============================================================================

#' Compute Distance Matrix for Sequence Data
#'
#' Computes dissimilarity matrix for categorical sequence data using various methods.
#' Includes traditional sequence analysis methods, advanced transition analysis, 
#' optimal matching with dynamic programming, and comprehensive stringdist algorithms.
#'
#' @param data A data frame or matrix where rows are sequences (cases/individuals) and 
#'   columns are time points. Each cell should contain a categorical state value.
#'   Missing values (NA) are handled appropriately by each method. Minimum 2 sequences required.
#' @param method Character string specifying the distance method. Options are:
#'   \itemize{
#'     \item \code{"euclidean"} - Fast numeric encoding (recommended for exploration)
#'     \item \code{"hamming"} - Position-wise differences
#'     \item \code{"lcs"} - Longest Common Subsequence similarity
#'     \item \code{"start_position"} - First occurrence positions
#'     \item \code{"transition"} - Advanced Markov transition analysis
#'     \item \code{"optimal_matching"} - Optimal matching with constant substitution costs
#'     \item \code{"osa"} - Optimal String Alignment (stringdist)
#'     \item \code{"lv"} - Levenshtein distance (stringdist)
#'     \item \code{"dl"} - Damerau-Levenshtein distance (stringdist)
#'     \item \code{"qgram"} - Q-gram distance (stringdist)
#'     \item \code{"jaro"} - Jaro distance (stringdist)
#'     \item \code{"jw"} - Jaro-Winkler distance (stringdist)
#'     \item \code{"cosine"} - Cosine distance on q-grams (stringdist)
#'     \item \code{"jaccard"} - Jaccard distance on q-grams (stringdist)
#'     \item \code{"soundex"} - Soundex phonetic distance (stringdist)
#'     \item \code{"lcs_stringdist"} - LCS via stringdist package
#'   }
#' @param q Integer. Q-gram size for qgram, cosine, and jaccard methods (default: 2)
#' @param p Numeric. Prefix scaling factor for Jaro-Winkler method (default: 0.1)
#' @param first_order_weight Numeric. Weight for basic state transitions (only for transition method, default: 0.3)
#' @param second_order_weight Numeric. Weight for higher-order patterns (only for transition method, default: 0.2)
#' @param time_weighted_weight Numeric. Weight for temporal emphasis (only for transition method, default: 0.2)
#' @param complexity_weight Numeric. Weight for pattern diversity (only for transition method, default: 0.15)
#' @param persistence_weight Numeric. Weight for state stability (only for transition method, default: 0.15)
#' @param substitution_cost Numeric. Cost for substituting one state for another (only for optimal_matching method, default: 1)
#' @param indel_cost Numeric. Cost for insertions and deletions (only for optimal_matching method, default: 1)
#'
#' @return A distance object of class 'dist' suitable for clustering functions.
#'
#' @examples
#' \dontrun{
#' # Example with categorical sequence data
#' data <- data.frame(
#'   T1 = c("A", "B", "A", "C"),
#'   T2 = c("B", "A", "B", "A"),
#'   T3 = c("C", "C", "A", "B")
#' )
#' 
#' # Compute different distance matrices
#' dist_euclidean <- compute_distance_matrix(data, "euclidean")
#' dist_levenshtein <- compute_distance_matrix(data, "lv")
#' dist_jaro_winkler <- compute_distance_matrix(data, "jw", p = 0.2)
#' dist_qgram <- compute_distance_matrix(data, "qgram", q = 3)
#' }
#'
#' @export
compute_distance_matrix <- function(data, 
                                   method = "euclidean",
                                   q = 2,
                                   p = 0.1,
                                   first_order_weight = 0.3,
                                   second_order_weight = 0.2,
                                   time_weighted_weight = 0.2,
                                   complexity_weight = 0.15,
                                   persistence_weight = 0.15,
                                   substitution_cost = 1,
                                   indel_cost = 1) {
  
  # Input validation
  if (!is.data.frame(data) && !is.matrix(data)) {
    stop("Error: 'data' must be a data frame or matrix. Received: ", class(data)[1])
  }
  
  if (nrow(data) < 2) {
    stop("Error: Data must contain at least 2 sequences for distance computation. Current sequences: ", nrow(data))
  }
  
  if (ncol(data) < 1) {
    stop("Error: Data must contain at least 1 time point. Current time points: ", ncol(data))
  }
  
  # Additional validation for parameter ranges
  if (!is.numeric(q) || q < 1) {
    stop("Error: 'q' parameter must be a positive integer. Received: ", q)
  }
  
  if (!is.numeric(p) || p < 0 || p > 1) {
    stop("Error: 'p' parameter must be between 0 and 1. Received: ", p)
  }
  
  if (!is.numeric(substitution_cost) || substitution_cost < 0) {
    stop("Error: 'substitution_cost' must be a non-negative number. Received: ", substitution_cost)
  }
  
  if (!is.numeric(indel_cost) || indel_cost < 0) {
    stop("Error: 'indel_cost' must be a non-negative number. Received: ", indel_cost)
  }
  
  # Expanded method validation
  valid_methods <- c("euclidean", "hamming", "lcs", "start_position", "transition",
                     "optimal_matching", "osa", "lv", "dl", "qgram", "jaro", "jw", "cosine", 
                     "jaccard", "soundex", "lcs_stringdist")
  
  method <- match.arg(method, valid_methods)
  
  # Compute distance based on method
  switch(method,
    "euclidean" = compute_euclidean_distance(data),
    "hamming" = compute_hamming_distance(data),
    "lcs" = compute_lcs_distance(data),
    "start_position" = compute_start_position_distance(data),
    "transition" = compute_transition_distance(data, 
                                              first_order_weight,
                                              second_order_weight,
                                              time_weighted_weight,
                                              complexity_weight,
                                              persistence_weight),
    "optimal_matching" = compute_optimal_matching_distance(data, substitution_cost, indel_cost),
    # Stringdist methods
    "osa" = compute_osa_distance(data),
    "lv" = compute_lv_distance(data),
    "dl" = compute_dl_distance(data),
    "qgram" = compute_qgram_distance(data, q),
    "jaro" = compute_jaro_distance(data),
    "jw" = compute_jw_distance(data, p),
    "cosine" = compute_cosine_distance(data, q),
    "jaccard" = compute_jaccard_distance(data, q),
    "soundex" = compute_soundex_distance(data),
    "lcs_stringdist" = compute_lcs_stringdist_distance(data)
  )
}

#' Cluster Sequences Using Distance Matrix
#'
#' Performs clustering on sequence data using specified distance and clustering methods.
#' Supports both hierarchical and partitioning clustering algorithms.
#'
#' @param data A data frame or matrix where rows are sequences and columns are time points
#' @param k Integer. Number of clusters to create
#' @param distance_method Character string specifying distance method (see \code{\link{compute_distance_matrix}})
#' @param clustering_method Character string specifying clustering method. Options are:
#'   \itemize{
#'     \item \code{"pam"} - Partitioning Around Medoids (default)
#'     \item \code{"ward.D2"} - Ward's method (recommended for hierarchical)
#'     \item \code{"ward.D"} - Original Ward's method
#'     \item \code{"complete"} - Complete linkage
#'     \item \code{"average"} - Average linkage (UPGMA)
#'     \item \code{"single"} - Single linkage
#'     \item \code{"mcquitty"}, \code{"median"}, \code{"centroid"} - Additional hierarchical methods
#'   }
#' @param ... Additional arguments passed to \code{\link{compute_distance_matrix}}
#'
#' @return A list containing:
#'   \item{assignments}{Integer vector of cluster assignments}
#'   \item{silhouette}{Silhouette score measuring clustering quality}
#'   \item{sizes}{Table of cluster sizes}
#'   \item{k}{Number of clusters}
#'   \item{method}{Clustering method used}
#'   \item{distance_matrix}{The computed distance matrix}
#'
#' @examples
#' \dontrun{
#' # Example clustering
#' data <- data.frame(
#'   T1 = c("A", "B", "A", "C", "A", "B"),
#'   T2 = c("B", "A", "B", "A", "C", "A"),
#'   T3 = c("C", "C", "A", "B", "B", "C")
#' )
#' 
#' # PAM clustering with euclidean distance (default)
#' result <- cluster_sequences(data, k = 2)
#' print(result$assignments)
#' print(result$silhouette)
#' }
#'
#' @export
cluster_sequences <- function(data, 
                             k, 
                             distance_method = "euclidean", 
                             clustering_method = "pam",
                             ...) {
  
  # Input validation
  if (!is.numeric(k) || length(k) != 1 || is.na(k) || k != round(k)) {
    stop("Error: 'k' must be a single integer value. Received: ", k)
  }
  
  if (k < 2 || k >= nrow(data)) {
    stop("Error: 'k' must be between 2 and ", nrow(data) - 1, ". Received: ", k, 
         " (for ", nrow(data), " sequences)")
  }
  
  clustering_method <- match.arg(clustering_method, 
                                c("pam", "ward.D2", "ward.D", "complete", "average", 
                                  "single", "mcquitty", "median", "centroid"))
  
  # Compute distance matrix
  dist_matrix <- compute_distance_matrix(data, distance_method, ...)
  
  # Perform clustering
  if (clustering_method == "pam") {
    clust_result <- cluster::pam(dist_matrix, k = k)
    assignments <- clust_result$clustering
    silhouette_score <- clust_result$silinfo$avg.width
  } else {
    hc <- stats::hclust(dist_matrix, method = clustering_method)
    assignments <- stats::cutree(hc, k = k)
    sil <- cluster::silhouette(assignments, dist_matrix)
    silhouette_score <- mean(sil[, 3])
  }
  
  cluster_sizes <- table(assignments)
  
  result <- list(
    assignments = assignments,
    silhouette = silhouette_score,
    sizes = cluster_sizes,
    k = k,
    method = clustering_method,
    distance_matrix = dist_matrix
  )
  
  # Add class for custom printing
  class(result) <- c("tna_cluster_result", "list")
  
  return(result)
}

#' Find Clusters for Specified Parameters
#'
#' Performs clustering for a specified range of k values without automatic optimization.
#' Returns results for all k values to allow manual inspection and selection.
#'
#' @param data A data frame or matrix where rows are sequences and columns are time points
#' @param distance_method Character string specifying distance method (see \code{\link{compute_distance_matrix}})
#' @param k_range Integer vector specifying the range of k values to test (default: 2:6)
#' @param clustering_method Character string specifying clustering method (see \code{\link{cluster_sequences}})
#' @param ... Additional arguments passed to \code{\link{compute_distance_matrix}}
#'
#' @return A list containing:
#'   \item{all_results}{List of clustering results for all tested k values}
#'   \item{distance_matrix}{The computed distance matrix}
#'   \item{clustering_method}{The clustering method used}
#'
#' @examples
#' \dontrun{
#' # Test multiple k values
#' data <- data.frame(
#'   T1 = c("A", "B", "A", "C", "A", "B", "C", "A"),
#'   T2 = c("B", "A", "B", "A", "C", "A", "A", "C"),
#'   T3 = c("C", "C", "A", "B", "B", "C", "B", "B")
#' )
#' 
#' results <- find_clusters_range(data, "euclidean", 2:4, "pam")
#' # Manually inspect results$all_results$k2, results$all_results$k3, etc.
#' }
#'
#' @export
find_clusters_range <- function(data, 
                                distance_method = "euclidean", 
                                k_range = 2:6, 
                                clustering_method = "pam",
                                ...) {
  
  # Input validation
  if (!is.numeric(k_range) || any(is.na(k_range)) || any(k_range != round(k_range))) {
    stop("Error: 'k_range' must contain only integer values. Received: ", paste(k_range, collapse = ", "))
  }
  
  k_range <- sort(unique(k_range))
  if (any(k_range < 2) || any(k_range >= nrow(data))) {
    invalid_k <- k_range[k_range < 2 | k_range >= nrow(data)]
    stop("Error: All values in k_range must be between 2 and ", nrow(data) - 1, 
         ". Invalid values: ", paste(invalid_k, collapse = ", "))
  }
  
  # Compute distance matrix once
  dist_matrix <- compute_distance_matrix(data, distance_method, ...)
  
  results <- list()
  
  for (k in k_range) {
    # Perform clustering
    if (clustering_method == "pam") {
      clust_result <- cluster::pam(dist_matrix, k = k)
      assignments <- clust_result$clustering
      silhouette_score <- clust_result$silinfo$avg.width
    } else {
      hc <- stats::hclust(dist_matrix, method = clustering_method)
      assignments <- stats::cutree(hc, k = k)
      sil <- cluster::silhouette(assignments, dist_matrix)
      silhouette_score <- mean(sil[, 3])
    }
    
    cluster_sizes <- table(assignments)
    
    result <- list(
      assignments = assignments,
      silhouette = silhouette_score,
      sizes = cluster_sizes,
      k = k,
      method = clustering_method
    )
    
    results[[paste0("k", k)]] <- result
  }
  
  list(
    all_results = results,
    distance_matrix = dist_matrix,
    clustering_method = clustering_method
  )
}

#' Compare Different Clustering Methods
#'
#' Compares multiple clustering algorithms using the same distance matrix and number of clusters.
#' Returns a comprehensive dataframe with clustering quality metrics for comparison.
#'
#' @param data A data frame or matrix where rows are sequences and columns are time points
#' @param distance_method Character string specifying distance method (see \code{\link{compute_distance_matrix}})
#' @param k Integer. Number of clusters to create
#' @param clustering_methods Character vector of clustering methods to compare (default: c("pam", "ward.D2", "complete", "average"))
#' @param ... Additional arguments passed to \code{\link{compute_distance_matrix}}
#'
#' @return A dataframe containing clustering quality metrics:
#'   \item{method}{Clustering method name}
#'   \item{silhouette}{Average silhouette score}
#'   \item{largest_cluster}{Size of largest cluster}
#'   \item{smallest_cluster}{Size of smallest cluster}
#'   \item{cluster_balance}{Ratio of smallest to largest cluster (balance measure)}
#'   \item{total_clusters}{Number of clusters}
#'   \item{cluster_variance}{Variance in cluster sizes}
#'   \item{assignments}{List column containing cluster assignments}
#'   \item{cluster_sizes}{List column containing cluster sizes as vectors}
#'   \item{cluster_proportions}{List column containing cluster proportions as vectors}
#'
#' @examples
#' \dontrun{
#' # Compare clustering methods
#' data <- data.frame(
#'   T1 = c("A", "B", "A", "C", "A", "B"),
#'   T2 = c("B", "A", "B", "A", "C", "A"),
#'   T3 = c("C", "C", "A", "B", "B", "C")
#' )
#' 
#' comparison <- compare_clustering_methods(data, "euclidean", 2)
#' print(comparison)
#' }
#'
#' @export
compare_clustering_methods <- function(data, 
                                      distance_method = "euclidean", 
                                      k = 3,
                                      clustering_methods = c("pam", "ward.D2"),
                                      ...) {
  
  # Input validation
  if (!is.numeric(k) || length(k) != 1 || is.na(k) || k != round(k)) {
    stop("Error: 'k' must be a single integer value. Received: ", k)
  }
  
  if (k < 2 || k >= nrow(data)) {
    stop("Error: 'k' must be between 2 and ", nrow(data) - 1, ". Received: ", k, 
         " (for ", nrow(data), " sequences)")
  }
  
  valid_methods <- c("pam", "ward.D2", "ward.D", "complete", "average", 
                     "single", "mcquitty", "median", "centroid")
  
  if (!all(clustering_methods %in% valid_methods)) {
    invalid_methods <- clustering_methods[!clustering_methods %in% valid_methods]
    stop("Error: Invalid clustering method(s): ", paste(invalid_methods, collapse = ", "), 
         ". Valid options are: ", paste(valid_methods, collapse = ", "))
  }
  
  # Compute distance matrix once
  dist_matrix <- compute_distance_matrix(data, distance_method, ...)
  
  # Initialize results dataframe
  results_df <- data.frame(
    method = character(),
    silhouette = numeric(),
    largest_cluster = integer(),
    smallest_cluster = integer(),
    cluster_balance = numeric(),
    total_clusters = integer(),
    cluster_variance = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Store assignments, sizes, and proportions separately
  assignments_list <- list()
  sizes_list <- list()
  proportions_list <- list()
  
  for (method in clustering_methods) {
    # Perform clustering
    if (method == "pam") {
      clust_result <- cluster::pam(dist_matrix, k = k)
      assignments <- clust_result$clustering
      silhouette_score <- clust_result$silinfo$avg.width
    } else {
      hc <- stats::hclust(dist_matrix, method = method)
      assignments <- stats::cutree(hc, k = k)
      sil <- cluster::silhouette(assignments, dist_matrix)
      silhouette_score <- mean(sil[, 3])
    }
    
    cluster_sizes <- table(assignments)
    largest_cluster <- max(cluster_sizes)
    smallest_cluster <- min(cluster_sizes)
    cluster_balance <- smallest_cluster / largest_cluster
    total_clusters <- length(cluster_sizes)
    cluster_variance <- var(as.numeric(cluster_sizes))
    
    # Calculate sizes and proportions as vectors
    sizes_vector <- as.numeric(cluster_sizes)
    proportions_vector <- round(sizes_vector / sum(sizes_vector), 2)
    
    # Add to results dataframe
    results_df <- rbind(results_df, data.frame(
      method = method,
      silhouette = silhouette_score,
      largest_cluster = largest_cluster,
      smallest_cluster = smallest_cluster,
      cluster_balance = cluster_balance,
      total_clusters = total_clusters,
      cluster_variance = cluster_variance,
      stringsAsFactors = FALSE
    ))
    
    # Store assignments, sizes, and proportions separately
    assignments_list[[method]] <- assignments
    sizes_list[[method]] <- sizes_vector
    proportions_list[[method]] <- proportions_vector
  }
  
  # Add assignments, sizes, and proportions as list columns
  results_df$assignments <- I(assignments_list[results_df$method])
  results_df$cluster_sizes <- I(sizes_list[results_df$method])
  results_df$cluster_proportions <- I(proportions_list[results_df$method])
  
  # Order by silhouette score (descending)
  results_df <- results_df[order(results_df$silhouette, decreasing = TRUE), ]
  rownames(results_df) <- NULL
  
  return(results_df)
}

#' Analyze Multiple Distance Methods for Sequences
#'
#' Computes distance matrices using multiple methods for comparison purposes.
#' This function is useful for exploring which distance measure works best for your data.
#'
#' @param data A data frame or matrix where rows are sequences and columns are time points
#' @param methods Character vector of distance methods to compute (default: c("euclidean", "lcs", "transition"))
#' @param ... Additional arguments passed to \code{\link{compute_distance_matrix}} (applies to all methods)
#'
#' @return A named list of distance matrices, one for each method
#'
#' @examples
#' \dontrun{
#' # Analyze multiple distance methods
#' data <- data.frame(
#'   T1 = c("A", "B", "A", "C"),
#'   T2 = c("B", "A", "B", "A"),
#'   T3 = c("C", "C", "A", "B")
#' )
#' 
#' distances <- analyze_distances(data, c("euclidean", "lcs", "hamming"))
#' # Use distances$euclidean, distances$lcs, etc.
#' }
#'
#' @export
analyze_distances <- function(data, 
                             methods = c("euclidean", "lcs", "transition"),
                             ...) {
  
  # Input validation
  valid_methods <- c("euclidean", "hamming", "lcs", "start_position", "transition",
                     "optimal_matching", "osa", "lv", "dl", "qgram", "jaro", "jw", "cosine", 
                     "jaccard", "soundex", "lcs_stringdist")
  
  if (!all(methods %in% valid_methods)) {
    invalid_methods <- methods[!methods %in% valid_methods]
    stop("Error: Invalid distance method(s): ", paste(invalid_methods, collapse = ", "), 
         ". Valid options are: ", paste(valid_methods, collapse = ", "))
  }
  
  results <- list()
  
  for (method in methods) {
    results[[method]] <- compute_distance_matrix(data, method, ...)
  }
  
  results
}

#' Complete Sequence Analysis and Clustering
#'
#' Performs comprehensive analysis by testing multiple distance and clustering methods.
#' Returns all results without automatic selection to allow manual inspection.
#' Supports parallel processing for faster computation.
#'
#' @param data A data frame or matrix where rows are sequences and columns are time points
#' @param k_range Integer vector specifying the range of k values to test (default: 2:5)
#' @param distance_methods Character vector of distance methods to test (default: c("euclidean", "lcs", "transition"))
#' @param clustering_methods Character vector of clustering methods to test (default: c("pam", "ward.D2"))
#' @param n_cores Integer. Number of cores for parallel processing. If NULL (default),
#'   uses sequential processing. If specified, automatically enables parallel processing.
#' @param min_cluster_size Integer. Minimum acceptable cluster size for filtering (default: 1)
#' @param balance_threshold Numeric. Minimum balance ratio for cluster quality assessment (default: 0.1)  
#' @param na_action Character. How to handle NAs: "ignore" (default), "remove", "error"
#' @param ... Additional arguments passed to \code{\link{compute_distance_matrix}}
#'
#' @return A list containing:
#'   \item{all_results}{Nested list of all results organized by distance method and clustering method}
#'   \item{summary}{Dataframe summarizing all combinations with key metrics}
#'
#' @examples
#' \dontrun{
#' # Complete analysis
#' data <- data.frame(
#'   T1 = c("A", "B", "A", "C", "A", "B", "C", "A"),
#'   T2 = c("B", "A", "B", "A", "C", "A", "A", "C"),
#'   T3 = c("C", "C", "A", "B", "B", "C", "B", "B")
#' )
#' 
#' # Sequential analysis
#' analysis <- cluster_complete_analysis(data)
#' 
#' # Parallel analysis (faster for multiple distance methods)
#' analysis_parallel <- cluster_complete_analysis(data, n_cores = 4)
#' print(analysis$summary)
#' }
#'
#' @export
cluster_complete_analysis <- function(data, 
                                      k_range = 2:5,
                                      distance_methods = c("euclidean", "lcs", "transition"),
                                      clustering_methods = c("pam", "ward.D2"),
                                      n_cores = NULL,
                                      min_cluster_size = 1,
                                      balance_threshold = 0.1,
                                      na_action = c("ignore", "remove", "error"),
                                      ...) {
  
  # Helper function: Setup parallel processing with robust error handling
  .setup_parallel_cluster <- function(n_cores = NULL, verbose = TRUE) {
    if (is.null(n_cores)) {
      if (verbose) cat("Using sequential processing\n")
      return(NULL)
    }
    
    if (!requireNamespace("parallel", quietly = TRUE)) {
      if (verbose) cat("Warning: parallel package not available, using sequential processing\n")
      return(NULL)
    }
    
    tryCatch({
      available_cores <- parallel::detectCores()
      if (is.na(available_cores) || available_cores < 2) {
        if (verbose) cat("Only 1 core detected, using sequential processing\n")
        return(NULL)
      }
      
      if (is.null(n_cores)) {
        n_cores <- max(1, available_cores - 1)
      } else {
        n_cores <- min(n_cores, available_cores - 1)
      }
      
      if (n_cores < 2) {
        if (verbose) cat("Less than 2 cores available for parallel processing, using sequential\n")
        return(NULL)
      }
      
      if (.Platform$OS.type == "windows") {
        cluster <- parallel::makeCluster(n_cores, type = "PSOCK")
      } else {
        cluster <- parallel::makeCluster(n_cores, type = "FORK")
      }
      
      test_result <- parallel::parLapply(cluster, 1:2, function(x) x^2)
      if (!identical(test_result, list(1, 4))) {
        parallel::stopCluster(cluster)
        if (verbose) cat("Parallel cluster test failed, using sequential processing\n")
        return(NULL)
      }
      
      if (verbose) {
        cat(sprintf("Parallel processing enabled: %d cores (%s)\n", 
                    n_cores, 
                    if(.Platform$OS.type == "windows") "PSOCK" else "FORK"))
      }
      
      return(list(
        cluster = cluster,
        n_cores = n_cores,
        type = if(.Platform$OS.type == "windows") "PSOCK" else "FORK"
      ))
      
    }, error = function(e) {
      if (verbose) {
        cat(sprintf("Parallel setup failed (%s), using sequential processing\n", 
                    e$message))
      }
      return(NULL)
    })
  }
  
  # Helper function: Clean up parallel cluster
  .cleanup_parallel_cluster <- function(cluster_info, verbose = TRUE) {
    if (is.null(cluster_info)) return()
    
    tryCatch({
      if (!is.null(cluster_info$cluster)) {
        parallel::stopCluster(cluster_info$cluster)
        if (verbose) cat("Parallel cluster stopped\n")
      }
    }, error = function(e) {
      if (verbose) cat("Warning: Error stopping parallel cluster:", e$message, "\n")
    })
  }
  
  # Helper function: Robust parallel lapply
  .robust_parallel_lapply <- function(X, FUN, cluster_info, ..., verbose = TRUE) {
    if (is.null(cluster_info)) {
      return(lapply(X, FUN, ...))
    }
    
    tryCatch({
      parallel::parLapply(cluster_info$cluster, X, FUN, ...)
    }, error = function(e) {
      if (verbose) {
        cat(sprintf("Parallel execution failed (%s), falling back to sequential\n", 
                    e$message))
      }
      lapply(X, FUN, ...)
    })
  }
  
  # Setup parallel processing
  cluster_info <- .setup_parallel_cluster(n_cores = n_cores, verbose = TRUE)
  on.exit(.cleanup_parallel_cluster(cluster_info, verbose = TRUE), add = TRUE)
  
  # Prepare worker function for processing each distance method
  process_distance_method <- function(distance_method) {
    method_results <- list()
    method_summary <- data.frame(
      distance_method = character(),
      clustering_method = character(),
      k = integer(),
      silhouette = numeric(),
      largest_cluster = integer(),
      smallest_cluster = integer(),
      cluster_balance = numeric(),
      stringsAsFactors = FALSE
    )
    method_sizes_list <- list()
    method_proportions_list <- list()
    
    for (clustering_method in clustering_methods) {
      result <- find_clusters_range(data, distance_method, k_range, clustering_method, ...)
      method_results[[clustering_method]] <- result
      
      # Extract summary information for each k
      for (k_name in names(result$all_results)) {
        k_result <- result$all_results[[k_name]]
        cluster_sizes <- k_result$sizes
        largest_cluster <- max(cluster_sizes)
        smallest_cluster <- min(cluster_sizes)
        cluster_balance <- smallest_cluster / largest_cluster
        
        # Calculate sizes and proportions as vectors
        sizes_vector <- as.numeric(cluster_sizes)
        proportions_vector <- round(sizes_vector / sum(sizes_vector), 2)
        
        # Create unique key for this combination
        combo_key <- paste(distance_method, clustering_method, k_result$k, sep = "_")
        method_sizes_list[[combo_key]] <- sizes_vector
        method_proportions_list[[combo_key]] <- proportions_vector
        
        method_summary <- rbind(method_summary, data.frame(
          distance_method = distance_method,
          clustering_method = clustering_method,
          k = k_result$k,
          silhouette = k_result$silhouette,
          largest_cluster = largest_cluster,
          smallest_cluster = smallest_cluster,
          cluster_balance = cluster_balance,
          stringsAsFactors = FALSE
        ))
      }
    }
    
    return(list(
      results = method_results,
      summary = method_summary,
      sizes_list = method_sizes_list,
      proportions_list = method_proportions_list
    ))
  }
  
  # Execute parallel or sequential processing
  if (!is.null(cluster_info)) {
    cat("Processing distance methods in parallel...\n")
    all_method_results <- .robust_parallel_lapply(distance_methods, process_distance_method, cluster_info, verbose = TRUE)
  } else {
    all_method_results <- lapply(distance_methods, process_distance_method)
  }
  
  # Combine results
  all_results <- list()
  summary_df <- data.frame(
    distance_method = character(),
    clustering_method = character(),
    k = integer(),
    silhouette = numeric(),
    largest_cluster = integer(),
    smallest_cluster = integer(),
    cluster_balance = numeric(),
    stringsAsFactors = FALSE
  )
  sizes_summary_list <- list()
  proportions_summary_list <- list()
  
  for (i in seq_along(distance_methods)) {
    distance_method <- distance_methods[i]
    method_data <- all_method_results[[i]]
    
    all_results[[distance_method]] <- method_data$results
    summary_df <- rbind(summary_df, method_data$summary)
    sizes_summary_list <- c(sizes_summary_list, method_data$sizes_list)
    proportions_summary_list <- c(proportions_summary_list, method_data$proportions_list)
  }
  
  # Order summary by silhouette score
  summary_df <- summary_df[order(summary_df$silhouette, decreasing = TRUE), ]
  rownames(summary_df) <- NULL
  
  # Add sizes and proportions as list columns
  combo_keys <- paste(summary_df$distance_method, summary_df$clustering_method, summary_df$k, sep = "_")
  summary_df$cluster_sizes <- I(sizes_summary_list[combo_keys])
  summary_df$cluster_proportions <- I(proportions_summary_list[combo_keys])
  
  list(
    all_results = all_results,
    summary = summary_df
  )
}

#' Print Method for TNA Cluster Results
#'
#' Custom print method that displays clustering results in a clean, readable format.
#'
#' @param x A tna_cluster_result object
#' @param ... Additional arguments (ignored)
#' @method print tna_cluster_result
#' @export
print.tna_cluster_result <- function(x, ...) {
  cat("=== TNA Clustering Result ===\n")
  cat("Method:", x$method, "\n")
  cat("Number of clusters (k):", x$k, "\n")
  cat("Silhouette score:", round(x$silhouette, 4), "\n")
  cat("Cluster sizes:", paste(as.numeric(x$sizes), collapse = ", "), "\n")
  
  # Calculate and display proportions
  proportions <- round(as.numeric(x$sizes) / sum(as.numeric(x$sizes)), 2)
  cat("Cluster proportions:", paste(proportions, collapse = ", "), "\n")
  
  cat("\nNote: Cluster assignments and distance matrix are stored in the object.\n")
  cat("Access them with result$assignments and result$distance_matrix\n")
  
  invisible(x)
}

#' Find Optimal Clusters with Advanced Suggestions
#'
#' Finds optimal number of clusters considering multiple criteria including
#' cluster balance, silhouette scores, and algorithm performance.
#'
#' @param data A data frame of sequences
#' @param distance_method Character. Distance method to use (default: "euclidean")
#' @param k_range Numeric vector. Range of k values to test (default: 2:6)
#' @param clustering_methods Character vector. Clustering methods to test
#' @param min_cluster_size Integer. Minimum acceptable cluster size (default: 1)
#' @param balance_threshold Numeric. Minimum balance ratio to consider acceptable (default: 0.1)
#' @param silhouette_weight Numeric. Weight for silhouette score in scoring (default: 0.4)
#' @param balance_weight Numeric. Weight for balance score in scoring (default: 0.3)
#' @param size_weight Numeric. Weight for size appropriateness in scoring (default: 0.3)
#' @param na_action Character. How to handle NAs: "ignore", "remove", "error" (default: "ignore")
#' @param ... Additional arguments passed to clustering functions
#'
#' @return A list containing optimal suggestions and detailed comparison
#' @export
find_optimal_clusters <- function(data,
                                 distance_method = "euclidean",
                                 k_range = 2:6,
                                 clustering_methods = c("pam", "ward.D2"),
                                 min_cluster_size = 1,
                                 balance_threshold = 0.1,
                                 silhouette_weight = 0.4,
                                 balance_weight = 0.3,
                                 size_weight = 0.3,
                                 na_action = c("ignore", "remove", "error"),
                                 ...) {
  
  na_action <- match.arg(na_action)
  
  # Handle NAs based on na_action
  if (na_action == "error") {
    if (any(is.na(data))) {
      stop("Data contains NAs and na_action='error' was specified")
    }
  } else if (na_action == "remove") {
    # Remove rows with any NAs
    complete_rows <- complete.cases(data)
    if (sum(complete_rows) < nrow(data)) {
      warning(sprintf("Removing %d rows with NAs", nrow(data) - sum(complete_rows)))
      data <- data[complete_rows, , drop = FALSE]
    }
  }
  # If na_action == "ignore", pass data as-is to underlying functions
  
  n_obs <- nrow(data)
  results_list <- list()
  comparison_df <- data.frame(
    k = integer(),
    clustering_method = character(),
    silhouette = numeric(),
    min_cluster_size = integer(),
    max_cluster_size = integer(),
    balance_ratio = numeric(),
    acceptable_balance = logical(),
    composite_score = numeric(),
    recommendation = character(),
    stringsAsFactors = FALSE
  )
  
  cat("=== Optimal Cluster Analysis ===\n")
  cat(sprintf("Distance method: %s\n", distance_method))
  cat(sprintf("Data: %d observations\n", n_obs))
  cat(sprintf("Testing k in {%s}\n", paste(range(k_range), collapse = "-")))
  cat(sprintf("Balance threshold: %.2f\n", balance_threshold))
  cat("\n")
  
  for (clustering_method in clustering_methods) {
    cat(sprintf("Testing %s clustering...\n", clustering_method))
    
    method_results <- find_clusters_range(data, distance_method, k_range, clustering_method, ...)
    results_list[[clustering_method]] <- method_results
    
    for (k_name in names(method_results$all_results)) {
      result <- method_results$all_results[[k_name]]
      k_val <- result$k
      
      # Calculate metrics
      cluster_sizes <- as.numeric(result$sizes)
      min_size <- min(cluster_sizes)
      max_size <- max(cluster_sizes)
      balance_ratio <- min_size / max_size
      silhouette_score <- result$silhouette
      
      # Check if clusters meet minimum size requirement
      acceptable_sizes <- all(cluster_sizes >= min_cluster_size)
      acceptable_balance <- balance_ratio >= balance_threshold
      
      # Calculate composite score
      # Normalize silhouette score (can be negative, so shift to 0-1 range)
      sil_normalized <- (silhouette_score + 1) / 2  # Shift from [-1,1] to [0,1]
      
      # Balance score (already 0-1)
      balance_score <- balance_ratio
      
      # Size appropriateness (penalize too many or too few clusters)
      ideal_k <- round(sqrt(n_obs / 2))  # Rough heuristic for ideal k
      size_score <- 1 - abs(k_val - ideal_k) / max(k_range)
      size_score <- max(0, size_score)  # Ensure non-negative
      
      composite_score <- (silhouette_weight * sil_normalized + 
                         balance_weight * balance_score + 
                         size_weight * size_score)
      
      # Generate recommendation
      recommendation <- "Poor"
      if (acceptable_sizes && acceptable_balance && silhouette_score > 0.5) {
        recommendation <- "Excellent"
      } else if (acceptable_sizes && acceptable_balance && silhouette_score > 0.3) {
        recommendation <- "Good"
      } else if (acceptable_sizes && silhouette_score > 0.3) {
        recommendation <- "Acceptable"
      } else if (acceptable_sizes) {
        recommendation <- "Marginal"
      }
      
      comparison_df <- rbind(comparison_df, data.frame(
        k = k_val,
        clustering_method = clustering_method,
        silhouette = round(silhouette_score, 3),
        min_cluster_size = min_size,
        max_cluster_size = max_size,
        balance_ratio = round(balance_ratio, 3),
        acceptable_balance = acceptable_balance,
        composite_score = round(composite_score, 3),
        recommendation = recommendation,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  # Sort by composite score
  comparison_df <- comparison_df[order(comparison_df$composite_score, decreasing = TRUE), ]
  
  # Find best overall
  best_overall <- comparison_df[1, ]
  
  # Find best for each criterion
  best_silhouette <- comparison_df[which.max(comparison_df$silhouette), ]
  best_balance <- comparison_df[which.max(comparison_df$balance_ratio), ]
  acceptable_options <- comparison_df[comparison_df$acceptable_balance, ]
  
  cat("=== Results Summary ===\n")
  print(comparison_df, row.names = FALSE)
  
  cat("\n=== Recommendations ===\n")
  cat(sprintf("[Best Overall] k=%d with %s (score=%.3f, sil=%.3f, bal=%.3f)\n", 
              best_overall$k, best_overall$clustering_method, 
              best_overall$composite_score, best_overall$silhouette, best_overall$balance_ratio))
  
  cat(sprintf("[Best Silhouette] k=%d with %s (silhouette=%.3f)\n", 
              best_silhouette$k, best_silhouette$clustering_method, best_silhouette$silhouette))
  
  cat(sprintf("[Best Balance] k=%d with %s (balance=%.3f)\n", 
              best_balance$k, best_balance$clustering_method, best_balance$balance_ratio))
  
  if (nrow(acceptable_options) > 0) {
    cat(sprintf("[OK] %d options meet balance threshold (>= %.2f)\n", 
                nrow(acceptable_options), balance_threshold))
  } else {
    cat(sprintf("[WARN] No options meet balance threshold (>= %.2f)\n", balance_threshold))
    cat("Consider: lowering balance_threshold, using fewer clusters, or different distance method\n")
  }
  
  # Algorithm suggestions
  cat("\n=== Algorithm Performance ===\n")
  algo_summary <- aggregate(composite_score ~ clustering_method, comparison_df, mean)
  algo_summary <- algo_summary[order(algo_summary$composite_score, decreasing = TRUE), ]
  for (i in 1:nrow(algo_summary)) {
    algo <- algo_summary$clustering_method[i]
    score <- algo_summary$composite_score[i]
    cat(sprintf("- %s: average score %.3f\n", algo, score))
  }
  
  return(list(
    best_overall = best_overall,
    best_silhouette = best_silhouette,
    best_balance = best_balance,
    comparison = comparison_df,
    acceptable_options = acceptable_options,
    algorithm_ranking = algo_summary,
    all_results = results_list,
    parameters = list(
      distance_method = distance_method,
      balance_threshold = balance_threshold,
      weights = c(silhouette = silhouette_weight, 
                 balance = balance_weight, 
                 size = size_weight)
    )
  ))
} 
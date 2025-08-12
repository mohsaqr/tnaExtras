# ==============================================================================
# SPECIALIZED SCRIPT FOR WEIGHTED SEQUENCE CLUSTERING
# ==============================================================================
# This script provides a focused, self-contained workflow for performing
# advanced sequence analysis using mandatory weights. It is not a general-purpose
# tool; it is designed specifically for calculating distances with either
# static weights or a dynamic "attention decay" model, and then using those
# distances for clustering.
#
# Required packages: cluster
# You can install it with: install.packages("cluster")
# ==============================================================================


# ==============================================================================
# PART 1: DEDICATED WEIGHTED DISTANCE FUNCTION
# ==============================================================================

#' Compute Weighted Sequence Distances
#'
#' This function calculates sequence distances using mandatory weights.
#' It can operate in two weighted modes for "euclidean" and "hamming" methods:
#' 1. Statically Weighted: `time_weights` and `state_weights` are used directly.
#' 2. Dynamically Weighted: A `decay_rate` is provided, which uses the static
#'    weights as a base to generate dynamic weights via an attention decay model.
#'
#' @param data The sequence data as a data frame or matrix.
#' @param method The distance metric to use: "euclidean" or "hamming".
#' @param time_weights A numeric vector of static weights for each time point.
#' @param state_weights A matrix of static weights for each state, per sequence.
#' @param decay_rate (Optional) If provided, this activates the dynamic attention
#'                   decay model. It's a positive number controlling how quickly
#'                   influence fades.
#' @return A distance object of class 'dist'.
compute_weighted_distance <- function(data,
                                      method,
                                      time_weights,
                                      state_weights,
                                      decay_rate = NULL) {

  # --- 1. Input Validation and Setup ---
  char_matrix <- as.matrix(data)
  n_seq <- nrow(char_matrix)
  n_time <- ncol(char_matrix)
  all_states <- sort(unique(as.vector(char_matrix[!is.na(char_matrix)])))
  n_states <- length(all_states)
  
  if (!method %in% c("euclidean", "hamming")) {
    stop("This function is designed for weighted 'euclidean' or 'hamming' distances.")
  }
  if (missing(time_weights) || missing(state_weights)) {
    stop("'time_weights' and 'state_weights' are required arguments.")
  }

  # --- 2. Generate the Final Weight Matrix ---
  state_weights <- state_weights[, all_states, drop = FALSE]
  final_weights <- NULL

  if (!is.null(decay_rate)) {
    # --- DYNAMIC MODEL: Use the attention decay model ---
    if (decay_rate < 0) stop("decay_rate must be non-negative.")
    
    state_indices <- matrix(match(char_matrix, all_states), nrow = n_seq)
    initial_influence <- state_weights[cbind(rep(1:n_seq, times = n_time), as.vector(state_indices))]
    initial_influence <- matrix(initial_influence, nrow = n_seq, ncol = n_time) * rep(time_weights, each = n_seq)
    initial_influence[is.na(char_matrix)] <- 0
    
    attention_landscape <- array(0, dim = c(n_seq, n_time, n_states), dimnames = list(NULL, NULL, all_states))
    for (i in 1:n_seq) {
      for (t in 1:n_time) {
        if (initial_influence[i, t] == 0) next
        current_state <- char_matrix[i, t]
        for (k in t:n_time) {
          decayed_influence <- initial_influence[i, t] * exp(-decay_rate * (k - t))
          attention_landscape[i, k, current_state] <- attention_landscape[i, k, current_state] + decayed_influence
        }
      }
    }
    
    final_weights <- matrix(0, nrow = n_seq, ncol = n_time)
    for (i in 1:n_seq) {
      for (t in 1:n_time) {
        current_state <- char_matrix[i, t]
        if (!is.na(current_state)) final_weights[i, t] <- attention_landscape[i, t, current_state]
      }
    }
    
  } else {
    # --- STATIC MODEL: Just combine the provided static weights ---
    state_indices <- matrix(match(char_matrix, all_states), nrow = n_seq)
    final_weights <- state_weights[cbind(rep(1:n_seq, times = n_time), as.vector(state_indices))]
    final_weights <- matrix(final_weights, nrow = n_seq, ncol = n_time) * rep(time_weights, each = n_seq)
    final_weights[is.na(char_matrix)] <- 0
  }

  # --- 3. Compute Final Distance ---
  if (method == "euclidean") {
    state_mapping <- seq_along(all_states); names(state_mapping) <- all_states
    numeric_matrix <- matrix(state_mapping[char_matrix], nrow = n_seq, ncol = n_time)
    numeric_matrix[is.na(numeric_matrix)] <- 0
    
    weighted_matrix <- numeric_matrix * final_weights
    return(stats::dist(weighted_matrix, method = "euclidean"))
    
  } else { # Hamming
    dist_matrix <- matrix(0, nrow = n_seq, ncol = n_seq)
    for (i in 1:(n_seq - 1)) {
      for (j in (i + 1):n_seq) {
        is_different <- char_matrix[i, ] != char_matrix[j, ]
        is_different[is.na(char_matrix[i, ]) & is.na(char_matrix[j, ])] <- FALSE
        
        avg_weights <- (final_weights[i, ] + final_weights[j, ]) / 2
        total_weighted_diff <- sum(avg_weights[is_different], na.rm = TRUE)
        dist_matrix[i, j] <- dist_matrix[j, i] <- total_weighted_diff
      }
    }
    return(as.dist(dist_matrix))
  }
}


# ==============================================================================
# PART 2: WEIGHTED CLUSTERING FUNCTION
# ==============================================================================

#' Cluster Sequences Using a Weighted Distance Matrix
#'
#' Performs clustering on sequence data using the dedicated weighted distance
#' function and a specified clustering algorithm.
#'
#' @param data A data frame or matrix where rows are sequences.
#' @param k Integer. The number of clusters to create.
#' @param distance_method Character string: "euclidean" or "hamming".
#' @param clustering_method Character string specifying clustering method (e.g., "pam", "ward.D2").
#' @param ... Required and optional arguments passed directly to `compute_weighted_distance`,
#'            such as `time_weights`, `state_weights`, and `decay_rate`.
#' @return A list containing clustering results.
cluster_weighted_sequences <- function(data, 
                                       k, 
                                       distance_method, 
                                       clustering_method = "pam",
                                       ...) {
  
  # --- 1. Input Validation ---
  if (missing(k) || missing(distance_method)) stop("'k' and 'distance_method' are required.")
  if (!is.numeric(k) || k < 2 || k >= nrow(data)) {
    stop("Error: 'k' must be an integer between 2 and ", nrow(data) - 1)
  }
  
  # --- 2. Compute Weighted Distance Matrix ---
  cat("Computing weighted distance matrix with method:", distance_method, "...\n")
  dist_matrix <- compute_weighted_distance(data, method = distance_method, ...)
  
  # --- 3. Perform Clustering ---
  cat("Performing clustering with method:", clustering_method, "...\n")
  
  if (clustering_method == "pam") {
    pam_result <- cluster::pam(dist_matrix, k = k)
    assignments <- pam_result$clustering
    silhouette_score <- pam_result$silinfo$avg.width
  } else {
    hc <- stats::hclust(dist_matrix, method = clustering_method)
    assignments <- stats::cutree(hc, k = k)
    sil <- cluster::silhouette(assignments, dist_matrix)
    silhouette_score <- mean(sil[, 3])
  }
  
  # --- 4. Compile and Return Results ---
  result <- list(
    assignments = assignments,
    silhouette = silhouette_score,
    sizes = table(assignments),
    k = k,
    distance_method = distance_method,
    clustering_method = clustering_method,
    distance_matrix = dist_matrix
  )
  
  cat("Clustering complete.\n")
  return(result)
}


# ==============================================================================
# PART 3: DEMONSTRATION
# ==============================================================================
cat("\n\n=== Specialized Weighted Clustering Demonstration ===\n\n")

# --- Example Data ---
tna_data_example <- data.frame(
  Time1 = c("A", "B", "C", "A", "A", "B"),
  Time2 = c("B", "C", NA, "B", "B", "A"),
  Time3 = c("C", "A", "B", "C", "C", "C"),
  Time4 = c("A", "B", "C", "A", "A", "B")
)

# --- Example Weights (Required) ---
state_weight_matrix <- matrix(c(
  2.0, 1.0, 1.5,  # Seq 1
  1.0, 2.5, 1.0,  # Seq 2
  1.5, 1.0, 2.0,  # Seq 3
  1.0, 1.0, 1.0,  # Seq 4
  2.5, 1.0, 1.0,  # Seq 5
  1.0, 2.0, 1.5   # Seq 6
), nrow = 6, byrow = TRUE)
colnames(state_weight_matrix) <- c("A", "B", "C")

time_weight_vector <- c(0.5, 1.0, 1.5, 2.0) # Emphasize later time points

# --- Perform Clustering with Dynamic Decay Weights ---
# We will find 3 clusters using a hamming distance that is dynamically weighted
# by our attention decay model.
clustering_result <- cluster_weighted_sequences(
  data = tna_data_example,
  k = 3,
  distance_method = "hamming",
  clustering_method = "pam",
  # These arguments are passed directly to compute_weighted_distance:
  time_weights = time_weight_vector,
  state_weights = state_weight_matrix,
  decay_rate = 0.75 # A fairly rapid decay
)

# --- Display Results ---
cat("\n--- Clustering Results ---\n")
cat("Number of clusters (k):", clustering_result$k, "\n")
cat("Clustering Algorithm:", clustering_result$clustering_method, "\n")
cat("Distance Metric:", clustering_result$distance_method, "(with dynamic decay weights)\n")
cat("Silhouette Score:", round(clustering_result$silhouette, 4), "\n\n")
cat("Cluster Assignments:\n")
print(clustering_result$assignments)
cat("\nCluster Sizes:\n")
print(clustering_result$sizes)

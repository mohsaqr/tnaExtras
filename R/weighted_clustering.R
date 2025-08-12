# ==============================================================================
# WEIGHTED SEQUENCE DISTANCES AND CLUSTERING
# ==============================================================================

#' Compute Weighted Sequence Distances
#'
#' Calculates sequence distances using mandatory weights. Supports two modes
#' for \code{"euclidean"} and \code{"hamming"} distance methods:
#' \itemize{
#'   \item Statically weighted: uses provided \code{time_weights} and \code{state_weights} directly
#'   \item Dynamically weighted: when \code{decay_rate} is provided, applies an
#'         attention decay model over time using the static weights as base
#' }
#'
#' @param data A data frame or matrix of sequences (rows = sequences, columns = time points)
#' @param method The distance metric to use: one of \code{"euclidean"} or \code{"hamming"}
#' @param time_weights Numeric vector of static weights for each time point (length = ncol(data))
#' @param state_weights Numeric matrix of static weights per state and sequence
#'   (nrow = nrow(data), ncol = number of unique states in \code{data}); columns must be named
#'   with the state labels present in \code{data}
#' @param decay_rate Optional non-negative numeric. If provided, enables the dynamic
#'   attention decay model controlling how quickly influence fades forward in time
#'
#' @return A distance object of class \code{dist}
#'
#' @examples
#' \dontrun{
#' data <- data.frame(
#'   T1 = c("A", "B", "C"),
#'   T2 = c("B", "C", "A"),
#'   T3 = c("C", "A", "B")
#' )
#' state_w <- matrix(c(1,2,1, 2,1,2, 1,1,2), nrow = 3, byrow = TRUE)
#' colnames(state_w) <- c("A","B","C")
#' time_w <- c(0.5, 1, 1.5)
#'
#' # Static weights
#' d1 <- compute_weighted_distance(data, method = "hamming",
#'                                 time_weights = time_w, state_weights = state_w)
#'
#' # Dynamic weights with decay
#' d2 <- compute_weighted_distance(data, method = "euclidean",
#'                                 time_weights = time_w, state_weights = state_w,
#'                                 decay_rate = 0.5)
#' }
#'
#' @export
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
  if (length(time_weights) != n_time) {
    stop("'time_weights' length (", length(time_weights), ") must equal number of time points (", n_time, ")")
  }
  if (!is.matrix(state_weights) || nrow(state_weights) != n_seq) {
    stop("'state_weights' must be a matrix with nrow equal to number of sequences (", n_seq, ")")
  }
  if (is.null(colnames(state_weights))) {
    stop("'state_weights' must have column names matching state labels in data")
  }

  # --- 2. Generate the Final Weight Matrix ---
  # Align state_weights columns to the states present in data
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
    return(stats::as.dist(dist_matrix))
  }
}


#' Cluster Sequences Using a Weighted Distance Matrix
#'
#' Performs clustering on sequence data using \code{compute_weighted_distance}
#' and a specified clustering algorithm.
#'
#' @param data A data frame or matrix where rows are sequences
#' @param k Integer. The number of clusters to create (2 .. nrow(data)-1)
#' @param distance_method Character string: one of \code{"euclidean"}, \code{"hamming"}
#' @param clustering_method Character string specifying clustering method (e.g., \code{"pam"}, \code{"ward.D2"})
#' @param ... Required and optional arguments passed directly to \code{compute_weighted_distance},
#'   such as \code{time_weights}, \code{state_weights}, and \code{decay_rate}
#'
#' @return A list containing clustering results: \code{assignments}, \code{silhouette}, \code{sizes},
#'   \code{k}, \code{distance_method}, \code{clustering_method}, and \code{distance_matrix}
#'
#' @examples
#' \dontrun{
#' data <- data.frame(
#'   T1 = c("A","B","C","A"),
#'   T2 = c("B","C", NA, "B"),
#'   T3 = c("C","A","B","C")
#' )
#' state_w <- matrix(c(2,1,1.5, 1,2.5,1, 1.5,1,2, 1,1,1), nrow = 4, byrow = TRUE)
#' colnames(state_w) <- c("A","B","C")
#' time_w <- c(0.5, 1, 1.5)
#' res <- cluster_weighted_sequences(data, k = 2, distance_method = "hamming",
#'                                   clustering_method = "pam",
#'                                   time_weights = time_w,
#'                                   state_weights = state_w,
#'                                   decay_rate = 0.75)
#' }
#'
#' @export
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
  dist_matrix <- compute_weighted_distance(data, method = distance_method, ...)

  # --- 3. Perform Clustering ---
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

  class(result) <- c("tna_cluster_result", "list")
  return(result)
}


#' Weighted Clustering (Alias)
#'
#' Alias for \code{cluster_weighted_sequences} matching user-facing naming.
#'
#' @inheritParams cluster_weighted_sequences
#' @inherit cluster_weighted_sequences return
#' @export
weighted_clustering <- function(data,
                                k,
                                distance_method,
                                clustering_method = "pam",
                                ...) {
  cluster_weighted_sequences(
    data = data,
    k = k,
    distance_method = distance_method,
    clustering_method = clustering_method,
    ...
  )
}


# ==============================================================================
# LINK PREDICTION FOR TNA NETWORKS
# ==============================================================================
#
# Functions for predicting missing or potential transitions in TNA networks
# using established link prediction algorithms.
#
# Functions:
# - predict_links(): Main function for link prediction
# - print.link_prediction(): Print method for results
# - summary.link_prediction(): Summary method for results
# - plot.link_prediction(): Visualization of predicted links
#
# ==============================================================================

#' Predict Links in a TNA Network
#'
#' Given a TNA model (transition network), predict which transitions are likely
#' to exist but haven't been observed, or score the strength of potential
#' connections between states.
#'
#' @param model A TNA model object (from tna package) or a transition matrix.
#'   Can also be a list containing a 'weights' or 'transition_matrix' element.
#' @param method Character string specifying the link prediction algorithm to use.
#'   Options are:
#'   \itemize{
#'     \item \code{"common_neighbors"} - Count shared successors/predecessors (fast baseline)
#'     \item \code{"jaccard"} - Normalized neighbor overlap
#'     \item \code{"adamic_adar"} - Weighted by inverse log degree
#'     \item \code{"preferential_attachment"} - Product of out-degree x in-degree
#'     \item \code{"katz"} - Weighted sum of all paths (global structure)
#'     \item \code{"hitting_time"} - Expected random walk steps
#'   }
#' @param top_n Integer. Return top N predicted links. Default is 10.
#' @param threshold Numeric. Minimum score threshold (alternative to top_n).
#'   If provided, overrides top_n.
#' @param exclude_existing Logical. If TRUE (default), only predict missing links.
#' @param damping Numeric. Damping factor for Katz/random walk methods. Default is 0.85.
#' @param max_iter Integer. Maximum iterations for iterative methods. Default is 100.
#' @param verbose Logical. Print progress messages. Default is TRUE.
#'
#' @return An object of class \code{link_prediction} containing:
#'   \item{predictions}{Data frame with columns: from, to, score, rank, exists}
#'   \item{method}{The method used}
#'   \item{network_stats}{List with n_states, n_edges, density}
#'   \item{model}{The original model}
#'
#' @details
#' Link prediction methods work differently:
#' \itemize{
#'   \item \strong{common_neighbors}: Simple and fast. Counts shared neighbors
#'     (successors and predecessors) between node pairs.
#'   \item \strong{jaccard}: Normalizes common neighbors by total neighbors.
#'     Better for comparing nodes with different degrees.
#'   \item \strong{adamic_adar}: Weights shared neighbors by the inverse log of
#'     their degree. Rare shared neighbors are more informative.
#'   \item \strong{preferential_attachment}: Product of out-degree and in-degree.
#'     Based on the idea that high-degree nodes attract more connections.
#'   \item \strong{katz}: Considers all paths between nodes, weighted by path length.
#'     Captures global network structure.
#'   \item \strong{hitting_time}: Based on random walk hitting times.
#'     Lower hitting time suggests stronger connection.
#' }
#'
#' @examples
#' \dontrun{
#' # Create a simple transition matrix
#' P <- matrix(c(0.7, 0.2, 0.1,
#'               0.3, 0.5, 0.2,
#'               0.1, 0.3, 0.6), nrow = 3, byrow = TRUE)
#' rownames(P) <- colnames(P) <- c("A", "B", "C")
#'
#' # Predict missing links using Katz index
#' pred <- predict_links(P, method = "katz", top_n = 5)
#' print(pred)
#'
#' # Compare methods
#' pred_aa <- predict_links(P, method = "adamic_adar")
#' pred_cn <- predict_links(P, method = "common_neighbors")
#' }
#'
#' @export
predict_links <- function(model,
                          method = c("katz", "adamic_adar", "jaccard",
                                    "common_neighbors", "hitting_time",
                                    "preferential_attachment"),
                          top_n = 10,
                          threshold = NULL,
                          exclude_existing = TRUE,
                          damping = 0.85,
                          max_iter = 100,
                          verbose = TRUE) {

  method <- match.arg(method)

  # Extract transition matrix from model (using shared function from markov_analysis.R)
  P <- extract_transition_matrix(model)

  if (is.null(P)) {
    stop("Could not extract transition matrix from model. ",
         "Please provide a matrix, a list with 'weights' or 'transition_matrix', ",
         "or a tna model object.")
  }

  # Validate matrix
  if (!is.matrix(P) || nrow(P) != ncol(P)) {
    stop("Transition matrix must be a square matrix.")
  }

  if (nrow(P) < 2) {
    stop("Transition matrix must have at least 2 states.")
  }

  # Ensure state names exist
  states <- rownames(P)
  if (is.null(states)) {
    states <- paste0("S", seq_len(nrow(P)))
    rownames(P) <- states
    colnames(P) <- states
  }

  # Create binary adjacency matrix
  A <- (P > 0) * 1
  n_states <- nrow(A)
  n_edges <- sum(A)
  max_edges <- n_states * n_states  # Including self-loops for transitions
  density <- n_edges / max_edges

  if (verbose) {
    cat("Link Prediction for TNA Network\n")
    cat("================================\n")
    cat(sprintf("States: %d\n", n_states))
    cat(sprintf("Existing edges: %d\n", n_edges))
    cat(sprintf("Network density: %.3f\n", density))
    cat(sprintf("Method: %s\n", method))
    cat("\n")
  }

  # Compute link prediction scores based on method
  scores <- switch(method,
    "common_neighbors" = .compute_common_neighbors(A),
    "jaccard" = .compute_jaccard(A),
    "adamic_adar" = .compute_adamic_adar(A),
    "preferential_attachment" = .compute_preferential_attachment(A),
    "katz" = .compute_katz(A, damping, max_iter, verbose),
    "hitting_time" = .compute_hitting_time(P, max_iter, verbose)
  )

  # Create predictions data frame
  predictions <- .create_predictions_df(scores, A, states, exclude_existing,
                                         top_n, threshold, verbose)

  # Build result object
  result <- list(
    predictions = predictions,
    method = method,
    network_stats = list(
      n_states = n_states,
      n_edges = n_edges,
      density = density
    ),
    model = model,
    scores_matrix = scores,
    adjacency = A,
    transition_matrix = P
  )

  class(result) <- c("link_prediction", "list")
  return(result)
}

# ==============================================================================
# INTERNAL HELPER FUNCTIONS
# ==============================================================================
# Note: extract_transition_matrix() is now a shared function in markov_analysis.R

#' Compute Common Neighbors score
#' @param A Binary adjacency matrix
#' @return Score matrix
#' @keywords internal
.compute_common_neighbors <- function(A) {
  n <- nrow(A)
  scores <- matrix(0, n, n)

  # For directed networks: shared successors + shared predecessors
  # Shared successors: A[i,] & A[j,]
  # Shared predecessors: A[,i] & A[,j]

  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      # Shared successors (both i and j connect to same nodes)
      shared_succ <- sum(A[i, ] & A[j, ])
      # Shared predecessors (same nodes connect to both i and j)
      shared_pred <- sum(A[, i] & A[, j])
      scores[i, j] <- shared_succ + shared_pred
    }
  }

  rownames(scores) <- rownames(A)
  colnames(scores) <- colnames(A)
  return(scores)
}

#' Compute Jaccard similarity score
#' @param A Binary adjacency matrix
#' @return Score matrix
#' @keywords internal
.compute_jaccard <- function(A) {
  n <- nrow(A)
  scores <- matrix(0, n, n)

  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      # Union and intersection for successors
      union_succ <- sum(A[i, ] | A[j, ])
      inter_succ <- sum(A[i, ] & A[j, ])

      # Union and intersection for predecessors
      union_pred <- sum(A[, i] | A[, j])
      inter_pred <- sum(A[, i] & A[, j])

      # Combined Jaccard
      total_union <- union_succ + union_pred
      total_inter <- inter_succ + inter_pred

      if (total_union > 0) {
        scores[i, j] <- total_inter / total_union
      }
    }
  }

  rownames(scores) <- rownames(A)
  colnames(scores) <- colnames(A)
  return(scores)
}

#' Compute Adamic-Adar score
#' @param A Binary adjacency matrix
#' @return Score matrix
#' @keywords internal
.compute_adamic_adar <- function(A) {
  n <- nrow(A)
  scores <- matrix(0, n, n)

  # Compute degrees
  out_degree <- rowSums(A)
  in_degree <- colSums(A)
  total_degree <- out_degree + in_degree

  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      # Shared successors weighted by 1/log(degree)
      shared_succ <- which(A[i, ] == 1 & A[j, ] == 1)
      aa_succ <- 0
      for (k in shared_succ) {
        if (total_degree[k] > 1) {
          aa_succ <- aa_succ + 1 / log(total_degree[k])
        }
      }

      # Shared predecessors weighted by 1/log(degree)
      shared_pred <- which(A[, i] == 1 & A[, j] == 1)
      aa_pred <- 0
      for (k in shared_pred) {
        if (total_degree[k] > 1) {
          aa_pred <- aa_pred + 1 / log(total_degree[k])
        }
      }

      scores[i, j] <- aa_succ + aa_pred
    }
  }

  rownames(scores) <- rownames(A)
  colnames(scores) <- colnames(A)
  return(scores)
}

#' Compute Preferential Attachment score
#' @param A Binary adjacency matrix
#' @return Score matrix
#' @keywords internal
.compute_preferential_attachment <- function(A) {
  n <- nrow(A)

  # Out-degree and in-degree
  out_degree <- rowSums(A)
  in_degree <- colSums(A)

  # PA score: out_degree(i) * in_degree(j)
  scores <- outer(out_degree, in_degree, "*")

  rownames(scores) <- rownames(A)
  colnames(scores) <- colnames(A)
  return(scores)
}

#' Compute Katz index score
#' @param A Binary adjacency matrix
#' @param damping Damping factor (beta)
#' @param max_iter Maximum iterations
#' @param verbose Print progress
#' @return Score matrix
#' @keywords internal
.compute_katz <- function(A, damping = 0.85, max_iter = 100, verbose = TRUE) {
  n <- nrow(A)

  # Check spectral radius to ensure convergence
  # damping should be less than 1/max_eigenvalue
  ev <- eigen(A, only.values = TRUE)$values
  spectral_radius <- max(abs(ev))

  if (spectral_radius > 0) {
    max_damping <- 0.99 / spectral_radius
    if (damping >= max_damping) {
      damping <- max_damping * 0.9
      if (verbose) {
        cat(sprintf("Adjusted damping to %.4f for convergence\n", damping))
      }
    }
  }

  # Iterative computation: S = sum_{l=1}^{inf} beta^l * A^l
  # Equivalent to S = beta*A + beta^2*A^2 + ...
  # Can be computed as S = (I - beta*A)^(-1) - I

  I_mat <- diag(n)
  scores <- NULL

  # Try matrix inversion approach first
  inv_result <- tryCatch({
    inv_mat <- solve(I_mat - damping * A)
    inv_mat - I_mat
  }, error = function(e) {
    NULL
  })

  if (!is.null(inv_result)) {
    scores <- inv_result
  } else {
    # Fall back to iterative approach
    if (verbose) {
      cat("Using iterative Katz computation...\n")
    }

    scores <- matrix(0, n, n)
    A_power <- A
    beta_power <- damping

    for (iter in seq_len(max_iter)) {
      scores <- scores + beta_power * A_power
      A_power <- A_power %*% A
      beta_power <- beta_power * damping

      # Early stopping if contribution is negligible
      if (max(abs(beta_power * A_power)) < 1e-10) {
        break
      }
    }
  }

  rownames(scores) <- rownames(A)
  colnames(scores) <- colnames(A)
  return(scores)
}

#' Compute Hitting Time based score
#' @param P Transition probability matrix (row-stochastic)
#' @param max_iter Maximum iterations
#' @param verbose Print progress
#' @return Score matrix (lower values = stronger predicted link)
#' @keywords internal
.compute_hitting_time <- function(P, max_iter = 100, verbose = TRUE) {
  n <- nrow(P)

  # Ensure P is row-stochastic
  row_sums <- rowSums(P)
  if (any(row_sums == 0)) {
    # Handle absorbing states by adding small self-loop
    zero_rows <- which(row_sums == 0)
    P[zero_rows, zero_rows] <- 1
    row_sums <- rowSums(P)
  }
  P <- P / row_sums

  # Initialize hitting time matrix
  # H[i,j] = expected number of steps to reach j from i
  H <- matrix(Inf, n, n)
  diag(H) <- 0  # Zero steps to reach yourself

  # Iterative computation
  # H[i,j] = 1 + sum_k P[i,k] * H[k,j] for k != j
  # We iterate until convergence

  for (iter in seq_len(max_iter)) {
    H_old <- H

    for (i in seq_len(n)) {
      for (j in seq_len(n)) {
        if (i != j) {
          # Expected hitting time from i to j
          # = 1 + sum over all k of P[i,k] * H[k,j]
          # but for k=j: if we step directly to j, we're done (contributes P[i,j] * 0)

          h_ij <- 1
          for (k in seq_len(n)) {
            if (P[i, k] > 0) {
              if (k == j) {
                # Stepping directly to target contributes 0 additional steps
                h_ij <- h_ij + P[i, k] * 0
              } else {
                h_ij <- h_ij + P[i, k] * H[k, j]
              }
            }
          }
          H[i, j] <- h_ij
        }
      }
    }

    # Check convergence
    max_change <- max(abs(H - H_old)[is.finite(H) & is.finite(H_old)])
    if (!is.na(max_change) && max_change < 1e-6) {
      if (verbose) {
        cat(sprintf("Hitting time converged after %d iterations\n", iter))
      }
      break
    }
  }

  # Convert to similarity score (inverse relationship)
  # Higher score = more likely link
  # Handle Inf values
  max_finite <- max(H[is.finite(H)])
  H[!is.finite(H)] <- max_finite * 2

  # Invert: lower hitting time = higher score
  scores <- max(H) - H

  rownames(scores) <- rownames(P)
  colnames(scores) <- colnames(P)
  return(scores)
}

#' Create predictions data frame from score matrix
#' @param scores Score matrix
#' @param A Adjacency matrix
#' @param states State names
#' @param exclude_existing Exclude existing edges
#' @param top_n Return top N
#' @param threshold Score threshold
#' @param verbose Print progress
#' @return Data frame of predictions
#' @keywords internal
.create_predictions_df <- function(scores, A, states, exclude_existing,
                                    top_n, threshold, verbose) {
  n <- length(states)

  # Create all pairs
  from <- rep(states, each = n)
  to <- rep(states, times = n)
  score <- as.vector(scores)
  exists <- as.vector(A) > 0

  predictions <- data.frame(
    from = from,
    to = to,
    score = score,
    exists = exists,
    stringsAsFactors = FALSE
  )

  # Exclude existing edges if requested
  if (exclude_existing) {
    predictions <- predictions[!predictions$exists, ]
  }

  # Remove self-loops from predictions (usually not meaningful for link prediction)
  predictions <- predictions[predictions$from != predictions$to, ]

  # Sort by score (descending)
  predictions <- predictions[order(predictions$score, decreasing = TRUE), ]

  # Apply threshold or top_n
  if (!is.null(threshold)) {
    predictions <- predictions[predictions$score >= threshold, ]
    if (verbose && nrow(predictions) > 0) {
      cat(sprintf("Found %d links with score >= %.4f\n", nrow(predictions), threshold))
    }
  } else {
    n_pred <- min(top_n, nrow(predictions))
    predictions <- predictions[seq_len(n_pred), ]
    if (verbose && nrow(predictions) > 0) {
      cat(sprintf("Returning top %d predicted links\n", nrow(predictions)))
    }
  }

  # Add rank
  predictions$rank <- seq_len(nrow(predictions))

  # Reset row names
  rownames(predictions) <- NULL

  return(predictions)
}

# ==============================================================================
# S3 METHODS
# ==============================================================================

#' Print Link Prediction Results
#'
#' @param x A link_prediction object
#' @param ... Additional arguments (ignored)
#' @return Invisibly returns x
#' @export
print.link_prediction <- function(x, ...) {
  cat("=== Link Prediction Results ===\n\n")

  cat(sprintf("Method: %s\n", x$method))
  cat(sprintf("Network: %d states, %d edges (density: %.3f)\n",
              x$network_stats$n_states,
              x$network_stats$n_edges,
              x$network_stats$density))
  cat("\n")

  if (nrow(x$predictions) == 0) {
    cat("No predictions found.\n")
  } else {
    cat(sprintf("Top %d Predicted Links:\n", nrow(x$predictions)))
    cat(strrep("-", 45), "\n")

    # Format output
    for (i in seq_len(min(10, nrow(x$predictions)))) {
      pred <- x$predictions[i, ]
      status <- if (pred$exists) " [exists]" else ""
      cat(sprintf("%2d. %s -> %s (score: %.4f)%s\n",
                  pred$rank, pred$from, pred$to, pred$score, status))
    }

    if (nrow(x$predictions) > 10) {
      cat(sprintf("... and %d more predictions\n", nrow(x$predictions) - 10))
    }
  }

  invisible(x)
}

#' Summary of Link Prediction Results
#'
#' @param object A link_prediction object
#' @param ... Additional arguments (ignored)
#' @return Invisibly returns summary statistics
#' @export
summary.link_prediction <- function(object, ...) {
  cat("=== Link Prediction Summary ===\n\n")

  # Network statistics
  cat("Network Statistics:\n")
  cat(sprintf("  States: %d\n", object$network_stats$n_states))
  cat(sprintf("  Existing edges: %d\n", object$network_stats$n_edges))
  cat(sprintf("  Network density: %.4f\n", object$network_stats$density))
  cat(sprintf("  Possible edges: %d\n",
              object$network_stats$n_states^2 - object$network_stats$n_states))
  cat("\n")

  # Method details
  cat("Prediction Method:\n")
  cat(sprintf("  Algorithm: %s\n", object$method))
  cat("\n")

  # Prediction statistics
  if (nrow(object$predictions) > 0) {
    cat("Prediction Statistics:\n")
    cat(sprintf("  Predictions returned: %d\n", nrow(object$predictions)))
    cat(sprintf("  Score range: [%.4f, %.4f]\n",
                min(object$predictions$score),
                max(object$predictions$score)))
    cat(sprintf("  Mean score: %.4f\n", mean(object$predictions$score)))
    cat(sprintf("  Median score: %.4f\n", median(object$predictions$score)))

    # Unique states in predictions
    unique_from <- length(unique(object$predictions$from))
    unique_to <- length(unique(object$predictions$to))
    cat(sprintf("  Unique source states: %d\n", unique_from))
    cat(sprintf("  Unique target states: %d\n", unique_to))
  } else {
    cat("No predictions found.\n")
  }

  # Score matrix statistics
  if (!is.null(object$scores_matrix)) {
    cat("\nScore Matrix Statistics:\n")
    non_diag <- object$scores_matrix[row(object$scores_matrix) != col(object$scores_matrix)]
    cat(sprintf("  Full matrix range: [%.4f, %.4f]\n",
                min(non_diag), max(non_diag)))
    cat(sprintf("  Full matrix mean: %.4f\n", mean(non_diag)))
  }

  invisible(list(
    network = object$network_stats,
    method = object$method,
    n_predictions = nrow(object$predictions),
    score_stats = if (nrow(object$predictions) > 0) {
      list(
        min = min(object$predictions$score),
        max = max(object$predictions$score),
        mean = mean(object$predictions$score),
        median = median(object$predictions$score)
      )
    } else NULL
  ))
}

#' Plot Link Prediction Results
#'
#' Visualizes predicted links on the network. Existing edges are shown in gray,
#' predicted edges in red with varying intensity based on prediction score.
#'
#' @param x A link_prediction object
#' @param max_predictions Maximum number of predictions to show (default: 10)
#' @param show_existing Logical. Show existing edges in background (default: TRUE)
#' @param edge_width_range Numeric vector of length 2. Range for edge widths (default: c(0.5, 3))
#' @param layout Character. Layout algorithm: "circle" (default) or "fruchterman"
#' @param vertex_size Numeric. Size of vertices (default: 15)
#' @param vertex_color Character. Color of vertices (default: "lightblue")
#' @param existing_color Character. Color for existing edges (default: "gray80")
#' @param predicted_color Character. Color for predicted edges (default: "red")
#' @param ... Additional arguments passed to plot
#' @return Invisibly returns the plot coordinates
#' @export
plot.link_prediction <- function(x,
                                  max_predictions = 10,
                                  show_existing = TRUE,
                                  edge_width_range = c(0.5, 3),
                                  layout = c("circle", "fruchterman"),
                                  vertex_size = 15,
                                  vertex_color = "lightblue",
                                  existing_color = "gray80",
                                  predicted_color = "red",
                                  ...) {

  layout <- match.arg(layout)

  # Get states
  states <- rownames(x$adjacency)
  n <- length(states)

  # Compute layout
  if (layout == "circle") {
    angles <- seq(0, 2 * pi, length.out = n + 1)[1:n]
    coords <- cbind(cos(angles), sin(angles))
  } else {
    # Fruchterman-Reingold layout (simplified)
    coords <- .fruchterman_reingold_layout(x$adjacency)
  }

  rownames(coords) <- states

  # Setup plot
  old_par <- par(mar = c(1, 1, 3, 1))
  on.exit(par(old_par), add = TRUE)

  plot(NULL, xlim = c(-1.5, 1.5), ylim = c(-1.5, 1.5),
       xlab = "", ylab = "", axes = FALSE, asp = 1,
       main = sprintf("Link Predictions (%s)", x$method), ...)

  # Draw existing edges (background)
  if (show_existing) {
    A <- x$adjacency
    for (i in seq_len(n)) {
      for (j in seq_len(n)) {
        if (A[i, j] > 0 && i != j) {
          .draw_curved_edge(coords[i, ], coords[j, ],
                           color = existing_color, width = 0.5)
        }
      }
    }
  }

  # Draw predicted edges
  preds <- x$predictions[seq_len(min(max_predictions, nrow(x$predictions))), ]

  if (nrow(preds) > 0) {
    # Normalize scores for edge width
    score_range <- range(preds$score)
    if (score_range[2] > score_range[1]) {
      normalized <- (preds$score - score_range[1]) / (score_range[2] - score_range[1])
    } else {
      normalized <- rep(1, nrow(preds))
    }

    widths <- edge_width_range[1] + normalized * (edge_width_range[2] - edge_width_range[1])

    # Color gradient based on score
    for (i in seq_len(nrow(preds))) {
      from_idx <- which(states == preds$from[i])
      to_idx <- which(states == preds$to[i])

      if (length(from_idx) > 0 && length(to_idx) > 0) {
        # Color intensity based on rank
        alpha <- 1 - (i - 1) / max(nrow(preds), 1) * 0.7
        color <- adjustcolor(predicted_color, alpha.f = alpha)

        .draw_curved_edge(coords[from_idx, ], coords[to_idx, ],
                         color = color, width = widths[i], curved = 0.2)
      }
    }
  }

  # Draw vertices
  symbols(coords[, 1], coords[, 2],
          circles = rep(vertex_size / 100, n),
          inches = FALSE, add = TRUE,
          bg = vertex_color, fg = "black")

  # Add labels
  text(coords[, 1], coords[, 2], labels = states, cex = 0.8)

  # Add legend
  legend("bottomright",
         legend = c("Existing edge", "Predicted link"),
         col = c(existing_color, predicted_color),
         lwd = c(1, 2),
         cex = 0.7,
         bty = "n")

  invisible(coords)
}

#' Draw a curved edge between two points
#' @param from Coordinates of start point
#' @param to Coordinates of end point
#' @param color Edge color
#' @param width Edge width
#' @param curved Curvature amount
#' @keywords internal
.draw_curved_edge <- function(from, to, color = "gray", width = 1, curved = 0.1) {
  # Simple curved line using quadratic bezier approximation
  mid <- (from + to) / 2

  # Perpendicular offset for curve
  diff <- to - from
  perp <- c(-diff[2], diff[1])
  perp <- perp / sqrt(sum(perp^2)) * curved

  control <- mid + perp

  # Draw using lines (approximation with multiple segments)
  t_seq <- seq(0, 1, length.out = 20)
  x_curve <- (1 - t_seq)^2 * from[1] + 2 * (1 - t_seq) * t_seq * control[1] + t_seq^2 * to[1]
  y_curve <- (1 - t_seq)^2 * from[2] + 2 * (1 - t_seq) * t_seq * control[2] + t_seq^2 * to[2]

  lines(x_curve, y_curve, col = color, lwd = width)

  # Draw arrowhead at end
  arrow_len <- 0.05
  angle <- atan2(to[2] - y_curve[length(y_curve) - 1],
                 to[1] - x_curve[length(x_curve) - 1])

  arrows_x <- c(to[1] - arrow_len * cos(angle - pi/6),
                to[1] - arrow_len * cos(angle + pi/6))
  arrows_y <- c(to[2] - arrow_len * sin(angle - pi/6),
                to[2] - arrow_len * sin(angle + pi/6))

  lines(c(to[1], arrows_x[1]), c(to[2], arrows_y[1]), col = color, lwd = width)
  lines(c(to[1], arrows_x[2]), c(to[2], arrows_y[2]), col = color, lwd = width)
}

#' Simple Fruchterman-Reingold layout
#' @param A Adjacency matrix
#' @param iterations Number of iterations
#' @return Coordinate matrix
#' @keywords internal
.fruchterman_reingold_layout <- function(A, iterations = 50) {
  n <- nrow(A)

  # Initialize with random positions
  set.seed(42)  # For reproducibility
  coords <- matrix(runif(n * 2, -1, 1), ncol = 2)

  # Layout parameters
  area <- 4
  k <- sqrt(area / n)  # Optimal distance
  temp <- 1

  for (iter in seq_len(iterations)) {
    # Calculate repulsive forces
    disp <- matrix(0, n, 2)

    for (i in seq_len(n)) {
      for (j in seq_len(n)) {
        if (i != j) {
          delta <- coords[i, ] - coords[j, ]
          dist <- sqrt(sum(delta^2))
          if (dist > 0) {
            # Repulsive force
            disp[i, ] <- disp[i, ] + (delta / dist) * (k^2 / dist)
          }
        }
      }
    }

    # Calculate attractive forces
    for (i in seq_len(n)) {
      for (j in seq_len(n)) {
        if (A[i, j] > 0 || A[j, i] > 0) {
          delta <- coords[i, ] - coords[j, ]
          dist <- sqrt(sum(delta^2))
          if (dist > 0) {
            # Attractive force
            disp[i, ] <- disp[i, ] - (delta / dist) * (dist^2 / k)
          }
        }
      }
    }

    # Update positions
    for (i in seq_len(n)) {
      disp_len <- sqrt(sum(disp[i, ]^2))
      if (disp_len > 0) {
        coords[i, ] <- coords[i, ] + (disp[i, ] / disp_len) * min(disp_len, temp)
      }
    }

    # Cool down
    temp <- temp * 0.95

    # Bound positions
    coords[, 1] <- pmax(pmin(coords[, 1], 1), -1)
    coords[, 2] <- pmax(pmin(coords[, 2], 1), -1)
  }

  return(coords)
}

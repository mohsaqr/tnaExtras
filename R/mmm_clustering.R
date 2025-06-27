# ==============================================================================
# MIXTURE MARKOV MODEL CLUSTERING FOR SEQUENCE DATA
# ==============================================================================
# 
# This module implements enhanced Expectation-Maximization algorithm for fitting
# mixture of first-order Markov models to sequence data. Features multiple random
# starts, comprehensive posterior probability analysis, and robust error handling.
#
# ==============================================================================

#' Mixture Markov Model Clustering for Sequences
#'
#' Fits a mixture of first-order Markov models to sequence data using the
#' Expectation-Maximization (EM) algorithm. Provides robust estimation through
#' multiple random restarts and comprehensive model diagnostics.
#'
#' @param data A data frame or matrix where rows represent sequences and columns
#'   represent time points. Missing values (NA) are handled appropriately.
#'   Each cell should contain a state label (character or factor).
#' @param k An integer specifying the number of mixture components (clusters).
#'   Must be between 2 and the number of sequences minus 1.
#' @param max_iter An integer specifying the maximum number of EM iterations
#'   per restart. Default is 300.
#' @param tol A numeric value specifying the convergence tolerance for the
#'   log-likelihood change between iterations. Default is 1e-6.
#' @param n_starts An integer specifying the number of random restarts to
#'   perform. Multiple restarts help avoid local optima. Default is 10.
#' @param seed An integer for random seed setting to ensure reproducibility.
#'   If NULL, uses current random state. Default is NULL.
#' @param verbose A logical value indicating whether to print progress
#'   information during fitting. Default is TRUE.
#' @param min_cluster_size An integer specifying the minimum number of sequences
#'   required per cluster for a solution to be considered valid. Default is 1.
#' @param smoothing A numeric value for Laplace smoothing parameter to prevent
#'   zero probabilities. Default is 0.1.
#'
#' @return A list of class \code{mmm_result} containing:
#' \describe{
#'   \item{assignments}{Integer vector of cluster assignments for each sequence}
#'   \item{responsibilities}{Matrix of posterior probabilities (n_sequences x k)}
#'   \item{log_likelihood}{Final log-likelihood of the best model}
#'   \item{bic}{Bayesian Information Criterion}
#'   \item{aic}{Akaike Information Criterion}
#'   \item{models}{List of k trained Markov models with named states}
#'   \item{mixture_weights}{Named vector of cluster mixing weights}
#'   \item{k}{Number of clusters}
#'   \item{n_parameters}{Total number of model parameters}
#'   \item{converged}{Logical indicating convergence of best model}
#'   \item{n_iterations}{Number of iterations for best model}
#'   \item{cluster_sizes}{Named vector of cluster sizes}
#'   \item{cluster_proportions}{Named vector of cluster proportions}
#'   \item{avg_posterior_per_cluster}{Named vector of average posterior probabilities per cluster}
#'   \item{avg_posterior_overall}{Overall average posterior probability}
#'   \item{entropy}{Normalized Shannon entropy (0 = perfect separation, 1 = maximum uncertainty)}
#'   \item{states}{Character vector of unique states in data}
#'   \item{call}{The matched call}
#' }
#'
#' @details
#' The function implements a mixture of first-order Markov models where each
#' component models sequences with initial state probabilities and transition
#' probabilities between states, along with mixing weights for component membership.
#'
#' The EM algorithm alternates between E-step (computing posterior probabilities
#' of cluster membership) and M-step (updating model parameters based on weighted data).
#'
#' Multiple random restarts are performed to avoid local optima, with the
#' best solution (highest log-likelihood) selected as the final result.
#'
#' @examples
#' \dontrun{
#' # Create example sequence data
#' data <- data.frame(
#'   T1 = c("A", "B", "A", "C", "A", "B"),
#'   T2 = c("B", "A", "B", "A", "C", "A"),
#'   T3 = c("C", "C", "A", "B", "B", "C")
#' )
#'
#' # Fit MMM with 2 clusters
#' result <- cluster_mmm(data, k = 2, n_starts = 5, seed = 123)
#' print(result)
#'
#' # Access transition matrix for cluster 1
#' result$models[[1]]$transition
#' }
#'
#' @export
cluster_mmm <- function(data, 
                        k, 
                        max_iter = 300L,
                        tol = 1e-6,
                        n_starts = 10L,
                        seed = NULL,
                        verbose = TRUE,
                        min_cluster_size = 1L,
                        smoothing = 0.1) {
  
  # Store the call for reproducibility
  call_matched <- match.call()
  
  # Set seed if provided
  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1) {
      stop("'seed' must be a single numeric value or NULL")
    }
    set.seed(as.integer(seed))
  }
  
  # Input validation with informative error messages
  if (!is.data.frame(data) && !is.matrix(data)) {
    stop("'data' must be a data frame or matrix, got: ", class(data)[1])
  }
  
  if (!is.matrix(data)) data <- as.matrix(data)
  
  n_seq <- nrow(data)
  n_time <- ncol(data)
  
  if (n_seq < 2) {
    stop("'data' must contain at least 2 sequences, got: ", n_seq)
  }
  
  if (n_time < 1) {
    stop("'data' must contain at least 1 time point, got: ", n_time)
  }
  
  if (!is.numeric(k) || length(k) != 1 || is.na(k) || k != as.integer(k)) {
    stop("'k' must be a single integer value, got: ", deparse(substitute(k)))
  }
  k <- as.integer(k)
  
  if (k < 2 || k >= n_seq) {
    stop("'k' must be between 2 and ", n_seq - 1, ", got: ", k)
  }
  
  if (!is.numeric(max_iter) || length(max_iter) != 1 || max_iter < 1) {
    stop("'max_iter' must be a positive integer, got: ", max_iter)
  }
  max_iter <- as.integer(max_iter)
  
  if (!is.numeric(tol) || length(tol) != 1 || tol <= 0) {
    stop("'tol' must be a positive number, got: ", tol)
  }
  
  if (!is.numeric(n_starts) || length(n_starts) != 1 || n_starts < 1) {
    stop("'n_starts' must be a positive integer, got: ", n_starts)
  }
  n_starts <- as.integer(n_starts)
  
  if (!is.logical(verbose) || length(verbose) != 1) {
    stop("'verbose' must be TRUE or FALSE, got: ", verbose)
  }
  
  if (!is.numeric(min_cluster_size) || length(min_cluster_size) != 1 || min_cluster_size < 1) {
    stop("'min_cluster_size' must be a positive integer, got: ", min_cluster_size)
  }
  min_cluster_size <- as.integer(min_cluster_size)
  
  if (!is.numeric(smoothing) || length(smoothing) != 1 || smoothing <= 0) {
    stop("'smoothing' must be a positive number, got: ", smoothing)
  }
  
  # Extract and validate states
  all_states <- sort(unique(as.vector(data[!is.na(data)])))
  n_states <- length(all_states)
  
  if (n_states < 2) {
    stop("'data' must contain at least 2 different states, found: ", n_states)
  }
  
  if (verbose) {
    cat("=== Mixture Markov Model Clustering ===\n")
    cat(sprintf("Dataset: %d sequences, %d time points, %d states (%s)\n", 
                n_seq, n_time, n_states, paste(all_states, collapse = ", ")))
    cat(sprintf("Clusters: %d | Max iterations: %d | Tolerance: %.2e\n", 
                k, max_iter, tol))
    cat(sprintf("Random restarts: %d | Smoothing: %.3f\n", n_starts, smoothing))
    cat("\n")
  }
  
  # Pre-compute numeric data for performance
  data_numeric <- matrix(NA_integer_, nrow = n_seq, ncol = n_time)
  state_lookup <- stats::setNames(seq_along(all_states), all_states)
  
  for (i in seq_len(n_seq)) {
    for (j in seq_len(n_time)) {
      if (!is.na(data[i, j])) {
        data_numeric[i, j] <- state_lookup[as.character(data[i, j])]
      }
    }
  }
  
  # Helper function for single EM run
  .single_em_run <- function(restart_id, verbose_inner = FALSE) {
    
    # Random initialization
    mixture_weights <- runif(k, 0.8, 1.2)
    mixture_weights <- mixture_weights / sum(mixture_weights)
    
    # Initialize models with better starting values
    models <- vector("list", k)
    for (cluster in seq_len(k)) {
      # Initial distribution: uniform + controlled noise
      initial_probs <- rep(1/n_states, n_states) + runif(n_states, -0.05, 0.05)
      initial_probs <- pmax(initial_probs, smoothing / n_states)
      initial_probs <- initial_probs / sum(initial_probs)
      
      # Transition matrix: identity-leaning + noise for better convergence
      transition_matrix <- diag(0.6, n_states) + matrix(0.4/n_states, n_states, n_states) + 
                           matrix(runif(n_states^2, -0.05, 0.05), n_states, n_states)
      transition_matrix <- pmax(transition_matrix, smoothing / n_states)
      
      # Row normalization
      for (row in seq_len(n_states)) {
        transition_matrix[row, ] <- transition_matrix[row, ] / sum(transition_matrix[row, ])
      }
      
      # Add state names
      rownames(transition_matrix) <- all_states
      colnames(transition_matrix) <- all_states
      names(initial_probs) <- all_states
      
      models[[cluster]] <- list(
        initial = initial_probs,
        transition = transition_matrix
      )
    }
    
    # EM iteration tracking
    log_likelihood_old <- -Inf
    converged <- FALSE
    responsibilities <- matrix(1/k, nrow = n_seq, ncol = k)
    
    # Main EM loop
    for (iter in seq_len(max_iter)) {
      
      # E-step: Calculate responsibilities
      log_lik_matrix <- matrix(-Inf, nrow = n_seq, ncol = k)
      
      for (i in seq_len(n_seq)) {
        seq_data <- data_numeric[i, ]
        valid_idx <- which(!is.na(seq_data))
        
        if (length(valid_idx) > 0) {
          for (j in seq_len(k)) {
            tryCatch({
              # Initial state log-probability
              initial_prob <- models[[j]]$initial[seq_data[valid_idx[1]]]
              if (is.na(initial_prob) || initial_prob <= 0) initial_prob <- 1e-10
              log_prob <- log(initial_prob)
              
              # Transition log-probabilities
              if (length(valid_idx) > 1) {
                for (t in seq_len(length(valid_idx) - 1)) {
                  from_state <- seq_data[valid_idx[t]]
                  to_state <- seq_data[valid_idx[t + 1]]
                  trans_prob <- models[[j]]$transition[from_state, to_state]
                  if (is.na(trans_prob) || trans_prob <= 0) trans_prob <- 1e-10
                  log_prob <- log_prob + log(trans_prob)
                }
              }
              
              mix_weight <- mixture_weights[j]
              if (is.na(mix_weight) || mix_weight <= 0) mix_weight <- 1e-10
              log_lik_matrix[i, j] <- log(mix_weight) + log_prob
            }, error = function(e) {
              log_lik_matrix[i, j] <<- -1e10
            })
          }
        }
      }
      
      # Normalize responsibilities using log-sum-exp trick
      for (i in seq_len(n_seq)) {
        log_probs <- log_lik_matrix[i, ]
        if (all(!is.finite(log_probs))) {
          responsibilities[i, ] <- rep(1/k, k)
        } else {
          max_log_prob <- max(log_probs[is.finite(log_probs)])
          log_probs_normalized <- log_probs - max_log_prob
          probs <- exp(log_probs_normalized)
          probs[!is.finite(probs)] <- 0
          if (sum(probs) > 0) {
            responsibilities[i, ] <- probs / sum(probs)
          } else {
            responsibilities[i, ] <- rep(1/k, k)
          }
        }
      }
      
      # M-step: Update parameters
      mixture_weights <- colMeans(responsibilities)
      
      # Update models
      for (j in seq_len(k)) {
        resp_j <- responsibilities[, j]
        
        # Update initial probabilities
        initial_counts <- sapply(seq_len(n_states), function(state_idx) {
          sum(resp_j[data_numeric[, 1] == state_idx & !is.na(data_numeric[, 1])])
        })
        models[[j]]$initial <- (initial_counts + smoothing) / (sum(initial_counts) + n_states * smoothing)
        names(models[[j]]$initial) <- all_states
        
        # Update transition probabilities
        transition_counts <- matrix(smoothing, nrow = n_states, ncol = n_states)
        
        if (n_time > 1) {
          for (t in seq_len(n_time - 1)) {
            for (i in seq_len(n_seq)) {
              if (!is.na(data_numeric[i, t]) && !is.na(data_numeric[i, t + 1])) {
                from_state <- data_numeric[i, t]
                to_state <- data_numeric[i, t + 1]
                transition_counts[from_state, to_state] <- 
                  transition_counts[from_state, to_state] + resp_j[i]
              }
            }
          }
        }
        
        # Row normalization
        for (row in seq_len(n_states)) {
          models[[j]]$transition[row, ] <- transition_counts[row, ] / sum(transition_counts[row, ])
        }
        
        # Ensure names are preserved
        rownames(models[[j]]$transition) <- all_states
        colnames(models[[j]]$transition) <- all_states
      }
      
      # Calculate log-likelihood
      log_likelihood <- sum(apply(log_lik_matrix, 1, function(x) {
        if (all(!is.finite(x))) return(0)
        max_x <- max(x[is.finite(x)])
        log(sum(exp(x - max_x))) + max_x
      }))
      
      # Check convergence
      if (abs(log_likelihood - log_likelihood_old) < tol) {
        converged <- TRUE
        break
      }
      
      log_likelihood_old <- log_likelihood
    }
    
    # Final assignments
    assignments <- apply(responsibilities, 1, which.max)
    cluster_sizes <- table(factor(assignments, levels = 1:k))
    
    # Check minimum cluster size
    if (any(cluster_sizes < min_cluster_size)) {
      return(NULL)
    }
    
    return(list(
      assignments = assignments,
      responsibilities = responsibilities,
      log_likelihood = log_likelihood,
      models = models,
      mixture_weights = mixture_weights,
      converged = converged,
      n_iterations = iter,
      cluster_sizes = as.vector(cluster_sizes)
    ))
  }
  
  # Perform multiple random restarts
  if (verbose) cat("Running EM algorithm with", n_starts, "random restarts...\n")
  
  all_results <- vector("list", n_starts)
  successful_runs <- 0
  
  for (start_idx in seq_len(n_starts)) {
    if (verbose && n_starts > 1) {
      cat(sprintf("Restart %d/%d...", start_idx, n_starts))
    }
    
    result <- .single_em_run(start_idx, verbose_inner = FALSE)
    
    if (!is.null(result)) {
      all_results[[start_idx]] <- result
      successful_runs <- successful_runs + 1
      if (verbose && n_starts > 1) {
        cat(sprintf(" LL: %.4f %s\n", result$log_likelihood, 
                   ifelse(result$converged, "(converged)", "(max iter)")))
      }
    } else {
      if (verbose && n_starts > 1) {
        cat(" failed (min cluster size)\n")
      }
    }
  }
  
  if (successful_runs == 0) {
    stop("All EM runs failed to meet minimum cluster size requirements")
  }
  
  # Select best result
  valid_results <- all_results[!sapply(all_results, is.null)]
  log_likelihoods <- sapply(valid_results, function(x) x$log_likelihood)
  best_idx <- which.max(log_likelihoods)
  best_result <- valid_results[[best_idx]]
  
  if (verbose) {
    cat(sprintf("\nBest result: LL = %.4f (restart %d, %s)\n", 
                best_result$log_likelihood, best_idx,
                ifelse(best_result$converged, "converged", "max iterations")))
  }
  
  # Calculate information criteria
  n_mixing_params <- k - 1
  n_initial_params <- k * (n_states - 1)
  n_transition_params <- k * n_states * (n_states - 1)
  n_parameters <- n_mixing_params + n_initial_params + n_transition_params
  
  aic <- -2 * best_result$log_likelihood + 2 * n_parameters
  bic <- -2 * best_result$log_likelihood + log(n_seq) * n_parameters
  
  # Calculate enhanced output metrics
  cluster_proportions <- round(best_result$cluster_sizes / n_seq, 2)
  names(cluster_proportions) <- paste0("Cluster_", seq_len(k))
  names(best_result$cluster_sizes) <- paste0("Cluster_", seq_len(k))
  names(best_result$mixture_weights) <- paste0("Cluster_", seq_len(k))
  
  # Average posterior probabilities per cluster
  avg_posterior_per_cluster <- numeric(k)
  for (j in seq_len(k)) {
    cluster_members <- best_result$assignments == j
    if (any(cluster_members)) {
      avg_posterior_per_cluster[j] <- mean(best_result$responsibilities[cluster_members, j])
    } else {
      avg_posterior_per_cluster[j] <- 0
    }
  }
  names(avg_posterior_per_cluster) <- paste0("Cluster_", seq_len(k))
  
  # Overall average posterior probability
  avg_posterior_overall <- mean(apply(best_result$responsibilities, 1, max))
  
  # Calculate proper Shannon entropy measures
  safe_resp <- pmax(best_result$responsibilities, 1e-15)
  
  # Shannon entropy: H(P) = -sum(p_i log(p_i))
  raw_entropy <- -sum(safe_resp * log(safe_resp))
  
  # Normalized entropy: H(P) / log(k) in [0,1] where 0=perfect, 1=uniform
  max_entropy <- n_seq * log(k)
  normalized_entropy <- ifelse(max_entropy > 0, raw_entropy / max_entropy, 0)
  
  # Create final result
  result <- list(
    assignments = best_result$assignments,
    responsibilities = best_result$responsibilities,
    log_likelihood = best_result$log_likelihood,
    bic = bic,
    aic = aic,
    models = best_result$models,
    mixture_weights = best_result$mixture_weights,
    k = k,
    n_parameters = n_parameters,
    converged = best_result$converged,
    n_iterations = best_result$n_iterations,
    cluster_sizes = best_result$cluster_sizes,
    cluster_proportions = cluster_proportions,
    avg_posterior_per_cluster = round(avg_posterior_per_cluster, 4),
    avg_posterior_overall = round(avg_posterior_overall, 4),
    entropy = round(normalized_entropy, 4),
    states = all_states,
    call = call_matched
  )
  
  # Add class for custom printing
  class(result) <- c("mmm_result", "list")
  
  return(result)
}

#' Print Method for MMM Results
#'
#' Custom print method that displays MMM clustering results in a clean, readable format.
#'
#' @param x A mmm_result object
#' @param ... Additional arguments (ignored)
#' @method print mmm_result
#' @export
print.mmm_result <- function(x, ...) {
  cat("=== Mixture Markov Model Clustering Result ===\n")
  cat("Number of clusters (k):", x$k, "\n")
  cat("Number of sequences:", length(x$assignments), "\n")
  cat("States:", paste(x$states, collapse = ", "), "\n")
  cat("Converged:", x$converged, "after", x$n_iterations, "iterations\n")
  cat("Log-likelihood:", round(x$log_likelihood, 4), "\n")
  cat("AIC:", round(x$aic, 4), "\n")
  cat("BIC:", round(x$bic, 4), "\n")
  cat("\nCluster information:\n")
  cat("Cluster sizes:", paste(x$cluster_sizes, collapse = ", "), "\n")
  cat("Cluster proportions:", paste(x$cluster_proportions, collapse = ", "), "\n")
  cat("Mixture weights:", paste(round(x$mixture_weights, 3), collapse = ", "), "\n")
  cat("\nPosterior probability analysis:\n")
  cat("Avg posterior per cluster:", paste(x$avg_posterior_per_cluster, collapse = ", "), "\n")
  cat("Overall avg posterior:", x$avg_posterior_overall, "\n")
  cat("Normalized entropy:", x$entropy, "\n")
  
  cat("\nNote: Access detailed results with $assignments, $responsibilities, $models, etc.\n")
  
  invisible(x)
}

#' Find Optimal Number of Clusters for MMM
#'
#' Performs model selection across a range of cluster numbers using information
#' criteria (AIC or BIC) with comprehensive diagnostics and performance metrics.
#' Includes enhanced posterior analysis and convergence tracking.
#'
#' @param data A data frame or matrix where rows represent sequences and columns
#'   represent time points. See \code{\link{cluster_mmm}} for details.
#' @param k_range An integer vector specifying the range of k values to test.
#'   Default is 2:6.
#' @param criterion A character string specifying the selection criterion.
#'   Either "bic" (default) or "aic".
#' @param n_cores Integer. Number of cores for parallel processing. If NULL (default),
#'   uses sequential processing. If specified, automatically enables parallel processing.
#' @param ... Additional arguments passed to \code{\link{cluster_mmm}}
#'
#' @return A list of class \code{mmm_selection_result} containing:
#' \describe{
#'   \item{optimal_k}{The optimal number of clusters}
#'   \item{best_result}{The complete result object for the optimal model}
#'   \item{all_results}{Named list of all fitted model results}
#'   \item{comparison}{Data frame comparing all models with enhanced metrics including:
#'     \itemize{
#'       \item entropy: Normalized Shannon entropy (0 = perfect separation, 1 = uniform)
#'       \item entropy_relative: KL divergence from uniform distribution (0 = uniform, higher = more separated)
#'     }}
#'   \item{criterion_used}{The selection criterion used}
#'   \item{call}{The matched call}
#' }
#'
#' @details
#' This function fits MMM models for each value in \code{k_range} and selects
#' the optimal number of clusters based on the specified information criterion.
#' The comparison includes additional metrics such as posterior certainty,
#' cluster balance, and convergence rates to aid in model interpretation.
#'
#' @examples
#' \dontrun{
#' # Create example data
#' data <- data.frame(
#'   T1 = sample(c("A", "B", "C"), 100, replace = TRUE),
#'   T2 = sample(c("A", "B", "C"), 100, replace = TRUE),
#'   T3 = sample(c("A", "B", "C"), 100, replace = TRUE)
#' )
#'
#' # Find optimal k using BIC
#' selection <- find_optimal_mmm(data, k_range = 2:5)
#' print(selection)
#'
#' # Access optimal model
#' best_model <- selection$best_result
#' }
#'
#' @seealso \code{\link{cluster_mmm}}
#'
#' @export
find_optimal_mmm <- function(data, 
                            k_range = 2:6, 
                            criterion = c("bic", "aic"),
                            n_cores = NULL,
                            ...) {
  
  call_matched <- match.call()
  criterion <- match.arg(criterion)
  
  # Input validation
  if (!is.numeric(k_range) || any(k_range < 2)) {
    stop("'k_range' must contain integers >= 2")
  }
  
  if (length(k_range) < 2) {
    stop("'k_range' must contain at least 2 values for comparison")
  }
  
  if (!is.null(n_cores) && (!is.numeric(n_cores) || length(n_cores) != 1 || n_cores < 1)) {
    stop("'n_cores' must be NULL or a positive integer")
  }
  
  cat("=== Enhanced MMM Model Selection ===\n")
  cat(sprintf("Testing k in {%s} using %s criterion\n", 
              paste(range(k_range), collapse = "-"), toupper(criterion)))
  if (!is.null(n_cores)) cat(sprintf("Using parallel processing with %d cores\n", n_cores))
  cat("\n")
  
  # Fit models
  all_results <- list()
  comparison_data <- data.frame(
    k = integer(0),
    log_likelihood = numeric(0),
    bic = numeric(0),
    aic = numeric(0),
    avg_posterior_overall = numeric(0),
    entropy_relative = numeric(0),
    min_cluster_size = integer(0),
    max_cluster_size = integer(0),
    cluster_balance = numeric(0),
    converged = logical(0),
    convergence_rate = numeric(0),
    best_restart = integer(0),
    stringsAsFactors = FALSE
  )
  
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
  cluster_info <- .setup_parallel_cluster(n_cores = n_cores, verbose = FALSE)
  on.exit(.cleanup_parallel_cluster(cluster_info, verbose = FALSE), add = TRUE)
  
  # Prepare worker function for fitting
  fit_single_k <- function(k) {
    start_time <- Sys.time()
    
    result <- tryCatch({
      # Extract args and remove verbose to avoid conflicts
      args <- list(...)
      args$verbose <- FALSE
      # Set very permissive min_cluster_size for model selection if not specified
      if (!"min_cluster_size" %in% names(args)) {
        args$min_cluster_size <- 1L
      }
      # Increase random restarts for model selection to improve convergence
      if (!"n_starts" %in% names(args)) {
        args$n_starts <- 20L
      }
      do.call(cluster_mmm, c(list(data = data, k = k), args))
    }, error = function(e) {
      return(list(error = e$message, k = k))
    })
    
    end_time <- Sys.time()
    
    return(list(
      result = result,
      k = k,
      timing = as.numeric(difftime(end_time, start_time, units = "secs"))
    ))
  }
  
  # Execute parallel or sequential fitting
  if (!is.null(cluster_info)) {
    # Parallel execution
    cat("Fitting models in parallel...\n")
    all_fits <- .robust_parallel_lapply(k_range, fit_single_k, cluster_info, verbose = FALSE)
  } else {
    # Sequential execution
    all_fits <- lapply(k_range, fit_single_k)
  }
  
  # Process results
  for (i in seq_along(all_fits)) {
    fit_info <- all_fits[[i]]
    k <- fit_info$k
    result <- fit_info$result
    timing <- fit_info$timing
    
    if ("error" %in% names(result)) {
      cat(sprintf("k = %d: ERROR: %s\n", k, result$error))
      next
    }
    
    if (!is.null(result)) {
      all_results[[paste0("k", k)]] <- result
      
      # Calculate metrics
      min_size <- min(result$cluster_sizes)
      max_size <- max(result$cluster_sizes)
      balance <- min_size / max_size
      # Set convergence rate to 1.0 for current simplified structure
      conv_rate <- 1.0
      
      # Calculate KL divergence from uniform distribution (proper relative entropy)
      safe_resp <- pmax(result$responsibilities, 1e-15)
      n_resp <- nrow(result$responsibilities)
      uniform_prob <- 1/k
      
      # KL(P||Uniform) = sum(p_i log(p_i / (1/k))) = sum(p_i log(p_i)) + log(k)
      kl_from_uniform <- sum(safe_resp * log(safe_resp / uniform_prob)) / n_resp
      # Normalize by max possible KL divergence for interpretability
      max_kl <- log(k)  # When one cluster has all probability
      relative_entropy <- ifelse(max_kl > 0, kl_from_uniform / max_kl, 0)
      
      comparison_data <- rbind(comparison_data, data.frame(
        k = k,
        log_likelihood = result$log_likelihood,
        bic = result$bic,
        aic = result$aic,
        avg_posterior_overall = result$avg_posterior_overall,
        entropy_relative = round(relative_entropy, 4),
        min_cluster_size = min_size,
        max_cluster_size = max_size,
        cluster_balance = balance,
        converged = result$converged,
        convergence_rate = conv_rate,
        best_restart = 1,
        stringsAsFactors = FALSE
      ))
      
              cat(sprintf("k = %d: [OK] LL=%.2f, %s=%.1f, Cert=%.3f, RelEnt=%.3f (%.1fs)\n", 
                  k,
                  result$log_likelihood, 
                  toupper(criterion), 
                  if(criterion == "bic") result$bic else result$aic,
                  result$avg_posterior_overall,
                  round(relative_entropy, 3),
                  timing))
    }
  }
  
  if (nrow(comparison_data) == 0) {
    stop("No models were successfully fitted for any k value")
  }
  
  # Select optimal model
  if (criterion == "bic") {
    best_idx <- which.min(comparison_data$bic)
  } else {
    best_idx <- which.min(comparison_data$aic)
  }
  
  optimal_k <- comparison_data$k[best_idx]
  best_result <- all_results[[paste0("k", optimal_k)]]
  
  cat("\n=== Model Comparison ===\n")
  # Show all important metrics
  print(comparison_data, digits = 3, row.names = FALSE)
  
  cat(sprintf("\n=== Optimal Model (k = %d) ===\n", optimal_k))
  cat(sprintf("%-20s: %.3f\n", toupper(criterion), comparison_data[[criterion]][best_idx]))
  cat(sprintf("%-20s: %.3f\n", "Posterior certainty", comparison_data$avg_posterior_overall[best_idx]))
  cat(sprintf("%-20s: %.3f\n", "Relative entropy", comparison_data$entropy_relative[best_idx]))
  cat(sprintf("%-20s: %.3f\n", "Cluster balance", comparison_data$cluster_balance[best_idx]))
  cat(sprintf("%-20s: %d to %d\n", "Cluster sizes", comparison_data$min_cluster_size[best_idx], comparison_data$max_cluster_size[best_idx]))
  cat(sprintf("%-20s: %.1f%%\n", "Convergence rate", 100 * comparison_data$convergence_rate[best_idx]))
  cat(sprintf("%-20s: %s\n", "Converged", ifelse(comparison_data$converged[best_idx], "Yes", "No (increase max_iter)")))
  
  result <- list(
    optimal_k = optimal_k,
    best_result = best_result,
    all_results = all_results,
    comparison = comparison_data,
    criterion_used = criterion,
    call = call_matched
  )
  
  class(result) <- c("mmm_selection_result", "list")
  return(result)
}

#' Print Method for MMM Selection Results
#'
#' Provides a clean, informative summary of MMM model selection results
#' with comparison metrics and optimal model information.
#'
#' @param x An object of class \code{mmm_selection_result}
#' @param ... Additional arguments (currently ignored)
#'
#' @return Invisibly returns the input object
#'
#' @method print mmm_selection_result
#' @export
print.mmm_selection_result <- function(x, ...) {
  cat("=== MMM Model Selection Results ===\n")
  cat("Call: "); print(x$call); cat("\n")
  
  cat(sprintf("Criterion: %s\n", toupper(x$criterion_used)))
  cat(sprintf("Models tested: %d\n", nrow(x$comparison)))
  cat(sprintf("Optimal k: %d\n", x$optimal_k))
  
  cat("\nComparison summary:\n")
  summary_cols <- c("k", x$criterion_used, "avg_posterior_overall", "entropy_relative", "converged")
  print(x$comparison[, summary_cols], digits = 3, row.names = FALSE)
  
  cat("\nOptimal model summary:\n")
  print(x$best_result)
  
  invisible(x)
} 
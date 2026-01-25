# ==============================================================================
# MARKOV CHAIN ANALYSIS FUNCTIONS
# ==============================================================================
#
# Comprehensive analysis functions for Markov chains including:
# - Extraction utilities for transition matrices and distributions
# - Validation and normalization functions
# - Information-theoretic measures (entropy, Gini)
# - Structural analysis (mixing time, conductance)
# - Full chain analysis with visualization
# - Stability and persistence metrics
#
# ==============================================================================

#' Extract Transition Matrix from Various Input Types
#'
#' This function extracts a transition matrix from different input formats,
#' providing flexibility in how Markov chain data can be specified.
#'
#' @param x Input object containing transition matrix information. Can be:
#'   \itemize{
#'     \item A numeric matrix (returned as-is)
#'     \item A list with element named 'weights', 'transition_matrix', or 'P'
#'     \item A TNA model object
#'   }
#'
#' @return A numeric matrix representing the transition probabilities between states.
#'   Element [i,j] represents the probability of transitioning from state i to state j.
#'
#' @details
#' The function searches for transition matrix information in the following order:
#' \enumerate{
#'   \item If \code{x} is already a matrix, returns it directly
#'   \item If \code{x} is a list, searches for 'weights' field first (TNA convention)
#'   \item Then checks for 'transition_matrix' field
#'   \item Finally checks for 'P' field (common mathematical notation)
#' }
#'
#' @examples
#' # From a matrix
#' P <- matrix(c(0.7, 0.3, 0.4, 0.6), nrow = 2, byrow = TRUE)
#' extract_transition_matrix(P)
#'
#' # From a list with 'weights' field
#' model <- list(weights = P, inits = c(0.5, 0.5))
#' extract_transition_matrix(model)
#'
#' # From a list with 'P' field
#' model2 <- list(P = P)
#' extract_transition_matrix(model2)
#'
#' @export
extract_transition_matrix <- function(x) {
  if (is.matrix(x)) {
    return(x)
  } else if (is.list(x)) {
    if ("weights" %in% names(x)) {
      return(x$weights)
    } else if ("transition_matrix" %in% names(x)) {
      return(x$transition_matrix)
    } else if ("P" %in% names(x)) {
      return(x$P)
    }
  }
  stop("Cannot extract transition matrix from input object. ",
       "Provide a matrix or list with 'weights', 'transition_matrix', or 'P' field.",
       call. = FALSE)
}

#' Extract Initial Distribution from Various Input Types
#'
#' Extracts or generates an initial probability distribution from different input
#' formats, with automatic fallback to uniform distribution if needed.
#'
#' @param x Input object containing initial distribution information. Can be:
#'   \itemize{
#'     \item A numeric vector of probabilities (returned as-is)
#'     \item A list with element named 'inits', 'initial_distribution', or 'pi0'
#'     \item A TNA model object with 'inits' field
#'     \item A TNA model object with 'data' field (computes empirical distribution)
#'   }
#' @param n_states Integer. Number of states in the Markov chain. Used to generate
#'   a uniform distribution if initial distribution cannot be extracted from \code{x}.
#'   Default is NULL.
#'
#' @return A numeric vector of length \code{n_states} representing the initial
#'   probability distribution. Values sum to 1.
#'
#' @details
#' The function searches for initial distribution in the following order:
#' \enumerate{
#'   \item If \code{x} is a numeric vector, returns it directly
#'   \item If \code{x} is a list, searches for 'inits' field (TNA convention)
#'   \item If TNA model has 'data' field, computes empirical initial distribution
#'   \item Then checks for 'initial_distribution' field
#'   \item Then checks for 'pi0' field (mathematical notation)
#'   \item If none found and \code{n_states} is provided, returns uniform distribution
#' }
#'
#' For TNA models with sequence data, the empirical initial distribution is computed
#' from the first time point of all sequences.
#'
#' If the initial distribution cannot be found and \code{n_states} is NULL,
#' the function raises an error.
#'
#' @examples
#' # From a vector
#' init_dist <- c(0.5, 0.3, 0.2)
#' extract_initial_distribution(init_dist)
#'
#' # From a TNA model with inits
#' model <- list(
#'   weights = matrix(c(0.7, 0.3, 0.4, 0.6), 2, 2),
#'   inits = c(0.6, 0.4),
#'   labels = c("Active", "Inactive")
#' )
#' extract_initial_distribution(model)
#'
#' # Generate uniform distribution when not found
#' extract_initial_distribution(list(), n_states = 3)
#'
#' @export
extract_initial_distribution <- function(x, n_states = NULL) {
  if (is.vector(x) && is.numeric(x)) {
    return(x)
  } else if (is.list(x)) {
    # First check for explicit inits field
    if ("inits" %in% names(x)) {
      return(x$inits)
    }
    # If TNA model has data, compute empirical initial distribution
    else if ("data" %in% names(x) && !is.null(x$data)) {
      data_matrix <- x$data
      if (is.matrix(data_matrix) || is.data.frame(data_matrix)) {
        # Get first time point
        first_col <- data_matrix[, 1]
        # Count frequencies
        freq_table <- table(first_col)
        # Get number of states from weights if available
        if ("weights" %in% names(x)) {
          n_states <- nrow(x$weights)
          # Initialize with zeros
          init_dist <- numeric(n_states)
          # Fill in observed frequencies
          for (state_idx in seq_len(n_states)) {
            if (as.character(state_idx) %in% names(freq_table)) {
              init_dist[state_idx] <- freq_table[as.character(state_idx)]
            }
          }
          # Normalize
          init_dist <- init_dist / sum(init_dist)
          # Add names if labels available
          if ("labels" %in% names(x)) {
            names(init_dist) <- x$labels
          }
          return(init_dist)
        }
      }
    }
    else if ("initial_distribution" %in% names(x)) {
      return(x$initial_distribution)
    } else if ("pi0" %in% names(x)) {
      return(x$pi0)
    }
  }

  if (!is.null(n_states)) {
    warning("Could not extract initial distribution. Using uniform distribution.",
            call. = FALSE)
    return(rep(1 / n_states, n_states))
  }

  stop("Cannot extract initial distribution from input object. ",
       "Provide a vector, list with 'inits'/'initial_distribution'/'pi0' field, ",
       "or specify 'n_states' for uniform distribution.",
       call. = FALSE)
}

#' Extract State Names from Various Input Types
#'
#' Extracts state labels from input or generates default names if not available.
#'
#' @param x Input object. Can be:
#'   \itemize{
#'     \item A matrix with column names
#'     \item A list with a 'labels' field
#'     \item A TNA model object
#'   }
#' @param n_states Integer. Number of states in the Markov chain. Used to validate
#'   extracted names and generate default names if needed.
#'
#' @return A character vector of length \code{n_states} containing state labels.
#'   If no labels found, returns default names "S1", "S2", ..., "Sn".
#'
#' @details
#' The function searches for state names in the following order:
#' \enumerate{
#'   \item If \code{x} is a list with 'labels' field of correct length
#'   \item If \code{x} is a matrix with column names of correct length
#'   \item Otherwise generates default names "S1", "S2", etc.
#' }
#'
#' State names are automatically converted to character type.
#'
#' @examples
#' # From matrix column names
#' P <- matrix(c(0.7, 0.3, 0.4, 0.6), nrow = 2)
#' colnames(P) <- c("Active", "Inactive")
#' extract_state_names(P, n_states = 2)
#'
#' # From list with labels
#' model <- list(labels = c("Low", "Medium", "High"))
#' extract_state_names(model, n_states = 3)
#'
#' # Default names when not specified
#' extract_state_names(list(), n_states = 4)
#'
#' @export
extract_state_names <- function(x, n_states) {
  state_names <- NULL

  if (is.list(x) && "labels" %in% names(x)) {
    if (!is.null(x$labels) && length(x$labels) == n_states) {
      state_names <- as.character(x$labels)
    }
  } else if (is.matrix(x)) {
    if (!is.null(colnames(x)) && length(colnames(x)) == n_states) {
      state_names <- colnames(x)
    }
  }

  if (is.null(state_names)) {
    state_names <- paste0("S", 1:n_states)
  }

  return(state_names)
}

#' Check if a Matrix is Stochastic
#'
#' Validates whether a matrix satisfies the requirements of a stochastic
#' (row-stochastic) matrix for Markov chain analysis.
#'
#' @param P Numeric matrix to validate.
#' @param tol Numeric. Tolerance for numerical comparisons. Values within
#'   \code{tol} of 0 or 1 are considered equal. Default is 1e-10.
#'
#' @return Logical. \code{TRUE} if the matrix is stochastic, \code{FALSE} otherwise.
#'
#' @details
#' A row-stochastic matrix must satisfy three conditions:
#' \enumerate{
#'   \item All elements must be non-negative (>= 0)
#'   \item All elements must be <= 1
#'   \item Each row must sum to 1
#' }
#'
#' The \code{tol} parameter allows for small numerical errors that may arise
#' from floating-point arithmetic.
#'
#' @examples
#' # Valid stochastic matrix
#' P1 <- matrix(c(0.7, 0.3, 0.4, 0.6), nrow = 2, byrow = TRUE)
#' is_stochastic(P1)
#'
#' # Invalid: row sums don't equal 1
#' P2 <- matrix(c(0.5, 0.3, 0.4, 0.6), nrow = 2, byrow = TRUE)
#' is_stochastic(P2)
#'
#' # Invalid: negative values
#' P3 <- matrix(c(1.2, -0.2, 0.4, 0.6), nrow = 2, byrow = TRUE)
#' is_stochastic(P3)
#'
#' @export
is_stochastic <- function(P, tol = 1e-10) {
  all(P >= -tol) && all(P <= 1 + tol) && all(abs(rowSums(P) - 1) < tol)
}

#' Normalize a Matrix to Make it Stochastic
#'
#' Converts a non-negative matrix into a row-stochastic matrix by normalizing
#' each row to sum to 1.
#'
#' @param P Numeric matrix with non-negative elements to normalize.
#'
#' @return A row-stochastic matrix where each row sums to 1.
#'
#' @details
#' Each row is divided by its sum to ensure row sums equal 1. If any row
#' sums to zero, a warning is issued and that row is left unchanged (to avoid
#' division by zero).
#'
#' This function is useful when you have a matrix of counts, rates, or
#' unnormalized probabilities that need to be converted into proper transition
#' probabilities.
#'
#' @examples
#' # Normalize a count matrix
#' counts <- matrix(c(10, 5, 20, 30), nrow = 2, byrow = TRUE)
#' P <- normalize_to_stochastic(counts)
#' print(P)
#' rowSums(P)  # Each row sums to 1
#'
#' # Already normalized matrix remains unchanged
#' P_normalized <- matrix(c(0.7, 0.3, 0.4, 0.6), nrow = 2, byrow = TRUE)
#' normalize_to_stochastic(P_normalized)
#'
#' @export
normalize_to_stochastic <- function(P) {
  row_sums <- rowSums(P)
  if (any(row_sums == 0)) {
    warning("Some rows sum to zero; cannot normalize those rows.", call. = FALSE)
    row_sums[row_sums == 0] <- 1
  }
  P / row_sums
}

#' Calculate Shannon Entropy of a Probability Distribution
#'
#' Computes the Shannon entropy, a measure of uncertainty or randomness in a
#' probability distribution.
#'
#' @param prob_dist Numeric vector. A probability distribution (must sum to 1).
#' @param base Numeric. Logarithm base for entropy calculation. Use 2 for
#'   bits (default), \code{exp(1)} for nats, or 10 for hartleys/bans.
#' @param tol Numeric. Tolerance for treating probabilities as zero to avoid
#'   numerical issues with log(0). Default is 1e-10.
#'
#' @return Numeric. The Shannon entropy value. Returns 0 if all mass is on
#'   one outcome (zero uncertainty). Higher values indicate more uniform
#'   distributions (higher uncertainty).
#'
#' @details
#' Shannon entropy is defined as: \deqn{H(X) = -\sum_{i} p_i \log(p_i)}
#'
#' The function:
#' \itemize{
#'   \item Removes probabilities below \code{tol} and NA values
#'   \item Returns 0 if no positive probabilities remain
#'   \item Uses natural logarithm base internally, then converts
#' }
#'
#' Maximum entropy occurs when all outcomes are equally likely.
#'
#' @examples
#' # Uniform distribution (maximum entropy)
#' calculate_entropy(c(0.25, 0.25, 0.25, 0.25))  # 2 bits
#'
#' # Concentrated distribution (low entropy)
#' calculate_entropy(c(0.9, 0.05, 0.03, 0.02))
#'
#' # Binary distribution
#' calculate_entropy(c(0.5, 0.5))  # 1 bit
#'
#' # Using natural logarithm (nats)
#' calculate_entropy(c(0.5, 0.5), base = exp(1))
#'
#' @references
#' Shannon, C. E. (1948). A mathematical theory of communication.
#' Bell System Technical Journal, 27(3), 379-423.
#'
#' @seealso \code{\link{calculate_gini_impurity}} for alternative uncertainty measure
#'
#' @export
calculate_entropy <- function(prob_dist, base = 2, tol = 1e-10) {
  prob_dist_nz <- prob_dist[prob_dist > tol & !is.na(prob_dist)]
  if (length(prob_dist_nz) == 0) return(0.0)
  -sum(prob_dist_nz * log(prob_dist_nz, base = base))
}

#' Calculate Gini Impurity of a Probability Distribution
#'
#' Computes the Gini impurity, a measure of how often a randomly chosen
#' element would be incorrectly classified if randomly labeled according
#' to the distribution.
#'
#' @param prob_dist Numeric vector. A probability distribution (should sum to 1).
#'
#' @return Numeric. The Gini impurity value between 0 and 1. Returns \code{NA}
#'   if the distribution is empty. A value of 0 indicates perfect purity
#'   (all mass on one outcome), while higher values indicate more impurity.
#'
#' @details
#' Gini impurity is defined as: \deqn{G = 1 - \sum_{i} p_i^2}
#'
#' Maximum impurity occurs with a uniform distribution. For a distribution
#' over \eqn{k} outcomes, maximum impurity is \eqn{(k-1)/k}.
#'
#' NA values in the distribution are automatically removed before calculation.
#'
#' Gini impurity is commonly used in decision trees and is related to, but
#' different from, entropy.
#'
#' @examples
#' # Pure distribution (minimum impurity)
#' calculate_gini_impurity(c(1, 0, 0))  # 0
#'
#' # Uniform distribution (maximum impurity for 3 states)
#' calculate_gini_impurity(c(1/3, 1/3, 1/3))  # 0.667
#'
#' # Binary case
#' calculate_gini_impurity(c(0.5, 0.5))  # 0.5
#' calculate_gini_impurity(c(0.9, 0.1))  # 0.18
#'
#' @references
#' Gini, C. (1912). Variabilita e mutabilita.
#' Reprinted in Memorie di metodologica statistica (Ed. Pizetti E, Salvemini, T).
#' Rome: Libreria Eredi Virgilio Veschi (1955).
#'
#' @seealso \code{\link{calculate_entropy}} for alternative uncertainty measure
#'
#' @export
calculate_gini_impurity <- function(prob_dist) {
  prob_dist_clean <- prob_dist[!is.na(prob_dist)]
  if (length(prob_dist_clean) == 0) return(NA_real_)
  1.0 - sum(prob_dist_clean^2)
}

#' Calculate Mixing Time for a Markov Chain
#'
#' Computes the number of steps required for the chain to converge to its
#' stationary distribution from any starting state.
#'
#' @param P Numeric matrix. A stochastic transition matrix.
#' @param epsilon Numeric. Convergence threshold for total variation distance.
#'   The chain is considered mixed when the distance between the current
#'   distribution and equilibrium is less than \code{epsilon}. Default is 1e-10.
#' @param max_iter Integer. Maximum number of iterations to prevent infinite loops.
#'   Default is 1000.
#'
#' @return Integer. The mixing time (number of steps to convergence). If
#'   convergence is not reached within \code{max_iter} steps, returns
#'   \code{max_iter} with a warning.
#'
#' @details
#' Mixing time is the number of steps \eqn{t} such that:
#' \deqn{\max_{i,j} |P^t(i,j) - \pi(j)| < \epsilon}
#'
#' where \eqn{P^t} is the t-step transition matrix and \eqn{\pi} is the
#' stationary distribution.
#'
#' The function:
#' \enumerate{
#'   \item Computes the equilibrium distribution from the principal eigenvector
#'   \item Iteratively multiplies the transition matrix
#'   \item Checks convergence at each step
#' }
#'
#' Smaller mixing times indicate faster convergence.
#'
#' @examples
#' # Fast mixing chain
#' P_fast <- matrix(c(0.1, 0.9, 0.9, 0.1), nrow = 2, byrow = TRUE)
#' calculate_mixing_time(P_fast)
#'
#' # Slower mixing chain
#' P_slow <- matrix(c(0.9, 0.1, 0.1, 0.9), nrow = 2, byrow = TRUE)
#' calculate_mixing_time(P_slow)
#'
#' @export
calculate_mixing_time <- function(P, epsilon = 1e-10, max_iter = 1000) {
  n_states <- nrow(P)
  P_power <- P
  mixing_time <- 1

  # Get equilibrium distribution
  eigen_result <- eigen(t(P))
  principal_idx <- which.min(abs(eigen_result$values - 1))
  eq_dist <- Re(eigen_result$vectors[, principal_idx])
  eq_dist <- eq_dist / sum(eq_dist)

  eq_matrix <- matrix(eq_dist, n_states, n_states, byrow = TRUE)

  while(mixing_time < max_iter) {
    max_diff <- max(abs(P_power - eq_matrix))
    if (max_diff < epsilon) break
    P_power <- P_power %*% P
    mixing_time <- mixing_time + 1
  }

  if (mixing_time >= max_iter) {
    warning("Mixing time calculation reached maximum iterations")
  }

  return(mixing_time)
}

#' Calculate Conductance of a Markov Chain
#'
#' @param P Transition matrix
#' @param max_subset_size Maximum subset size to check (for efficiency)
#' @return Numeric conductance value
#' @export
calculate_conductance <- function(P, max_subset_size = NULL) {
  n_states <- nrow(P)

  if (is.null(max_subset_size)) {
    max_subset_size <- min(floor(n_states / 2), 10)
  }

  eigen_result <- eigen(t(P))
  principal_idx <- which.min(abs(eigen_result$values - 1))
  pi <- Re(eigen_result$vectors[, principal_idx])
  pi <- pi / sum(pi)

  min_conductance <- Inf

  for(i in 1:min(max_subset_size, n_states - 1)) {
    subsets <- utils::combn(1:n_states, i, simplify = FALSE)

    for(subset in subsets) {
      complement <- setdiff(1:n_states, subset)
      flow <- sum(pi[subset] * P[subset, complement])
      vol <- sum(pi[subset])

      if(vol > 0 && vol <= 0.5) {
        conductance <- flow / vol
        min_conductance <- min(min_conductance, conductance)
      }
    }
  }

  return(min_conductance)
}

#' Comprehensive Markov Chain Analysis
#'
#' Performs a complete analysis of a discrete-time Markov chain, computing
#' structural properties, equilibrium distributions, passage times, convergence
#' metrics, and time evolution dynamics.
#'
#' @param x Input object containing the Markov chain. Can be:
#'   \itemize{
#'     \item A numeric stochastic matrix (rows sum to 1)
#'     \item A list with fields 'weights', 'transition_matrix', or 'P' for the matrix
#'     \item A TNA model object
#'   }
#' @param initial_distribution Numeric vector. Initial probability distribution
#'   over states. If NULL (default), extracts from \code{x} or uses uniform
#'   distribution. Must sum to 1.
#' @param num_steps Integer. Number of time steps for evolution analysis.
#'   Default is 100. Larger values provide more detailed convergence tracking
#'   but increase computation time.
#' @param verbose Logical. If TRUE (default), prints a formatted summary of
#'   results to console.
#' @param calculate_all Logical. If TRUE (default), computes all properties
#'   including mixing time, conductance, and hitting probabilities. Set to
#'   FALSE for large chains (>20 states) to improve performance.
#' @param tol Numeric. Tolerance for numerical comparisons and convergence
#'   tests. Default is 1e-10.
#'
#' @return An object of class 'markov_analysis' (a list) containing:
#'   \describe{
#'     \item{n_states}{Number of states}
#'     \item{state_names}{Character vector of state labels}
#'     \item{transition_matrix}{The transition probability matrix}
#'     \item{initial_distribution}{The initial distribution used}
#'     \item{eigenvalues}{Complex vector of eigenvalues}
#'     \item{equilibrium_distribution}{Stationary/equilibrium distribution}
#'     \item{stationary_found}{Logical indicating if stationary distribution exists}
#'     \item{mean_first_passage_times}{Matrix of expected hitting times}
#'     \item{fundamental_matrix}{Fundamental matrix Z}
#'     \item{hitting_probabilities}{List of hitting probability vectors}
#'     \item{absorption_times}{Expected absorption times by state}
#'     \item{sojourn_times}{Expected time spent in each state per visit}
#'     \item{accessibility_matrix}{Logical matrix of state reachability}
#'     \item{communicating_classes}{List of communicating class memberships}
#'     \item{period}{Period of the chain (1 = aperiodic)}
#'     \item{is_irreducible}{Logical, TRUE if chain is irreducible}
#'     \item{is_regular}{Logical, TRUE if chain is regular}
#'     \item{is_reversible}{Logical, TRUE if chain satisfies detailed balance}
#'     \item{is_stochastic}{Logical, TRUE if input matrix is row-stochastic}
#'     \item{spectral_gap}{1 minus second-largest eigenvalue modulus}
#'     \item{slem}{Second-largest eigenvalue modulus}
#'     \item{mixing_time}{Number of steps to convergence}
#'     \item{conductance}{Minimum conductance (bottleneck measure)}
#'     \item{entropy_rate}{Long-run entropy rate}
#'     \item{num_steps}{Number of evolution steps computed}
#'     \item{distributions_over_time}{List of distributions at each time}
#'     \item{distance_to_stationary}{L1 distance to equilibrium over time}
#'     \item{change_in_distribution}{Change in distribution between steps}
#'     \item{entropy_over_time}{Entropy at each time step}
#'     \item{gini_over_time}{Gini impurity at each time step}
#'     \item{expected_transitions_over_time}{Expected number of transitions}
#'     \item{absorbing_chain_analysis}{List with absorbing state analysis}
#'   }
#'
#' @details
#' This function provides a one-stop analysis for Markov chains, computing:
#'
#' \strong{Structural Properties:}
#' \itemize{
#'   \item Irreducibility: Can all states reach each other?
#'   \item Periodicity: Does the chain have cyclic behavior?
#'   \item Communicating classes: Groups of mutually accessible states
#'   \item Reversibility: Does detailed balance hold?
#' }
#'
#' \strong{Equilibrium Analysis:}
#' \itemize{
#'   \item Stationary distribution (from eigenvector of P')
#'   \item Entropy rate of the stationary process
#' }
#'
#' \strong{Passage Time Analysis:}
#' \itemize{
#'   \item Mean first passage times between all state pairs
#'   \item Mean recurrence times (diagonal elements)
#'   \item Expected sojourn times
#' }
#'
#' \strong{Convergence Analysis:}
#' \itemize{
#'   \item Spectral gap (convergence rate indicator)
#'   \item Mixing time (steps to epsilon-convergence)
#'   \item Conductance (bottleneck measure)
#' }
#'
#' \strong{Time Evolution:}
#' \itemize{
#'   \item Distribution evolution from initial state
#'   \item Convergence to equilibrium
#'   \item Entropy and Gini impurity trajectories
#' }
#'
#' \strong{Absorbing Chain Analysis:}
#' If absorbing states exist, computes absorption probabilities and
#' expected absorption times from transient states.
#'
#' @note
#' For large chains (>20 states), set \code{calculate_all = FALSE} to skip
#' computationally intensive calculations (mixing time, conductance, hitting
#' probabilities).
#'
#' The function handles non-stochastic matrices gracefully but issues warnings.
#' Small negative values (< tol) are automatically corrected to zero.
#'
#' @examples
#' # Basic 2-state chain
#' P <- matrix(c(0.7, 0.3, 0.4, 0.6), nrow = 2, byrow = TRUE)
#' colnames(P) <- c("Healthy", "Sick")
#' results <- analyze_markov_chain(P, num_steps = 50)
#'
#' # Access specific results
#' results$equilibrium_distribution
#' results$mean_first_passage_times
#'
#' @references
#' Kemeny, J. G., & Snell, J. L. (1976). Finite Markov Chains.
#' Springer-Verlag New York.
#'
#' Norris, J. R. (1998). Markov Chains. Cambridge University Press.
#'
#' @seealso
#' \code{\link{plot.markov_analysis}} for visualization,
#' \code{\link{analyze_tna_model}} for TNA-specific wrapper
#'
#' @export
analyze_markov_chain <- function(x,
                                 initial_distribution = NULL,
                                 num_steps = 100,
                                 verbose = TRUE,
                                 calculate_all = TRUE,
                                 tol = 1e-10) {

  # Extract transition matrix
  P <- extract_transition_matrix(x)

  # Basic input validation
  if (!is.matrix(P) || !is.numeric(P))
    stop("Input must be a numeric matrix")
  if (nrow(P) != ncol(P))
    stop("Input must be a square matrix")

  n_states <- nrow(P)

  # Clean small negative values
  P[P < 0 & P > -tol] <- 0

  if (!is_stochastic(P, tol)) {
    warning("Input matrix is not stochastic. Consider normalizing.")
  }

  # Extract or set initial distribution
  if (is.null(initial_distribution)) {
    initial_distribution <- extract_initial_distribution(x, n_states)
  }

  # Extract state names
  state_names <- extract_state_names(x, n_states)
  colnames(P) <- rownames(P) <- state_names
  names(initial_distribution) <- state_names

  # Compute eigenvalues and eigenvectors
  eigen_result <- eigen(t(P))
  principal_idx <- which.min(abs(eigen_result$values - 1))

  # Calculate equilibrium distribution
  eq_dist <- Re(eigen_result$vectors[, principal_idx])
  eq_dist[eq_dist < 0 & eq_dist > -tol] <- 0
  eq_dist <- eq_dist / sum(eq_dist)
  names(eq_dist) <- state_names

  stationary_found <- abs(eigen_result$values[principal_idx] - 1) < tol &&
    all(eq_dist >= -tol)

  # Calculate mean first passage times
  mfpts <- matrix(NA, nrow = n_states, ncol = n_states)
  colnames(mfpts) <- rownames(mfpts) <- state_names

  for (target in 1:n_states) {
    if (n_states > 1) {
      Q <- P[-target, -target, drop = FALSE]
      I <- diag(n_states - 1)
      tryCatch({
        fund_matrix <- solve(I - Q)
        mfpts[-target, target] <- rowSums(fund_matrix)
      }, error = function(e) {
        warning("Could not calculate MFPT for target state ", target)
      })
    }
  }

  # Fill diagonal with mean recurrence times
  if (stationary_found) {
    for (i in 1:n_states) {
      if (eq_dist[i] > tol) {
        mfpts[i, i] <- 1 / eq_dist[i]
      }
    }
  }

  # Calculate fundamental matrix
  I <- diag(n_states)
  Z <- tryCatch({
    solve(I - P + matrix(eq_dist, n_states, n_states, byrow = TRUE))
  }, error = function(e) {
    warning("Could not calculate fundamental matrix")
    NULL
  })

  # Determine chain properties
  is_irreducible <- all(is.finite(mfpts) | row(mfpts) == col(mfpts))
  eigenvalues <- eigen_result$values
  spectral_gap <- 1 - max(abs(eigenvalues[-principal_idx]))

  # Calculate period
  period <- 1
  if (n_states > 1) {
    P_power <- P
    for (i in 2:min(n_states, 50)) {
      P_power <- P_power %*% P
      if (all(abs(diag(P_power)) > tol)) {
        period <- i
        break
      }
    }
  }

  is_regular <- period == 1 && is_irreducible

  # Calculate SLEM (Second Largest Eigenvalue Modulus)
  moduli <- Mod(eigenvalues)
  moduli_not_one <- moduli[abs(moduli - 1.0) > tol]
  slem <- if (length(moduli_not_one) > 0) max(moduli_not_one) else NA

  # Time evolution analysis
  distributions_over_time <- list()
  entropy_over_time <- numeric(num_steps + 1)
  gini_over_time <- numeric(num_steps + 1)
  expected_transitions_over_time <- numeric(num_steps + 1)
  distance_to_stationary <- numeric(num_steps + 1)
  change_in_distribution <- numeric(num_steps)

  current_dist <- initial_distribution
  prev_dist <- current_dist

  for (k in 0:num_steps) {
    distributions_over_time[[k + 1]] <- current_dist
    entropy_over_time[k + 1] <- calculate_entropy(current_dist, base = 2, tol = tol)
    gini_over_time[k + 1] <- calculate_gini_impurity(current_dist)

    # Expected transitions
    prob_leaving <- 1 - diag(P)
    expected_transitions_over_time[k + 1] <- sum(current_dist * prob_leaving)

    # Distance to stationary
    if (stationary_found) {
      distance_to_stationary[k + 1] <- sum(abs(current_dist - eq_dist))
    }

    # Change from previous step
    if (k > 0) {
      change_in_distribution[k] <- sum(abs(current_dist - prev_dist))
    }

    # Update for next iteration
    if (k < num_steps) {
      prev_dist <- current_dist
      current_dist <- as.vector(current_dist %*% P)
      names(current_dist) <- state_names
    }
  }

  # Sojourn times
  diag_vals <- diag(P)
  sojourn_times <- ifelse(abs(diag_vals - 1.0) < tol, Inf, 1 / (1 - diag_vals))
  names(sojourn_times) <- state_names

  # Absorbing state analysis
  absorbing_indices <- which(abs(diag_vals - 1.0) < tol)
  is_absorbing_chain <- length(absorbing_indices) > 0

  absorbing_analysis <- list(
    is_absorbing_chain = is_absorbing_chain,
    absorbing_states_indices = absorbing_indices,
    absorbing_states_names = if(is_absorbing_chain) state_names[absorbing_indices] else NULL,
    transient_states_indices = NULL,
    transient_states_names = NULL,
    Q_matrix = NULL,
    R_matrix = NULL,
    fundamental_matrix_N = NULL,
    expected_absorption_times = NULL,
    absorption_probabilities_B = NULL
  )

  if (is_absorbing_chain) {
    transient_idx <- setdiff(1:n_states, absorbing_indices)
    absorbing_analysis$transient_states_indices <- transient_idx
    absorbing_analysis$transient_states_names <- state_names[transient_idx]

    if (length(transient_idx) > 0) {
      Q <- P[transient_idx, transient_idx, drop = FALSE]
      R_mat <- P[transient_idx, absorbing_indices, drop = FALSE]

      absorbing_analysis$Q_matrix <- Q
      absorbing_analysis$R_matrix <- R_mat

      tryCatch({
        I_Q <- diag(length(transient_idx))
        N_fund <- solve(I_Q - Q)
        rownames(N_fund) <- colnames(N_fund) <- state_names[transient_idx]
        absorbing_analysis$fundamental_matrix_N <- N_fund

        # Expected absorption times
        exp_times <- N_fund %*% rep(1, ncol(N_fund))
        rownames(exp_times) <- state_names[transient_idx]
        colnames(exp_times) <- "Expected_Steps"
        absorbing_analysis$expected_absorption_times <- exp_times

        # Absorption probabilities
        if (ncol(R_mat) > 0) {
          B_probs <- N_fund %*% R_mat
          rownames(B_probs) <- state_names[transient_idx]
          colnames(B_probs) <- state_names[absorbing_indices]
          absorbing_analysis$absorption_probabilities_B <- B_probs
        }
      }, error = function(e) {
        warning("Could not compute absorbing chain metrics: ", e$message)
      })
    }
  }

  # Calculate optional expensive properties
  mixing_time <- NA
  conductance <- NA
  hitting_probs <- NULL

  if (calculate_all && n_states <= 20) {
    mixing_time <- tryCatch(
      calculate_mixing_time(P, epsilon = tol),
      error = function(e) NA
    )

    conductance <- tryCatch(
      calculate_conductance(P),
      error = function(e) NA
    )

    hitting_probs <- lapply(1:n_states, function(i) {
      tryCatch({
        A <- diag(n_states) - P
        A[i, ] <- 0
        A[i, i] <- 1
        b <- rep(0, n_states)
        b[i] <- 1
        solve(A, b)
      }, error = function(e) NULL)
    })
    names(hitting_probs) <- state_names
  }

  # Entropy rate
  entropy_rate <- NA
  if (stationary_found) {
    entropy_rate <- tryCatch({
      H_rate <- 0.0
      for (i in 1:n_states) {
        if (eq_dist[i] > tol) {
          for (j in 1:n_states) {
            if (P[i, j] > tol) {
              H_rate <- H_rate - (eq_dist[i] * P[i, j] * log2(P[i, j]))
            }
          }
        }
      }
      H_rate
    }, error = function(e) NA)
  }

  # Reversibility check
  reversible <- tryCatch({
    if (!stationary_found) return(NA)
    max_violation <- 0
    for (i in 1:n_states) {
      for (j in 1:n_states) {
        violation <- abs(eq_dist[i] * P[i, j] - eq_dist[j] * P[j, i])
        max_violation <- max(max_violation, violation)
      }
    }
    max_violation < tol
  }, error = function(e) NA)

  absorption_times <- if (!is.null(Z)) {
    tryCatch(colSums(Z), error = function(e) rep(NA, n_states))
  } else {
    rep(NA, n_states)
  }
  names(absorption_times) <- state_names

  # Calculate state accessibility matrix
  accessibility <- matrix(FALSE, n_states, n_states)
  colnames(accessibility) <- rownames(accessibility) <- state_names
  for (i in 1:n_states) {
    for (j in 1:n_states) {
      accessibility[i, j] <- is.finite(mfpts[i, j]) || i == j
    }
  }

  # Identify communicating classes
  classes <- list()
  unvisited <- 1:n_states
  while (length(unvisited) > 0) {
    current <- unvisited[1]
    class_members <- which(accessibility[current, ] & accessibility[, current])
    classes <- c(classes, list(state_names[class_members]))
    unvisited <- setdiff(unvisited, class_members)
  }

  # Create result object
  result <- list(
    # Basic info
    n_states = n_states,
    state_names = state_names,
    transition_matrix = P,
    initial_distribution = initial_distribution,

    # Eigenvalues and equilibrium
    eigenvalues = eigenvalues,
    equilibrium_distribution = eq_dist,
    stationary_found = stationary_found,

    # Passage times
    mean_first_passage_times = mfpts,
    fundamental_matrix = Z,
    hitting_probabilities = hitting_probs,
    absorption_times = absorption_times,
    sojourn_times = sojourn_times,

    # Structure
    accessibility_matrix = accessibility,
    communicating_classes = classes,
    period = period,
    is_irreducible = is_irreducible,
    is_regular = is_regular,
    is_reversible = reversible,
    is_stochastic = is_stochastic(P, tol),

    # Convergence
    spectral_gap = spectral_gap,
    slem = slem,
    mixing_time = mixing_time,
    conductance = conductance,

    # Information theory
    entropy_rate = entropy_rate,

    # Time evolution
    num_steps = num_steps,
    distributions_over_time = distributions_over_time,
    distance_to_stationary = distance_to_stationary,
    change_in_distribution = change_in_distribution,
    entropy_over_time = entropy_over_time,
    gini_over_time = gini_over_time,
    expected_transitions_over_time = expected_transitions_over_time,

    # Absorbing chain analysis
    absorbing_chain_analysis = absorbing_analysis
  )

  class(result) <- c("markov_analysis", "list")

  if (verbose) {
    print(result)
  }

  return(invisible(result))
}

#' Print method for markov_analysis objects
#'
#' @param x A markov_analysis object
#' @param ... Additional arguments (unused)
#' @export
print.markov_analysis <- function(x, ...) {
  cat("\n===================================================================\n")
  cat("           MARKOV CHAIN ANALYSIS RESULTS\n")
  cat("===================================================================\n")

  cat("\n-- BASIC INFORMATION\n")
  cat("   Number of states:", x$n_states, "\n")
  cat("   State names:", paste(x$state_names, collapse = ", "), "\n")
  cat("   Is stochastic:", x$is_stochastic, "\n")

  cat("\n-- EQUILIBRIUM DISTRIBUTION\n")
  if (x$stationary_found) {
    eq_df <- data.frame(
      State = names(x$equilibrium_distribution),
      Probability = round(x$equilibrium_distribution, 4)
    )
    print(eq_df, row.names = FALSE)
    cat("   Entropy:", round(calculate_entropy(x$equilibrium_distribution), 4), "bits\n")
  } else {
    cat("   Stationary distribution not found\n")
  }

  cat("\n-- STRUCTURAL PROPERTIES\n")
  cat("   Period:", x$period, "\n")
  cat("   Irreducible:", x$is_irreducible, "\n")
  cat("   Regular:", x$is_regular, "\n")
  cat("   Reversible:", if(is.na(x$is_reversible)) "Unknown" else x$is_reversible, "\n")
  cat("   Number of communicating classes:", length(x$communicating_classes), "\n")

  if (length(x$communicating_classes) > 1) {
    cat("   Communicating classes:\n")
    for (i in seq_along(x$communicating_classes)) {
      cat("     Class", i, ":", paste(x$communicating_classes[[i]], collapse = ", "), "\n")
    }
  }

  if (x$absorbing_chain_analysis$is_absorbing_chain) {
    cat("\n-- ABSORBING STATES\n")
    cat("   Absorbing states:", paste(x$absorbing_chain_analysis$absorbing_states_names, collapse = ", "), "\n")
    cat("   Transient states:", paste(x$absorbing_chain_analysis$transient_states_names, collapse = ", "), "\n")

    if (!is.null(x$absorbing_chain_analysis$expected_absorption_times)) {
      cat("   Expected absorption times:\n")
      print(round(x$absorbing_chain_analysis$expected_absorption_times, 2))
    }
  }

  cat("\n-- CONVERGENCE PROPERTIES\n")
  cat("   Spectral gap:", round(x$spectral_gap, 4), "\n")
  if (!is.na(x$slem)) {
    cat("   SLEM (2nd largest eigenvalue modulus):", round(x$slem, 4), "\n")
  }
  if (!is.na(x$mixing_time)) {
    cat("   Mixing time:", x$mixing_time, "\n")
  }
  if (!is.na(x$conductance)) {
    cat("   Conductance:", round(x$conductance, 4), "\n")
  }
  if (!is.na(x$entropy_rate)) {
    cat("   Entropy rate:", round(x$entropy_rate, 4), "bits\n")
  }

  cat("\n-- TIME EVOLUTION (", x$num_steps, " steps)\n", sep="")
  if (x$stationary_found) {
    final_dist <- sum(abs(x$distributions_over_time[[length(x$distributions_over_time)]] -
                            x$equilibrium_distribution))
    cat("   Final distance to stationary:", round(final_dist, 6), "\n")
  }
  cat("   Final entropy:", round(x$entropy_over_time[length(x$entropy_over_time)], 4), "bits\n")
  cat("   Final Gini impurity:", round(x$gini_over_time[length(x$gini_over_time)], 4), "\n")

  cat("\n-- MEAN FIRST PASSAGE TIMES\n")
  print(round(x$mean_first_passage_times, 2))

  cat("\n===================================================================\n")
  invisible(x)
}

#' Plot Markov chain analysis results
#'
#' @param x A markov_analysis object
#' @param type Type of plot ("distance", "distributions", "entropy", "eigenvalues", "mfpt")
#' @param ... Additional plotting parameters
#' @export
plot.markov_analysis <- function(x, type = "distance", ...) {

  if (type == "distance") {
    if (x$stationary_found && length(x$distance_to_stationary) > 0) {
      graphics::plot(0:(length(x$distance_to_stationary) - 1), x$distance_to_stationary,
           type = 'l', col = 'blue', lwd = 2,
           xlab = "Time Step", ylab = "L1 Distance to Stationary",
           main = "Convergence to Stationary Distribution",
           panel.first = graphics::grid(), ...)
    } else {
      graphics::plot(1, type="n", axes=FALSE, xlab="", ylab="",
           main="Distance to Stationary: Not Available")
      graphics::text(1, 1, "Stationary distribution not found", cex = 0.8)
    }

  } else if (type == "distributions") {
    if (length(x$distributions_over_time) > 0) {
      n_steps <- length(x$distributions_over_time)
      prob_matrix <- do.call(rbind, x$distributions_over_time)

      colors <- grDevices::rainbow(x$n_states)

      graphics::plot(0:(n_steps - 1), prob_matrix[, 1],
           type = 'l', col = colors[1], lwd = 2,
           xlab = "Time Step", ylab = "Probability",
           main = "Probability Distribution Over Time",
           ylim = c(0, 1),
           panel.first = graphics::grid(), ...)

      if (x$n_states > 1) {
        for (i in 2:x$n_states) {
          graphics::lines(0:(n_steps - 1), prob_matrix[, i], col = colors[i], lwd = 2)
        }
      }

      graphics::legend("topright", legend = x$state_names, col = colors, lwd = 2, cex = 0.7)
    }

  } else if (type == "entropy") {
    graphics::plot(0:(length(x$entropy_over_time) - 1), x$entropy_over_time,
         type = 'l', col = 'darkgreen', lwd = 2,
         xlab = "Time Step", ylab = "Entropy (bits)",
         main = "Entropy Over Time",
         panel.first = graphics::grid(), ...)

  } else if (type == "eigenvalues") {
    graphics::plot(Re(x$eigenvalues), Im(x$eigenvalues),
         pch = 19, col = 'steelblue', cex = 1.5,
         xlab = "Real Part", ylab = "Imaginary Part",
         main = "Eigenvalue Spectrum",
         panel.first = graphics::grid(), ...)
    graphics::abline(v = 0, h = 0, lty = 2, col = "gray")

    # Add unit circle
    theta <- seq(0, 2*pi, length.out = 100)
    graphics::lines(cos(theta), sin(theta), lty = 2, col = "red")

  } else if (type == "mfpt") {
    if (!is.null(x$mean_first_passage_times)) {
      graphics::heatmap(x$mean_first_passage_times,
              main = "Mean First Passage Times",
              col = grDevices::heat.colors(20),
              Rowv = NA, Colv = NA, scale = "none", ...)
    }
  }
}

#' Compare multiple Markov chains
#'
#' @param ... Named Markov chain matrices or analysis objects
#' @return Data frame comparing chain properties
#' @export
compare_markov_chains <- function(...) {
  chains <- list(...)
  chain_names <- names(chains)
  if (is.null(chain_names)) {
    chain_names <- paste0("Chain_", seq_along(chains))
  }

  results <- lapply(seq_along(chains), function(i) {
    analysis <- if (inherits(chains[[i]], "markov_analysis")) {
      chains[[i]]
    } else {
      analyze_markov_chain(chains[[i]], verbose = FALSE, calculate_all = FALSE)
    }

    data.frame(
      Chain = chain_names[i],
      States = analysis$n_states,
      Irreducible = analysis$is_irreducible,
      Period = analysis$period,
      Regular = analysis$is_regular,
      Reversible = if(is.na(analysis$is_reversible)) "Unknown" else analysis$is_reversible,
      SpectralGap = round(analysis$spectral_gap, 4),
      SLEM = if(is.na(analysis$slem)) NA else round(analysis$slem, 4),
      EntropyRate = if(is.na(analysis$entropy_rate)) NA else round(analysis$entropy_rate, 4),
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, results)
}

#' Simulate a Markov chain trajectory
#'
#' @param P Transition matrix or analysis object
#' @param n_steps Number of steps to simulate
#' @param initial_state Initial state (index or name)
#' @param return_path Return full path or just final state
#' @return Vector of visited states
#' @export
simulate_markov_chain <- function(P, n_steps, initial_state = 1, return_path = TRUE) {
  P <- extract_transition_matrix(P)
  n_states <- nrow(P)

  if (is.character(initial_state)) {
    initial_state <- which(colnames(P) == initial_state)
    if (length(initial_state) == 0) stop("State name not found")
  }

  states <- integer(n_steps + 1)
  states[1] <- initial_state

  for (i in 2:(n_steps + 1)) {
    states[i] <- sample(1:n_states, 1, prob = P[states[i-1], ])
  }

  if (!is.null(colnames(P))) {
    states <- colnames(P)[states]
  }

  if (return_path) {
    return(states)
  } else {
    return(states[length(states)])
  }
}

#' Calculate return times for a state
#'
#' @param P Transition matrix
#' @param state State index or name
#' @param n_simulations Number of simulations
#' @param max_steps Maximum steps per simulation
#' @return Vector of return times
#' @export
calculate_return_times <- function(P, state = 1, n_simulations = 1000, max_steps = 1000) {
  P <- extract_transition_matrix(P)

  if (is.character(state)) {
    state <- which(colnames(P) == state)
    if (length(state) == 0) stop("State name not found")
  }

  return_times <- numeric(n_simulations)

  for (sim in 1:n_simulations) {
    path <- simulate_markov_chain(P, max_steps, initial_state = state, return_path = TRUE)
    if (is.character(path)) {
      path <- match(path, colnames(P))
    }
    first_return <- which(path[-1] == state)[1]
    return_times[sim] <- if (is.na(first_return)) NA else first_return
  }

  return_times[!is.na(return_times)]
}

#' Generate multiple plots for Markov chain analysis
#'
#' @param x A markov_analysis object
#' @param ask Ask before each plot (interactive sessions)
#' @export
plot_all <- function(x, ask = interactive()) {
  if (!inherits(x, "markov_analysis")) {
    stop("Input must be a markov_analysis object")
  }

  old_ask <- graphics::par("ask")
  if (ask) graphics::par(ask = TRUE)

  tryCatch({
    # Plot 1: Distance to stationary
    if (x$stationary_found) {
      plot(x, type = "distance")
    }

    # Plot 2: Distributions over time
    plot(x, type = "distributions")

    # Plot 3: Entropy over time
    plot(x, type = "entropy")

    # Plot 4: Eigenvalue spectrum
    plot(x, type = "eigenvalues")

    # Plot 5: Mean first passage times heatmap
    if (!is.null(x$mean_first_passage_times) && any(is.finite(x$mean_first_passage_times))) {
      plot(x, type = "mfpt")
    }
  }, finally = {
    graphics::par(ask = old_ask)
  })
}

#' Analyze TNA (Temporal Network Analysis) Model
#'
#' Convenience wrapper around \code{\link{analyze_markov_chain}} specifically
#' designed for TNA model objects. Provides automatic extraction of transition
#' matrices, initial distributions, and state labels from TNA format.
#'
#' @param tna_model List. A TNA model object containing:
#'   \describe{
#'     \item{weights}{Required. Numeric matrix of transition probabilities}
#'     \item{inits}{Optional. Numeric vector of initial state probabilities}
#'     \item{labels}{Optional. Character vector of state names}
#'   }
#' @param num_steps Integer. Number of time steps for evolution analysis.
#'   Default is 100.
#' @param tol Numeric. Tolerance for numerical comparisons. Default is 1e-9.
#' @param plot_results Logical. If TRUE (default), generates comprehensive
#'   plots including convergence, distributions, entropy, eigenvalues, and
#'   mean first passage times. Plots are displayed interactively if in an
#'   interactive session.
#' @param verbose Logical. If TRUE (default), prints detailed formatted
#'   analysis results to console.
#' @param calculate_all Logical. If TRUE (default), computes all metrics
#'   including computationally intensive ones. Set to FALSE for large models.
#'
#' @return An object of class 'markov_analysis' containing all analysis
#'   results. See \code{\link{analyze_markov_chain}} for detailed structure.
#'
#' @details
#' This function:
#' \enumerate{
#'   \item Validates the TNA model structure
#'   \item Extracts transition matrix, initial distribution, and labels
#'   \item Calls \code{\link{analyze_markov_chain}} for core analysis
#'   \item Optionally generates comprehensive visualizations
#' }
#'
#' The function provides informative error messages if the model structure
#' is invalid or dimensions are inconsistent.
#'
#' @section TNA Model Format:
#' TNA models are lists with specific fields:
#' \itemize{
#'   \item \strong{weights}: The core transition matrix (required)
#'   \item \strong{inits}: Initial probability distribution (optional, defaults to uniform)
#'   \item \strong{labels}: Human-readable state names (optional, defaults to "S1", "S2", ...)
#' }
#'
#' @examples
#' # Create a simple TNA model
#' model <- list(
#'   weights = matrix(c(0.7, 0.3, 0.2, 0.8), nrow = 2, byrow = TRUE),
#'   inits = c(0.5, 0.5),
#'   labels = c("Active", "Inactive")
#' )
#'
#' # Analyze with plots
#' results <- analyze_tna_model(model, num_steps = 50)
#'
#' @seealso
#' \code{\link{analyze_markov_chain}} for general Markov chain analysis,
#' \code{\link{create_example_tna_model}} for example models
#'
#' @export
analyze_tna_model <- function(tna_model,
                              num_steps = 100,
                              tol = 1e-9,
                              plot_results = TRUE,
                              verbose = TRUE,
                              calculate_all = TRUE) {

  # Validate TNA model structure
  if (!is.list(tna_model)) {
    stop("tna_model must be a list")
  }

  if (!("weights" %in% names(tna_model))) {
    stop("tna_model must contain 'weights' (transition matrix)")
  }

  # Extract components
  P <- tna_model$weights

  # Get initial distribution
  initial_dist <- NULL
  if ("inits" %in% names(tna_model)) {
    initial_dist <- tna_model$inits
  }

  # Validate dimensions
  if (!is.matrix(P)) {
    stop("model$weights must be a matrix")
  }

  if (!is.null(initial_dist)) {
    if (!is.vector(initial_dist)) {
      stop("model$inits must be a vector")
    }
    if (length(initial_dist) != nrow(P)) {
      stop("Length of inits must match dimensions of weights matrix")
    }
  }

  # Print info
  if (verbose) {
    cat("\n===================================================================\n")
    cat("         ANALYZING TNA MODEL\n")
    cat("===================================================================\n")
    cat("Number of states:", nrow(P), "\n")

    if ("labels" %in% names(tna_model) && !is.null(tna_model$labels)) {
      cat("State labels:", paste(tna_model$labels, collapse = ", "), "\n")
    }
    cat("\n")
  }

  # Perform analysis
  results <- analyze_markov_chain(
    x = tna_model,
    initial_distribution = initial_dist,
    num_steps = num_steps,
    verbose = verbose,
    calculate_all = calculate_all,
    tol = tol
  )

  # Generate plots if requested
  if (plot_results) {
    if (verbose) {
      cat("\n===================================================================\n")
      cat("         GENERATING PLOTS\n")
      cat("===================================================================\n\n")
    }

    plot_all(results, ask = interactive())
  }

  return(results)
}

#' Summary method for markov_analysis objects
#'
#' @param object A markov_analysis object
#' @param ... Additional arguments (unused)
#' @export
summary.markov_analysis <- function(object, ...) {
  cat("\nMarkov Chain Summary\n")
  cat("====================\n\n")

  cat("States:", object$n_states, "\n")
  cat("Irreducible:", object$is_irreducible, "\n")
  cat("Aperiodic:", object$period == 1, "\n")
  cat("Regular:", object$is_regular, "\n")

  if (object$stationary_found) {
    cat("\nStationary Distribution:\n")
    print(round(object$equilibrium_distribution, 4))

    cat("\nConvergence Rate (Spectral Gap):", round(object$spectral_gap, 4), "\n")

    if (!is.na(object$mixing_time)) {
      cat("Mixing Time:", object$mixing_time, "steps\n")
    }
  }

  if (object$absorbing_chain_analysis$is_absorbing_chain) {
    cat("\nAbsorbing Chain Properties:\n")
    cat("Absorbing States:", paste(object$absorbing_chain_analysis$absorbing_states_names, collapse = ", "), "\n")
    cat("Transient States:", paste(object$absorbing_chain_analysis$transient_states_names, collapse = ", "), "\n")
  }

  cat("\nCommunicating Classes:", length(object$communicating_classes), "\n")

  invisible(object)
}

# ============================================================================
# EXAMPLE USAGE
# ============================================================================

#' Create Example TNA Model (Student Engagement)
#'
#' Creates a realistic TNA model object representing student engagement patterns
#' over time. The model has the same structure as actual TNA package output.
#'
#' @param n_students Integer. Number of student sequences to simulate. Default is 100.
#' @param n_timepoints Integer. Number of time points per sequence. Default is 20.
#' @param seed Integer. Random seed for reproducibility. Default is 123.
#'
#' @return A list with class "tna" containing:
#'   \describe{
#'     \item{weights}{3x3 transition probability matrix}
#'     \item{inits}{Named vector of initial state probabilities}
#'     \item{labels}{Character vector of state names}
#'     \item{data}{Matrix of sequences (numeric state values 1-3)}
#'   }
#'
#' @details
#' The model represents three engagement levels:
#' \itemize{
#'   \item \strong{Active}: Highly engaged students (state 1)
#'   \item \strong{Average}: Moderately engaged students (state 2)
#'   \item \strong{Disengaged}: Low engagement students (state 3)
#' }
#'
#' Transition patterns reflect realistic educational dynamics:
#' \itemize{
#'   \item Active students tend to stay active (high persistence)
#'   \item Average students can move in either direction
#'   \item Disengaged students slowly re-engage
#' }
#'
#' @examples
#' # Create default example model
#' model <- create_example_tna_model()
#' print(str(model))
#'
#' # Analyze the model
#' results <- analyze_tna_model(model, num_steps = 50)
#'
#' @export
create_example_tna_model <- function(n_students = 100, n_timepoints = 20, seed = 123) {
  set.seed(seed)

  # Realistic engagement transition matrix
  weights <- matrix(c(
    0.85, 0.12, 0.03,   # Active: mostly stays active
    0.25, 0.60, 0.15,   # Average: can go either way
    0.10, 0.30, 0.60    # Disengaged: slowly improves
  ), nrow = 3, byrow = TRUE)

  labels <- c("Active", "Average", "Disengaged")
  dimnames(weights) <- list(labels, labels)

  # Empirical initial distribution (from first observation)
  inits <- c(Active = 0.35, Average = 0.35, Disengaged = 0.30)

  # Simulate sequence data
  data_matrix <- matrix(NA, nrow = n_students, ncol = n_timepoints)
  colnames(data_matrix) <- paste0("T", seq_len(n_timepoints))

  for (i in seq_len(n_students)) {
    # Sample initial state
    current_state <- sample(1:3, 1, prob = inits)
    data_matrix[i, 1] <- current_state

    # Simulate trajectory
    for (t in 2:n_timepoints) {
      current_state <- sample(1:3, 1, prob = weights[current_state, ])
      data_matrix[i, t] <- current_state
    }
  }

  # Create TNA model object
  model <- list(
    weights = weights,
    inits = inits,
    labels = labels,
    data = data_matrix
  )

  class(model) <- "tna"
  attr(model, "type") <- "relative"

  return(model)
}

#' Create Example Absorbing Chain TNA Model (Course Completion)
#'
#' Creates a realistic TNA model with absorbing states representing student
#' progression through a course to either completion or dropout.
#'
#' @param n_students Integer. Number of student sequences to simulate. Default is 100.
#' @param n_timepoints Integer. Number of time points per sequence. Default is 20.
#' @param seed Integer. Random seed for reproducibility. Default is 456.
#'
#' @return A list with class "tna" containing:
#'   \describe{
#'     \item{weights}{4x4 transition probability matrix with absorbing states}
#'     \item{inits}{Named vector of initial state probabilities}
#'     \item{labels}{Character vector of state names}
#'     \item{data}{Matrix of sequences (numeric state values 1-4)}
#'   }
#'
#' @details
#' The model represents four course states:
#' \itemize{
#'   \item \strong{Struggling}: Students having difficulty (transient state 1)
#'   \item \strong{Progressing}: Students making steady progress (transient state 2)
#'   \item \strong{Completed}: Successfully finished course (absorbing state 3)
#'   \item \strong{Dropped}: Left the course (absorbing state 4)
#' }
#'
#' Once students reach "Completed" or "Dropped", they remain in that state.
#'
#' @examples
#' # Create example absorbing chain
#' model <- create_example_absorbing_chain()
#' print(str(model))
#'
#' # Analyze absorption properties
#' results <- analyze_tna_model(model, num_steps = 30)
#' results$absorbing_chain_analysis
#'
#' @export
create_example_absorbing_chain <- function(n_students = 100, n_timepoints = 20, seed = 456) {
  set.seed(seed)

  # Absorbing chain: Struggling, Progressing, Completed, Dropped
  weights <- matrix(c(
    0.50, 0.35, 0.05, 0.10,   # Struggling: may improve, complete, or drop
    0.15, 0.55, 0.25, 0.05,   # Progressing: likely to complete
    0.00, 0.00, 1.00, 0.00,   # Completed: absorbing state
    0.00, 0.00, 0.00, 1.00    # Dropped: absorbing state
  ), nrow = 4, byrow = TRUE)

  labels <- c("Struggling", "Progressing", "Completed", "Dropped")
  dimnames(weights) <- list(labels, labels)

  # Initial distribution (only transient states)
  inits <- c(Struggling = 0.40, Progressing = 0.60, Completed = 0.00, Dropped = 0.00)

  # Simulate sequence data
  data_matrix <- matrix(NA, nrow = n_students, ncol = n_timepoints)
  colnames(data_matrix) <- paste0("T", seq_len(n_timepoints))

  for (i in seq_len(n_students)) {
    # Sample initial state (only transient states)
    current_state <- sample(1:2, 1, prob = c(0.40, 0.60))
    data_matrix[i, 1] <- current_state

    # Simulate trajectory
    for (t in 2:n_timepoints) {
      current_state <- sample(1:4, 1, prob = weights[current_state, ])
      data_matrix[i, t] <- current_state
      # Once absorbed, stay in that state
      if (current_state %in% c(3, 4)) {
        data_matrix[i, t:n_timepoints] <- current_state
        break
      }
    }
  }

  # Create TNA model object
  model <- list(
    weights = weights,
    inits = inits,
    labels = labels,
    data = data_matrix
  )

  class(model) <- "tna"
  attr(model, "type") <- "relative"

  return(model)
}

# ============================================================================
# STATE STABILITY AND PERSISTENCE ANALYSIS
# ============================================================================

#' Calculate state persistence (self-transition probability)
#'
#' @param P Transition matrix or analysis object
#' @return Named vector of self-transition probabilities
#' @export
calculate_state_persistence <- function(P) {
  P <- extract_transition_matrix(P)
  persistence <- diag(P)
  if (!is.null(colnames(P))) {
    names(persistence) <- colnames(P)
  }
  return(persistence)
}

#' Calculate State Stability Index
#'
#' Computes comprehensive stability metrics for each state in a Markov chain,
#' combining persistence, recurrence, and stationary occupation measures.
#'
#' @param x Input object. Can be:
#'   \itemize{
#'     \item A transition matrix
#'     \item A markov_analysis object from \code{\link{analyze_markov_chain}}
#'   }
#'
#' @return A data frame with one row per state and the following columns:
#'   \describe{
#'     \item{State}{State name}
#'     \item{Persistence}{Self-transition probability P(i,i)}
#'     \item{Sojourn_Time}{Expected duration of stay per visit}
#'     \item{Recurrence_Time}{Mean time to return to state}
#'     \item{Stationary_Prob}{Long-run occupation probability}
#'     \item{Stability_Score}{Composite stability metric (higher = more stable)}
#'   }
#'
#' @details
#' State stability is assessed through multiple complementary metrics:
#'
#' \strong{Persistence:} The probability P(i,i) of remaining in the same state.
#' Higher values indicate states that resist leaving.
#'
#' \strong{Sojourn Time:} Expected number of consecutive time steps spent in
#' a state per visit, computed as 1/(1 - P(i,i)).
#'
#' \strong{Recurrence Time:} Mean number of steps to return to a state after
#' leaving it, computed as 1/pi(i) where pi is the stationary distribution.
#'
#' \strong{Stationary Probability:} Long-run proportion of time spent in the
#' state under equilibrium.
#'
#' \strong{Stability Score:} A composite metric defined as:
#' \deqn{Score = \frac{P(i,i) \times \pi(i)}{\sqrt{RecurrenceTime}}}
#'
#' This rewards states that persist, are frequently visited, and return quickly.
#'
#' @examples
#' # Basic example
#' P <- matrix(c(0.9, 0.1, 0.0,
#'               0.1, 0.7, 0.2,
#'               0.0, 0.3, 0.7), nrow = 3, byrow = TRUE)
#' colnames(P) <- c("Stable", "Medium", "Unstable")
#'
#' stability <- calculate_state_stability(P)
#' print(stability)
#'
#' @seealso
#' \code{\link{calculate_state_volatility}} for complementary volatility metrics,
#' \code{\link{stability_report}} for comprehensive stability analysis
#'
#' @export
calculate_state_stability <- function(x) {
  if (inherits(x, "markov_analysis")) {
    analysis <- x
    P <- analysis$transition_matrix
  } else {
    P <- extract_transition_matrix(x)
    analysis <- analyze_markov_chain(P, verbose = FALSE, calculate_all = FALSE)
  }

  state_names <- analysis$state_names
  n_states <- analysis$n_states

  # Self-transition probabilities (persistence)
  persistence <- diag(P)

  # Mean recurrence times (from MFPT diagonal)
  recurrence_times <- diag(analysis$mean_first_passage_times)

  # Stationary probabilities
  stationary_prob <- analysis$equilibrium_distribution

  # Expected sojourn time (how long you stay when you visit)
  sojourn <- analysis$sojourn_times

  # Stability score: composite measure (higher = more stable)
  # Formula: (persistence * stationary_prob) / sqrt(recurrence_time)
  # This rewards states that: persist, are visited often, and return quickly
  stability_score <- rep(NA, n_states)
  if (analysis$stationary_found) {
    stability_score <- (persistence * stationary_prob) / sqrt(pmax(recurrence_times, 1))
  }

  result <- data.frame(
    State = state_names,
    Persistence = round(persistence, 4),
    Sojourn_Time = round(sojourn, 2),
    Recurrence_Time = round(recurrence_times, 2),
    Stationary_Prob = round(stationary_prob, 4),
    Stability_Score = round(stability_score, 4),
    stringsAsFactors = FALSE
  )

  return(result)
}

#' Calculate state volatility (transition spread)
#'
#' Measures how dispersed the outgoing transitions are from each state.
#' Low entropy = concentrated transitions (stable), High entropy = dispersed (volatile)
#'
#' @param P Transition matrix or analysis object
#' @param tol Tolerance for zero probabilities
#' @return Data frame with volatility metrics
#' @export
calculate_state_volatility <- function(P, tol = 1e-10) {
  P <- extract_transition_matrix(P)
  state_names <- colnames(P)
  if (is.null(state_names)) {
    state_names <- paste0("S", seq_len(nrow(P)))
  }

  n_states <- nrow(P)
  entropy <- numeric(n_states)
  gini <- numeric(n_states)
  exit_probability <- numeric(n_states)
  effective_destinations <- numeric(n_states)

  for (i in 1:n_states) {
    row <- P[i, ]

    # Entropy of outgoing transitions (uncertainty in next state)
    entropy[i] <- calculate_entropy(row, base = 2, tol = tol)

    # Gini impurity (how concentrated the transitions are)
    gini[i] <- calculate_gini_impurity(row)

    # Exit probability (1 - self-loop)
    exit_probability[i] <- 1 - P[i, i]

    # Effective number of destinations (exponential of entropy)
    effective_destinations[i] <- 2^entropy[i]
  }

  result <- data.frame(
    State = state_names,
    Exit_Probability = round(exit_probability, 4),
    Transition_Entropy = round(entropy, 4),
    Transition_Gini = round(gini, 4),
    Effective_Destinations = round(effective_destinations, 2),
    stringsAsFactors = FALSE
  )

  return(result)
}

#' Calculate expected number of visits to each state
#'
#' Before absorption or over a finite time horizon, how many times
#' do we expect to visit each state?
#'
#' @param P Transition matrix or analysis object
#' @param initial_state Initial state or distribution
#' @param n_steps Time horizon (if NULL, uses fundamental matrix for infinite horizon)
#' @return Named vector of expected visits
#' @export
calculate_expected_visits <- function(P, initial_state = 1, n_steps = NULL) {
  P <- extract_transition_matrix(P)
  n_states <- nrow(P)
  state_names <- colnames(P)

  # Handle initial state
  if (length(initial_state) == 1) {
    if (is.character(initial_state)) {
      idx <- which(state_names == initial_state)
      if (length(idx) == 0) stop("State name not found")
      initial_state <- idx
    }
    initial_dist <- rep(0, n_states)
    initial_dist[initial_state] <- 1
  } else {
    initial_dist <- initial_state
  }

  if (is.null(n_steps)) {
    # Infinite horizon using fundamental matrix
    expected_visits <- tryCatch({
      I <- diag(n_states)
      # For ergodic chains, expected visits = (I - P + 1*pi')^-1
      analysis <- analyze_markov_chain(P, verbose = FALSE, calculate_all = FALSE)
      if (analysis$stationary_found) {
        eq_dist <- analysis$equilibrium_distribution
        Z <- solve(I - P + matrix(eq_dist, n_states, n_states, byrow = TRUE))
        as.vector(initial_dist %*% Z)
      } else {
        warning("Cannot calculate infinite horizon visits without stationary distribution",
                call. = FALSE)
        rep(NA_real_, n_states)
      }
    }, error = function(e) {
      warning("Could not calculate expected visits: ", e$message, call. = FALSE)
      rep(NA_real_, n_states)
    })
  } else {
    # Finite horizon: sum of occupation probabilities
    current_dist <- initial_dist
    occupation <- matrix(0, nrow = n_steps + 1, ncol = n_states)

    for (t in 0:n_steps) {
      occupation[t + 1, ] <- current_dist
      if (t < n_steps) {
        current_dist <- as.vector(current_dist %*% P)
      }
    }
    expected_visits <- colSums(occupation)
  }

  names(expected_visits) <- state_names
  return(expected_visits)
}

#' Calculate state transition resilience
#'
#' Measures how resistant each state is to leaving after perturbations.
#' Uses the sensitivity of staying probability to small changes in transition matrix.
#'
#' @param P Transition matrix or analysis object
#' @param epsilon Perturbation size
#' @return Data frame with resilience metrics
#' @export
calculate_state_resilience <- function(P, epsilon = 0.01) {
  P <- extract_transition_matrix(P)
  state_names <- colnames(P)
  n_states <- nrow(P)

  base_persistence <- diag(P)
  resilience <- numeric(n_states)

  # For each state, perturb outgoing transitions and measure change
  for (i in 1:n_states) {
    if (base_persistence[i] > epsilon) {
      # Create perturbed matrix: reduce self-loop, distribute to others
      P_perturbed <- P
      P_perturbed[i, i] <- P[i, i] - epsilon

      # Distribute epsilon proportionally to other states
      other_states <- setdiff(1:n_states, i)
      if (length(other_states) > 0) {
        distribution <- P[i, other_states] / sum(P[i, other_states])
        P_perturbed[i, other_states] <- P[i, other_states] + epsilon * distribution
      }

      # Calculate change in long-term occupation
      tryCatch({
        # Get stationary distributions
        eigen_base <- eigen(t(P))
        idx_base <- which.min(abs(eigen_base$values - 1))
        pi_base <- Re(eigen_base$vectors[, idx_base])
        pi_base <- pi_base / sum(pi_base)

        eigen_pert <- eigen(t(P_perturbed))
        idx_pert <- which.min(abs(eigen_pert$values - 1))
        pi_pert <- Re(eigen_pert$vectors[, idx_pert])
        pi_pert <- pi_pert / sum(pi_pert)

        # Resilience: smaller change = more resilient
        resilience[i] <- abs(pi_base[i] - pi_pert[i]) / epsilon
      }, error = function(e) {
        resilience[i] <- NA
      })
    } else {
      resilience[i] <- NA  # State has very low self-transition
    }
  }

  result <- data.frame(
    State = state_names,
    Base_Persistence = round(base_persistence, 4),
    Resilience_Index = round(1 / (1 + resilience), 4),  # Transform so higher = more resilient
    Sensitivity = round(resilience, 4),
    stringsAsFactors = FALSE
  )

  return(result)
}

#' Plot state stability analysis
#'
#' @param x Transition matrix or markov_analysis object
#' @param type Type of plot ("stability", "volatility", "persistence", "occupation")
#' @param state State to focus on (for persistence and occupation plots)
#' @param n_steps Number of steps for time-based plots
#' @export
plot_stability_analysis <- function(x, type = "stability", state = 1, n_steps = 100) {

  if (type == "stability") {
    stability <- calculate_state_stability(x)

    graphics::par(mfrow = c(2, 2))

    # Plot 1: Persistence
    graphics::barplot(stability$Persistence, names.arg = stability$State,
            main = "State Persistence (Self-Transition Probability)",
            ylab = "P(stay in state)", col = "steelblue", las = 2)
    graphics::abline(h = 0.5, lty = 2, col = "red")

    # Plot 2: Sojourn Time
    graphics::barplot(stability$Sojourn_Time, names.arg = stability$State,
            main = "Expected Sojourn Time",
            ylab = "Time steps", col = "darkgreen", las = 2)

    # Plot 3: Stationary Probability
    graphics::barplot(stability$Stationary_Prob, names.arg = stability$State,
            main = "Long-term Occupation (Stationary Distribution)",
            ylab = "Probability", col = "purple", las = 2)

    # Plot 4: Stability Score
    graphics::barplot(stability$Stability_Score, names.arg = stability$State,
            main = "Overall Stability Score",
            ylab = "Score (higher = more stable)", col = "orange", las = 2)

    graphics::par(mfrow = c(1, 1))

  } else if (type == "volatility") {
    volatility <- calculate_state_volatility(x)

    graphics::par(mfrow = c(2, 2))

    graphics::barplot(volatility$Exit_Probability, names.arg = volatility$State,
            main = "Exit Probability (1 - Self-Loop)",
            ylab = "P(leave state)", col = "coral", las = 2)

    graphics::barplot(volatility$Transition_Entropy, names.arg = volatility$State,
            main = "Transition Entropy (Uncertainty)",
            ylab = "Bits", col = "steelblue", las = 2)

    graphics::barplot(volatility$Effective_Destinations, names.arg = volatility$State,
            main = "Effective Number of Destinations",
            ylab = "Count", col = "darkgreen", las = 2)

    graphics::barplot(volatility$Transition_Gini, names.arg = volatility$State,
            main = "Transition Gini Impurity",
            ylab = "Impurity", col = "purple", las = 2)

    graphics::par(mfrow = c(1, 1))

  } else if (type == "persistence") {
    P <- extract_transition_matrix(x)
    state_names <- colnames(P)

    if (is.character(state)) {
      state_idx <- which(state_names == state)
      if (length(state_idx) == 0) stop("State name not found")
    } else {
      state_idx <- state
    }

    # Calculate occupation probability over time
    current_dist <- rep(0, nrow(P))
    current_dist[state_idx] <- 1
    persistence <- numeric(n_steps + 1)

    for (t in 0:n_steps) {
      persistence[t + 1] <- current_dist[state_idx]
      if (t < n_steps) {
        current_dist <- as.vector(current_dist %*% P)
      }
    }

    graphics::plot(0:n_steps, persistence, type = 'l', lwd = 2, col = 'blue',
         xlab = "Time Step", ylab = "P(still in starting state)",
         main = paste("State Persistence Trajectory for State", state),
         panel.first = graphics::grid())

    # Add exponential decay reference
    if (length(persistence) > 1 && persistence[2] > 0) {
      lambda <- -log(persistence[2])
      exp_decay <- exp(-lambda * (0:n_steps))
      graphics::lines(0:n_steps, exp_decay, lty = 2, col = 'red', lwd = 2)
      graphics::legend("topright", legend = c("Actual", "Exponential Decay"),
             col = c("blue", "red"), lty = c(1, 2), lwd = 2)
    }

  } else if (type == "occupation") {
    P <- extract_transition_matrix(x)
    state_names <- colnames(P)

    if (is.character(state)) {
      state_idx <- which(state_names == state)
      if (length(state_idx) == 0) stop("State name not found")
    } else {
      state_idx <- state
    }

    # Calculate occupation probability over time
    n_states <- nrow(P)
    current_dist <- rep(0, n_states)
    current_dist[state_idx] <- 1
    occupation <- matrix(0, nrow = n_steps + 1, ncol = n_states)
    colnames(occupation) <- state_names

    for (t in 0:n_steps) {
      occupation[t + 1, ] <- current_dist
      if (t < n_steps) {
        current_dist <- as.vector(current_dist %*% P)
      }
    }

    colors <- grDevices::rainbow(n_states)

    graphics::plot(0:n_steps, occupation[, 1], type = 'l', lwd = 2, col = colors[1],
         xlab = "Time Step", ylab = "Occupation Probability",
         main = paste("State Occupation Over Time (Starting from State", state, ")"),
         ylim = c(0, 1), panel.first = graphics::grid())

    if (n_states > 1) {
      for (i in 2:n_states) {
        graphics::lines(0:n_steps, occupation[, i], lwd = 2, col = colors[i])
      }
    }

    graphics::legend("topright", legend = colnames(occupation), col = colors, lwd = 2, cex = 0.7)
  }
}

#' Comprehensive stability report
#'
#' @param x Transition matrix or markov_analysis object
#' @param print_report Print detailed report
#' @return List with all stability metrics
#' @export
stability_report <- function(x, print_report = TRUE) {
  stability <- calculate_state_stability(x)
  volatility <- calculate_state_volatility(x)
  resilience <- calculate_state_resilience(x)

  # Merge all metrics
  report <- merge(stability, volatility, by = "State")
  report <- merge(report, resilience, by = "State")

  if (print_report) {
    cat("\n===================================================================\n")
    cat("           STATE STABILITY REPORT\n")
    cat("===================================================================\n\n")

    cat("-- MOST STABLE STATES (by Stability Score)\n")
    top_stable <- utils::head(report[order(-report$Stability_Score), ], 3)
    print(top_stable[, c("State", "Persistence", "Stability_Score", "Sojourn_Time")], row.names = FALSE)

    cat("\n-- MOST VOLATILE STATES (by Transition Entropy)\n")
    top_volatile <- utils::head(report[order(-report$Transition_Entropy), ], 3)
    print(top_volatile[, c("State", "Exit_Probability", "Transition_Entropy", "Effective_Destinations")], row.names = FALSE)

    cat("\n-- MOST RESILIENT STATES\n")
    top_resilient <- utils::head(report[order(-report$Resilience_Index), ], 3)
    print(top_resilient[, c("State", "Base_Persistence", "Resilience_Index")], row.names = FALSE)

    cat("\n===================================================================\n")
  }

  invisible(report)
}

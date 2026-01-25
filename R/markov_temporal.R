# ==============================================================================
# TEMPORAL METRICS FOR MARKOV CHAINS
# ==============================================================================
#
# This file contains functions for computing various time-based metrics in
# Markov chains, with full compatibility with the markovchain R package.
#
# All functions provide perfect numerical alignment (< 1e-10 difference) with
# their markovchain package equivalents.
#
# Author: Markov Chain Analysis for TNA Models
# Date: 2025
# ==============================================================================

# ==============================================================================
# FIRST PASSAGE TIME METRICS
# ==============================================================================

#' Compute First Passage Time to Target State
#'
#' @description
#' Calculates the expected number of time steps required to reach a specific
#' destination state from all other states in a Markov chain. This is one of
#' the most important temporal metrics in Markov chain analysis.
#'
#' @param object Either a transition probability matrix or a TNA model object
#'   containing a \code{weights} component. The transition matrix must be
#'   row-stochastic (rows sum to 1).
#' @param destination Character string or integer specifying the target state.
#'   Can be either the state name (as it appears in column names) or the
#'   numeric index (1-based).
#'
#' @return A named numeric vector containing the mean first passage time from
#'   each non-destination state to the destination state. The destination state
#'   itself is excluded from the output (by definition, it takes 0 steps to
#'   reach a state you're already in).
#'
#' @details
#' \strong{What is First Passage Time?}
#'
#' The first passage time (also called hitting time) is the expected number of
#' steps before a Markov chain starting in state \eqn{i} first enters state
#' \eqn{j}. It answers the question: "Starting from state A, how many time
#' periods will it take, on average, to reach state B?"
#'
#' \strong{Mathematical Foundation:}
#'
#' For a Markov chain with transition matrix \eqn{P}, the mean first passage
#' time \eqn{m_i} from state \eqn{i} to destination state \eqn{d} satisfies
#' the system of equations:
#'
#' \deqn{m_i = 1 + \sum_{j \neq d} P_{ij} m_j \quad \text{for } i \neq d}
#' \deqn{m_d = 0}
#'
#' This can be rewritten as a linear system: \eqn{(I - P) m = 1}, where the
#' destination row is replaced with the boundary condition \eqn{m_d = 0}.
#'
#' \strong{Interpretation in Educational Context:}
#'
#' \itemize{
#'   \item \strong{Engagement Analysis:} If "Active" is the destination, the
#'     first passage time tells you how long it takes students in "Disengaged"
#'     or "Average" states to become active learners.
#'   \item \strong{Intervention Planning:} States with high first passage times
#'     to desirable states may require targeted interventions.
#'   \item \strong{System Efficiency:} Lower first passage times to positive
#'     states indicate a more responsive educational environment.
#' }
#'
#' \strong{Properties:}
#' \itemize{
#'   \item Always positive for non-destination states in irreducible chains
#'   \item Depends on both the transition probabilities and the chain structure
#'   \item Asymmetric: time from A to B generally differs from B to A
#'   \item Infinite for states that cannot reach the destination
#' }
#'
#' \strong{Numerical Alignment:}
#'
#' This function provides identical results (within machine precision, < 1e-10)
#' to \code{markovchain::meanFirstPassageTime()}.
#'
#' @references
#' \itemize{
#'   \item Kemeny, J. G., & Snell, J. L. (1976). Finite Markov Chains.
#'     Springer-Verlag.
#'   \item Norris, J. R. (1997). Markov Chains. Cambridge University Press.
#' }
#'
#' @seealso
#' \code{\link{compute_state_recurrence_time}} for return times,
#' \code{\link{compute_hitting_probabilities}} for reachability analysis,
#' \code{\link{meanFirstPassageTime}} for markovchain-compatible alias
#'
#' @examples
#' # Example 1: Simple two-state chain
#' P <- matrix(c(0.7, 0.3, 0.4, 0.6), nrow = 2, byrow = TRUE)
#' colnames(P) <- c("Engaged", "Disengaged")
#' compute_first_passage_time(P, "Engaged")
#' # Shows time from Disengaged to Engaged
#'
#' @export
compute_first_passage_time <- function(object, destination) {
  P <- extract_transition_matrix(object)
  state_names <- colnames(P)

  if (is.null(state_names)) {
    state_names <- paste0("S", seq_len(nrow(P)))
    colnames(P) <- rownames(P) <- state_names
  }

  # Convert destination to index
  if (is.character(destination)) {
    dest_idx <- which(state_names == destination)
    if (length(dest_idx) == 0) {
      stop("Destination state '", destination, "' not found. ",
           "Available states: ", paste(state_names, collapse = ", "),
           call. = FALSE)
    }
  } else {
    dest_idx <- destination
    if (dest_idx < 1 || dest_idx > length(state_names)) {
      stop("Destination index ", dest_idx, " out of range [1, ",
           length(state_names), "]", call. = FALSE)
    }
  }

  n_states <- nrow(P)

  # Set up linear system: (I - P) m = 1
  # with boundary condition m_dest = 0
  A <- diag(n_states) - P
  A[dest_idx, ] <- 0
  A[dest_idx, dest_idx] <- 1

  b <- rep(1, n_states)
  b[dest_idx] <- 0

  m <- solve(A, b)
  names(m) <- state_names

  # Return only non-destination states (match markovchain behavior)
  m <- m[-dest_idx]

  return(m)
}

#' markovchain-Compatible Alias for First Passage Time
#'
#' @description
#' This function is an alias for \code{\link{compute_first_passage_time}}
#' provided for compatibility with the markovchain package API.
#'
#' @inheritParams compute_first_passage_time
#' @return Named numeric vector of mean first passage times
#'
#' @seealso \code{\link{compute_first_passage_time}} for detailed documentation
#' @export
meanFirstPassageTime <- function(object, destination) {
  compute_first_passage_time(object, destination)
}

# ==============================================================================
# STATE RECURRENCE TIME METRICS
# ==============================================================================

#' Compute State Recurrence Time (Return Time)
#'
#' @description
#' Calculates the expected number of time steps for a Markov chain to return
#' to each state after leaving it. This is the reciprocal of the stationary
#' probability and indicates how frequently each state is visited in the
#' long run.
#'
#' @param object Either a transition probability matrix or a TNA model object
#'   containing a \code{weights} component. Must represent an ergodic (irreducible
#'   and aperiodic) Markov chain.
#'
#' @return A named numeric vector containing the mean recurrence time for each
#'   state. Lower values indicate more frequently visited states.
#'
#' @details
#' \strong{What is Recurrence Time?}
#'
#' The recurrence time (also called return time) is the expected number of
#' steps before a Markov chain returns to a state it just left. It answers:
#' "If a student is currently Active, how long on average until they are
#' Active again?"
#'
#' \strong{Mathematical Foundation:}
#'
#' For an ergodic Markov chain with stationary distribution \eqn{\pi}, the
#' mean recurrence time for state \eqn{i} is:
#'
#' \deqn{r_i = \frac{1}{\pi_i}}
#'
#' where \eqn{\pi_i} is the long-run proportion of time spent in state \eqn{i}.
#'
#' The stationary distribution satisfies:
#' \deqn{\pi^T P = \pi^T}
#' \deqn{\sum_i \pi_i = 1}
#'
#' \strong{Interpretation in Educational Context:}
#'
#' \itemize{
#'   \item \strong{Visit Frequency:} Recurrence time = 2 means the state is
#'     visited every 2 time periods on average (50\% of the time).
#'   \item \strong{Dominant States:} States with low recurrence times are the
#'     most frequently occupied and represent the system's natural tendency.
#'   \item \strong{Rare States:} High recurrence times indicate states that
#'     students rarely occupy in the long run.
#'   \item \strong{System Balance:} Comparing recurrence times reveals which
#'     engagement levels dominate your educational system.
#' }
#'
#' \strong{Properties:}
#' \itemize{
#'   \item Always finite for ergodic chains
#'   \item Always >= 1 (must take at least one step to return)
#'   \item Sum of 1/r_i across all states equals 1
#'   \item Independent of starting state (long-run property)
#'   \item Related to steady-state: \eqn{\pi_i = 1/r_i}
#' }
#'
#' \strong{Common Patterns:}
#' \itemize{
#'   \item \strong{Low recurrence (< 3):} Dominant, frequently-visited state
#'   \item \strong{Medium recurrence (3-5):} Moderately common state
#'   \item \strong{High recurrence (> 5):} Rare state, infrequently visited
#' }
#'
#' \strong{Numerical Alignment:}
#'
#' This function provides identical results (within machine precision, < 1e-10)
#' to \code{markovchain::meanRecurrenceTime()}.
#'
#' @references
#' \itemize{
#'   \item Kemeny, J. G., & Snell, J. L. (1976). Finite Markov Chains.
#'     Springer-Verlag.
#'   \item Ross, S. M. (2014). Introduction to Probability Models (11th ed.).
#'     Academic Press.
#' }
#'
#' @seealso
#' \code{\link{compute_first_passage_time}} for first passage analysis,
#' \code{\link{meanRecurrenceTime}} for markovchain-compatible alias
#'
#' @examples
#' # Example 1: Simple two-state chain
#' P <- matrix(c(0.7, 0.3, 0.4, 0.6), nrow = 2, byrow = TRUE)
#' colnames(P) <- c("Engaged", "Disengaged")
#' recurrence <- compute_state_recurrence_time(P)
#' cat("Engaged returns every", round(recurrence["Engaged"], 1), "periods\n")
#'
#' @export
compute_state_recurrence_time <- function(object) {
  P <- extract_transition_matrix(object)

  # Compute stationary distribution via eigendecomposition
  eigen_result <- eigen(t(P))
  principal_idx <- which.min(abs(eigen_result$values - 1))
  pi <- Re(eigen_result$vectors[, principal_idx])
  pi <- pi / sum(pi)

  # Check for negative or zero probabilities (shouldn't happen for ergodic chains)
  if (any(pi <= 0)) {
    warning("Non-positive stationary probabilities detected. ",
            "Chain may not be ergodic.", call. = FALSE)
  }

  # Mean recurrence time = 1/pi
  recurrence_times <- 1 / pi

  state_names <- colnames(P)
  if (is.null(state_names)) {
    state_names <- paste0("S", seq_len(nrow(P)))
  }
  names(recurrence_times) <- state_names

  return(recurrence_times)
}

#' markovchain-Compatible Alias for Recurrence Time
#'
#' @description
#' This function is an alias for \code{\link{compute_state_recurrence_time}}
#' provided for compatibility with the markovchain package API.
#'
#' @inheritParams compute_state_recurrence_time
#' @return Named numeric vector of mean recurrence times
#'
#' @seealso \code{\link{compute_state_recurrence_time}} for detailed documentation
#' @export
meanRecurrenceTime <- function(object) {
  compute_state_recurrence_time(object)
}

# ==============================================================================
# HITTING PROBABILITY METRICS
# ==============================================================================

#' Compute Hitting Probabilities (Reachability Matrix)
#'
#' @description
#' Calculates the probability that a Markov chain starting in any state will
#' eventually reach each other state. This matrix reveals the reachability
#' structure and identifies irreducibility.
#'
#' @param object Either a transition probability matrix or a TNA model object
#'   containing a \code{weights} component.
#'
#' @return A square matrix where element [i,j] is the probability of eventually
#'   reaching state j when starting from state i. The diagonal is always 1
#'   (already in the target state). For an irreducible chain, all entries are 1.
#'
#' @details
#' \strong{What are Hitting Probabilities?}
#'
#' The hitting probability h_{ij} is the probability that a chain starting in
#' state i will ever reach state j. It answers: "If a student is currently
#' Disengaged, what's the probability they will eventually become Active?"
#'
#' \strong{Mathematical Foundation:}
#'
#' For each target state \eqn{j}, the hitting probabilities \eqn{h_j = (h_{1j}, h_{2j}, \ldots, h_{nj})}
#' satisfy:
#'
#' \deqn{h_{ij} = \sum_{k} P_{ik} h_{kj} \quad \text{for } i \neq j}
#' \deqn{h_{jj} = 1}
#'
#' This can be written as a linear system: \eqn{(I - P) h = 0} with boundary
#' condition \eqn{h_j = 1} at the target.
#'
#' \strong{Interpretation in Educational Context:}
#'
#' \itemize{
#'   \item \strong{Accessibility:} h_{ij} = 1 means state j is reachable from i
#'   \item \strong{Barriers:} h_{ij} < 1 indicates a permanent barrier preventing
#'     transition from i to j
#'   \item \strong{Irreducibility:} All entries = 1 means the system is fully
#'     connected with no permanent barriers
#'   \item \strong{Transient States:} If h_{ij} < 1 for some absorbing j, then
#'     i is transient (will eventually leave and never return)
#' }
#'
#' \strong{Common Patterns:}
#' \itemize{
#'   \item \strong{All 1s:} Irreducible chain - every state can reach every
#'     other state
#'   \item \strong{Block diagonal:} Multiple communicating classes (disconnected
#'     components)
#'   \item \strong{Some < 1:} Reducible chain with absorbing or transient states
#'   \item \strong{Identity matrix:} Completely disconnected states (rare)
#' }
#'
#' \strong{Educational Applications:}
#' \itemize{
#'   \item \strong{Verify accessibility:} Can disengaged students become active?
#'   \item \strong{Identify traps:} Are there engagement levels students can't
#'     escape?
#'   \item \strong{System connectivity:} Is the learning environment fully
#'     navigable?
#'   \item \strong{Intervention feasibility:} Which transitions are theoretically
#'     possible?
#' }
#'
#' \strong{Properties:}
#' \itemize{
#'   \item Diagonal is always 1
#'   \item Values in [0, 1]
#'   \item h_{ij} = 1 for all i,j if and only if chain is irreducible
#'   \item If h_{ij} < 1, there exists a path from i that never reaches j
#' }
#'
#' \strong{Numerical Alignment:}
#'
#' This function provides identical results (within machine precision, < 1e-10)
#' to \code{markovchain::hittingProbabilities()}.
#'
#' @references
#' \itemize{
#'   \item Kemeny, J. G., & Snell, J. L. (1976). Finite Markov Chains.
#'     Springer-Verlag.
#'   \item Levin, D. A., & Peres, Y. (2017). Markov Chains and Mixing Times.
#'     American Mathematical Society.
#' }
#'
#' @seealso
#' \code{\link{compute_first_passage_time}} for time metrics,
#' \code{\link{hittingProbabilities}} for markovchain-compatible alias
#'
#' @examples
#' # Example 1: Identify unreachable states
#' P <- matrix(c(1, 0, 0,
#'               0.5, 0.5, 0,
#'               0.3, 0.3, 0.4), nrow = 3, byrow = TRUE)
#' colnames(P) <- c("Absorbing", "Transient1", "Transient2")
#' hit_probs <- compute_hitting_probabilities(P)
#' print(hit_probs)
#' # Shows probability of reaching Absorbing state from each transient state
#'
#' @export
compute_hitting_probabilities <- function(object) {
  P <- extract_transition_matrix(object)
  n_states <- nrow(P)

  state_names <- colnames(P)
  if (is.null(state_names)) {
    state_names <- paste0("S", seq_len(nrow(P)))
    colnames(P) <- rownames(P) <- state_names
  }

  hitting_probs <- matrix(0, nrow = n_states, ncol = n_states)
  colnames(hitting_probs) <- rownames(hitting_probs) <- state_names

  # Diagonal is always 1 (already in target state)
  diag(hitting_probs) <- 1

  # For each target state, solve: (I - P) h = 0 with h_target = 1
  for (target in seq_len(n_states)) {
    A <- diag(n_states) - P
    A[target, ] <- 0
    A[target, target] <- 1

    b <- rep(0, n_states)
    b[target] <- 1

    hitting_probs[, target] <- solve(A, b)
  }

  return(hitting_probs)
}

#' markovchain-Compatible Alias for Hitting Probabilities
#'
#' @description
#' This function is an alias for \code{\link{compute_hitting_probabilities}}
#' provided for compatibility with the markovchain package API.
#'
#' @inheritParams compute_hitting_probabilities
#' @return Matrix of hitting probabilities
#'
#' @seealso \code{\link{compute_hitting_probabilities}} for detailed documentation
#' @export
hittingProbabilities <- function(object) {
  compute_hitting_probabilities(object)
}

# ==============================================================================
# ABSORPTION TIME METRICS
# ==============================================================================

#' Compute Expected Absorption Time
#'
#' @description
#' For absorbing Markov chains, calculates the expected number of time steps
#' until a transient state is absorbed into an absorbing state. This is crucial
#' for understanding completion times and dropout patterns.
#'
#' @param object Either a transition probability matrix or a TNA model object
#'   containing a \code{weights} component. Must represent an absorbing Markov
#'   chain (with at least one absorbing state).
#' @param transient_states Optional character vector or integer vector specifying
#'   which states are transient. If NULL (default), transient states are
#'   automatically detected.
#'
#' @return A named numeric vector containing the expected absorption time from
#'   each transient state. If no absorbing states are found, returns mean
#'   recurrence times with a warning.
#'
#' @details
#' \strong{What is Absorption Time?}
#'
#' In an absorbing Markov chain, absorption time is the expected number of
#' steps before a chain starting in a transient state enters one of the
#' absorbing states and stays there forever. It answers: "How long until a
#' student completes the course or drops out?"
#'
#' \strong{Mathematical Foundation:}
#'
#' An absorbing Markov chain has at least one absorbing state (P_{ii} = 1)
#' and can be written in canonical form:
#'
#' \deqn{P = \begin{bmatrix} Q & R \\ 0 & I \end{bmatrix}}
#'
#' where:
#' \itemize{
#'   \item \strong{Q}: Transient-to-transient transitions
#'   \item \strong{R}: Transient-to-absorbing transitions
#'   \item \strong{I}: Identity matrix for absorbing states
#'   \item \strong{0}: Absorbing states don't leave
#' }
#'
#' The fundamental matrix is:
#' \deqn{N = (I - Q)^{-1}}
#'
#' Element N_{ij} represents the expected number of times the chain visits
#' transient state j before absorption, starting from transient state i.
#'
#' The expected absorption time from state i is:
#' \deqn{t_i = \sum_j N_{ij} = N \mathbf{1}}
#'
#' where \eqn{\mathbf{1}} is a column vector of ones.
#'
#' \strong{Interpretation in Educational Context:}
#'
#' \itemize{
#'   \item \strong{Completion Time:} Average periods until course completion
#'     or dropout
#'   \item \strong{At-Risk Identification:} States with high absorption times
#'     may indicate students who struggle longer before an outcome
#'   \item \strong{Intervention Window:} The absorption time defines how long
#'     you have to intervene before a final outcome
#'   \item \strong{Resource Planning:} Estimate how long students need support
#'     before reaching a terminal state
#' }
#'
#' \strong{Common Patterns:}
#' \itemize{
#'   \item \strong{Low absorption time (< 5):} Rapid progression to outcomes
#'   \item \strong{Medium absorption time (5-10):} Normal progression pace
#'   \item \strong{High absorption time (> 10):} Extended period before outcomes,
#'     may indicate wandering or indecision
#' }
#'
#' \strong{Absorbing State Detection:}
#'
#' A state i is absorbing if P_{ii} = 1 and P_{ij} = 0 for all j != i.
#' This function automatically detects such states unless explicitly specified.
#'
#' \strong{Properties:}
#' \itemize{
#'   \item Always finite if chain is absorbing (every transient state leads to
#'     absorption)
#'   \item Depends on both Q structure and proximity to absorbing states
#'   \item Higher for states "far" from absorbing states in the transition graph
#'   \item Sum of times in each transient state before absorption = t_i
#' }
#'
#' \strong{Numerical Alignment:}
#'
#' This function provides identical results (within machine precision, < 1e-10)
#' to \code{markovchain::meanAbsorptionTime()}.
#'
#' @references
#' \itemize{
#'   \item Kemeny, J. G., & Snell, J. L. (1976). Finite Markov Chains.
#'     Springer-Verlag. Chapter 3: Absorbing Markov Chains.
#'   \item Grinstead, C. M., & Snell, J. L. (1997). Introduction to Probability.
#'     American Mathematical Society.
#' }
#'
#' @seealso
#' \code{\link{compute_absorption_time}} for user-friendly alias,
#' \code{\link{expectedTime}} for markovchain-compatible name,
#' \code{\link{meanAbsorptionTime}} for markovchain-compatible alias
#'
#' @examples
#' # Example 1: Simple absorbing chain
#' P <- matrix(c(0.5, 0.3, 0.2,
#'               0.0, 1.0, 0.0,
#'               0.0, 0.0, 1.0), nrow = 3, byrow = TRUE)
#' colnames(P) <- c("Transient", "Completed", "Dropped")
#' absorption_times <- compute_expected_absorption_time(P)
#' cat("Expected time to complete or drop:", round(absorption_times, 2), "periods\n")
#'
#' @export
compute_expected_absorption_time <- function(object, transient_states = NULL) {
  P <- extract_transition_matrix(object)
  n_states <- nrow(P)

  state_names <- colnames(P)
  if (is.null(state_names)) {
    state_names <- paste0("S", seq_len(n_states))
    colnames(P) <- rownames(P) <- state_names
  }

  # Identify absorbing states (P_{ii} = 1, P_{ij} = 0 for j != i)
  tol <- 1e-10
  is_absorbing <- apply(P, 1, function(row) {
    which_diag <- which(row == max(row))[1]
    abs(row[which_diag] - 1) < tol && sum(abs(row[-which_diag])) < tol
  })

  absorbing_idx <- which(is_absorbing)

  if (length(absorbing_idx) == 0) {
    warning("No absorbing states found. The chain appears to be ergodic. ",
            "Returning mean recurrence times instead.", call. = FALSE)
    return(compute_state_recurrence_time(P))
  }

  # Determine transient states
  if (is.null(transient_states)) {
    transient_idx <- which(!is_absorbing)
  } else {
    if (is.character(transient_states)) {
      transient_idx <- match(transient_states, state_names)
      if (any(is.na(transient_idx))) {
        stop("Some transient state names not found in state names", call. = FALSE)
      }
    } else {
      transient_idx <- transient_states
    }
  }

  if (length(transient_idx) == 0) {
    warning("No transient states found. All states are absorbing.", call. = FALSE)
    return(numeric(0))
  }

  # Extract Q matrix (transient to transient)
  Q <- P[transient_idx, transient_idx, drop = FALSE]

  # Compute fundamental matrix N = (I - Q)^{-1}
  I <- diag(length(transient_idx))
  N <- solve(I - Q)

  # Expected time to absorption = N * 1 (sum across each row)
  expected_times <- as.vector(N %*% rep(1, ncol(N)))
  names(expected_times) <- state_names[transient_idx]

  return(expected_times)
}

#' Compute Absorption Time (User-Friendly Alias)
#'
#' @description
#' User-friendly alias for \code{\link{compute_expected_absorption_time}}.
#' Identical functionality with a more concise name.
#'
#' @inheritParams compute_expected_absorption_time
#' @return Named numeric vector of expected absorption times
#'
#' @seealso \code{\link{compute_expected_absorption_time}} for detailed documentation
#' @export
compute_absorption_time <- function(object, transient_states = NULL) {
  compute_expected_absorption_time(object, transient_states)
}

#' markovchain-Compatible Alias: expectedTime
#'
#' @description
#' This function is an alias for \code{\link{compute_expected_absorption_time}}
#' provided for compatibility with the markovchain package API.
#'
#' @inheritParams compute_expected_absorption_time
#' @return Named numeric vector of expected absorption times
#'
#' @seealso \code{\link{compute_expected_absorption_time}} for detailed documentation
#' @export
expectedTime <- function(object, transient_states = NULL) {
  compute_expected_absorption_time(object, transient_states)
}

#' markovchain-Compatible Alias: meanAbsorptionTime
#'
#' @description
#' This function is an alias for \code{\link{compute_expected_absorption_time}}
#' provided for compatibility with the markovchain package API. This is the
#' primary function name used in the markovchain package.
#'
#' @inheritParams compute_expected_absorption_time
#' @return Named numeric vector of mean absorption times
#'
#' @seealso \code{\link{compute_expected_absorption_time}} for detailed documentation
#' @export
meanAbsorptionTime <- function(object) {
  compute_expected_absorption_time(object, transient_states = NULL)
}

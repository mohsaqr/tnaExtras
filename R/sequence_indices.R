# ==============================================================================
# SEQUENCE INDICES TOOLKIT - COMPLEXITY SYSTEMS VERSION
# ==============================================================================
# 
# A comprehensive toolkit for computing sequence-level indices and statistics
# with focus on complexity systems, Markov chain dynamics, stability measures,
# and favorable/unfavorable state analysis.
#
# Main function: compute_sequence_indices()
#
# ==============================================================================

#' Compute Comprehensive Sequence Indices - Complexity Systems Approach
#'
#' Calculates a wide range of sequence-level statistics and indices that 
#' characterize the structure, dynamics, and complexity of individual sequences
#' with emphasis on system stability, Markov chain properties, and state favorability.
#' Supports both data.frame input and group_tna objects from the tna package.
#'
#' @param data Data frame with sequences in wide format (columns are time points) OR a group_tna object
#' @param group_col Column name or index containing group information (optional). 
#'   Ignored if data is a group_tna object.
#' @param all_possible_states Vector of all possible states (optional, will be inferred if not provided)
#' @param favorable_states Vector of states considered "favorable" for system analysis (optional)
#' @param min_length Minimum sequence length to include (default: 1)
#' @param return_group_summary Whether to return group-level summaries (default: FALSE)
#'
#' @return Data frame with sequence indices, or list with individual and group summaries
#'
#' @examples
#' \dontrun{
#' # With data.frame
#' data <- read.csv("sequence_data.csv")
#' # Define favorable states for your system
#' favorable_states <- c("consensus", "cohesion", "synthesis", "coregulate")
#' indices <- compute_sequence_indices(data, group_col = "Group", favorable_states = favorable_states)
#' 
#' # With group_tna object
#' # group_tna_obj <- tna::create_group_tna(...)
#' # indices <- compute_sequence_indices(group_tna_obj, favorable_states = favorable_states)
#' }
#'
#' @export
compute_sequence_indices <- function(data, 
                                   group_col = NULL, 
                                   all_possible_states = NULL,
                                   favorable_states = NULL,
                                   min_length = 1,
                                   return_group_summary = FALSE) {
  
  # =====================================================================
  # INPUT VALIDATION AND GROUP_TNA SUPPORT
  # =====================================================================
  
  # Check if input is a group_tna object
  if (is_group_tna(data)) {
    cat("Detected group_tna object, converting to tnaExtras format...\n")
    converted <- convert_group_tna(data)
    data <- converted$data
    group_col <- converted$group_col
    group_info <- converted$group_info
    
    cat("Successfully converted group_tna object:\n")
    cat("  Label:", group_info$label %||% "Unknown", "\n")
    cat("  Groups:", paste(group_info$levels, collapse = ", "), "\n")
    cat("  Total sequences:", nrow(data), "\n")
    
    # Use state labels from group_tna if available and not specified
    if (is.null(all_possible_states) && !is.null(group_info$state_labels)) {
      all_possible_states <- group_info$state_labels
      cat("  Using state labels from group_tna:", paste(all_possible_states, collapse = ", "), "\n")
    }
  } else {
    group_info <- NULL
  }
  
  # Input validation
  if (!is.data.frame(data)) {
    stop("data must be a data frame or group_tna object")
  }
  
  if (ncol(data) < 2) {
    stop("data must have at least 2 columns")
  }
  
  # Handle group column
  has_groups <- !is.null(group_col)
  if (has_groups) {
    if (is.character(group_col)) {
      if (!group_col %in% names(data)) {
        stop("Group column '", group_col, "' not found in data")
      }
      group_idx <- which(names(data) == group_col)
    } else if (is.numeric(group_col)) {
      group_idx <- group_col
      group_col <- names(data)[group_idx]
    } else {
      stop("group_col must be either a column name or column index")
    }
    
    groups <- data[[group_col]]
    sequence_data <- data[, -group_idx, drop = FALSE]
  } else {
    groups <- NULL
    sequence_data <- data
  }
  
  # Determine all possible states if not provided
  if (is.null(all_possible_states)) {
    all_values <- unlist(sequence_data)
    all_possible_states <- unique(all_values[!is.na(all_values) & all_values != ""])
    all_possible_states <- sort(all_possible_states)
  }
  
  n_possible_states <- length(all_possible_states)
  
  # Initialize results data frame - CORRECTED BASED ON USER FEEDBACK
  n_sequences <- nrow(sequence_data)
  results <- data.frame(
    sequence_id = 1:n_sequences,
    
    # Essential structure measures
    sequence_length = numeric(n_sequences),
    valid_observations = numeric(n_sequences),
    valid_proportion = numeric(n_sequences),  # Keep valid_proportion instead of missing measures
    
    # Core complexity measures
    unique_states_visited = numeric(n_sequences),
    mean_spell_duration = numeric(n_sequences),
    longitudinal_entropy = numeric(n_sequences),
    simpson_diversity = numeric(n_sequences),  # Keep simpson_diversity as requested
    
    # System dynamics (corrected based on feedback)
    self_loop_tendency = numeric(n_sequences),   # Keep self_loop_tendency instead of change_volatility
    transition_rate = numeric(n_sequences),      # Keep transition_rate instead of others
    transition_complexity = numeric(n_sequences),   # Keep this, it's unique
    
    # Initial state influence (most meaningful ones)
    initial_state_persistence = numeric(n_sequences),
    initial_state_influence_decay = numeric(n_sequences),
    
    # Feedback loops (most important)
    cyclic_feedback_strength = numeric(n_sequences),  # Keep this, remove others
    
    # Markov chain measures removed as requested
    
    # Temporal measures
    first_state = character(n_sequences),
    last_state = character(n_sequences),
    attractor_state = character(n_sequences),
    attractor_strength = numeric(n_sequences),  # Will fix the calculation
    
    # Emergent patterns
    emergent_state = character(n_sequences),
    emergent_state_persistence = numeric(n_sequences),
    emergent_state_proportion = numeric(n_sequences),
    
    # Advanced measures
    integrative_potential = numeric(n_sequences),
    complexity_index = numeric(n_sequences),
    
    stringsAsFactors = FALSE
  )
  
  # Add group column if present
  if (has_groups) {
    results[[group_col]] <- groups
  }
  
  # Add essential favorable state measures if favorable states are specified
  if (!is.null(favorable_states)) {
    results$proportion_favorable_states <- numeric(n_sequences)
    results$favorable_state_stability <- numeric(n_sequences)  # Keep only essential measures
  }
  
  # Compute indices for each sequence
  for (i in 1:n_sequences) {
    seq_values <- as.character(sequence_data[i, ])
    indices <- compute_single_sequence_indices(
      seq_values, 
      all_possible_states, 
      favorable_states
    )
    
    # Fill in the results
    for (col_name in names(indices)) {
      if (col_name %in% names(results)) {
        results[i, col_name] <- indices[[col_name]]
      }
    }
  }
  
  # Filter by minimum length
  if (min_length > 1) {
    valid_sequences <- results$valid_observations >= min_length
    results <- results[valid_sequences, ]
    
    if (has_groups) {
      groups <- groups[valid_sequences]
    }
  }
  
  # Return results
  if (return_group_summary && has_groups) {
    group_summary <- compute_group_summaries(results, groups)
    return(list(
      individual_indices = results,
      group_summaries = group_summary,
      parameters = list(
        n_sequences = nrow(results),
        n_possible_states = n_possible_states,
        all_possible_states = all_possible_states,
        favorable_states = favorable_states,
        min_length = min_length,
        group_tna_info = group_info  # Store original group_tna metadata
      )
    ))
  } else {
    # Add group_tna_info to simple results if from group_tna object
    if (!is.null(group_info)) {
      attr(results, "group_tna_info") <- group_info
    }
    return(results)
  }
}

#' Compute indices for a single sequence - Complexity Systems Approach
#'
#' @param seq_values Character vector representing a single sequence
#' @param all_possible_states Vector of all possible states
#' @param favorable_states Vector of favorable states for system analysis
#' @return List of computed indices
compute_single_sequence_indices <- function(seq_values, all_possible_states, favorable_states = NULL) {
  
  # Basic sequence processing - fix NA handling
  sequence_length <- length(seq_values)
  valid_states <- seq_values[!is.na(seq_values) & seq_values != "" & seq_values != "NA"]
  valid_observations <- length(valid_states)
  
  # Handle empty sequences
  if (valid_observations == 0) {
    return(create_empty_indices_list(sequence_length, valid_observations))
  }
  
  # Core analysis - CORRECTED BASED ON FEEDBACK
  spells <- compute_spells(valid_states)
  unique_states_visited <- length(spells$unique_states)
  mean_spell_duration <- mean(spells$durations)
  valid_proportion <- valid_observations / sequence_length
  
  # Diversity measures - keep both entropy and simpson as requested
  state_freqs <- table(valid_states)
  state_props <- state_freqs / valid_observations
  longitudinal_entropy <- -sum(state_props * log(state_props))
  simpson_diversity <- 1 - sum(state_props^2)
  
  # Transition analysis for complexity measures
  transitions <- compute_transitions(valid_states)
  
  # CORRECTED COMPLEXITY MEASURES BASED ON FEEDBACK
  
  # Self-loop tendency - keep this instead of change_volatility
  feedback_measures <- compute_feedback_loops(valid_states, transitions$matrix)
  self_loop_tendency <- feedback_measures$self_loop_tendency
  
  # Transition rate - keep this instead of proportion_transitions
  transition_rate <- ifelse(valid_observations <= 1, 0, 
                          transitions$count / (valid_observations - 1))
  
  # Transition Complexity - accounts for favorable/unfavorable
  transition_complexity <- compute_transition_complexity(valid_states, transitions$matrix, favorable_states)
  
  # ESSENTIAL INITIAL STATE MEASURES
  initial_measures <- compute_initial_state_influence(valid_states)
  
  # Markov measures removed as requested
  
  # Temporal measures - FIX LAST STATE ISSUE
  first_state <- valid_states[1]
  last_state <- valid_states[length(valid_states)]  # Now correctly gets last valid state, not NA
  attractor_info <- compute_attractor_state(valid_states)
  attractor_state <- attractor_info$state
  attractor_strength <- attractor_info$strength
  
  # Emergent state measures
  emergent_info <- compute_emergent_state(valid_states)
  
  # Advanced measures
  integrative_potential <- compute_integrative_potential(valid_states, favorable_states)
  complexity_index <- compute_complexity_index(valid_states, length(spells$durations), transitions$count)
  
  # Create corrected results list - BASED ON USER FEEDBACK
  indices <- list(
    sequence_length = sequence_length,
    valid_observations = valid_observations,
    valid_proportion = valid_proportion,
    
    unique_states_visited = unique_states_visited,
    mean_spell_duration = mean_spell_duration,
    longitudinal_entropy = longitudinal_entropy,
    simpson_diversity = simpson_diversity,
    
    # System dynamics (corrected based on feedback)
    self_loop_tendency = self_loop_tendency,
    transition_rate = transition_rate,
    transition_complexity = transition_complexity,
    
    # Initial state influence (essential only)
    initial_state_persistence = initial_measures$persistence,
    initial_state_influence_decay = initial_measures$influence_decay,
    
    # Feedback loops (most important)
    cyclic_feedback_strength = feedback_measures$cyclic_feedback_strength,
    
    first_state = first_state,
    last_state = last_state,
    attractor_state = attractor_state,
    attractor_strength = attractor_strength,  # Now fixed to show meaningful values
    
    emergent_state = emergent_info$state,
    emergent_state_persistence = emergent_info$persistence,
    emergent_state_proportion = emergent_info$proportion,
    
    integrative_potential = integrative_potential,
    complexity_index = complexity_index
  )
  
  # Add essential favorable state measures if favorable states specified
  if (!is.null(favorable_states)) {
    favorable_measures <- compute_favorable_state_measures(valid_states, favorable_states)
    # Only keep essential measures to avoid redundancy
    indices$proportion_favorable_states <- favorable_measures$proportion_favorable_states
    indices$favorable_state_stability <- favorable_measures$favorable_state_stability
  }
  
  return(indices)
}

# ==============================================================================
# COMPLEXITY SYSTEMS HELPER FUNCTIONS
# ==============================================================================

#' Compute system stability - how stable the system is overall
compute_system_stability <- function(states) {
  if (length(states) <= 1) return(1)
  
  # Stability = 1 - (transitions / max_possible_transitions)
  # Higher values = more stable system
  changes <- sum(states[-1] != states[-length(states)])
  max_changes <- length(states) - 1
  stability <- 1 - (changes / max_changes)
  
  return(stability)
}

#' Compute change volatility - how much the system changes
compute_change_volatility <- function(states) {
  if (length(states) <= 1) return(0)
  
  # Volatility = proportion of time points where state changes
  # Higher values = more volatile system
  changes <- sum(states[-1] != states[-length(states)])
  volatility <- changes / (length(states) - 1)
  
  return(volatility)
}

#' Compute transition complexity considering favorable/unfavorable states
compute_transition_complexity <- function(states, transition_matrix, favorable_states = NULL) {
  if (length(states) <= 1) return(0)
  
  # Basic complexity: diversity of transitions used
  n_unique_states <- length(unique(states))
  max_possible_transitions <- n_unique_states * (n_unique_states - 1)
  
  if (max_possible_transitions == 0) return(0)
  
  observed_transitions <- sum(transition_matrix > 0)
  basic_complexity <- observed_transitions / max_possible_transitions
  
  # If favorable states specified, weight by favorability
  if (!is.null(favorable_states)) {
    # Count favorable-to-unfavorable and unfavorable-to-favorable transitions
    # These indicate system instability
    favorable_states_in_seq <- intersect(rownames(transition_matrix), favorable_states)
    unfavorable_states_in_seq <- setdiff(rownames(transition_matrix), favorable_states)
    
    if (length(favorable_states_in_seq) > 0 && length(unfavorable_states_in_seq) > 0) {
      # Transitions from favorable to unfavorable (bad)
      fav_to_unfav <- sum(transition_matrix[favorable_states_in_seq, unfavorable_states_in_seq, drop = FALSE])
      # Transitions from unfavorable to favorable (good)
      unfav_to_fav <- sum(transition_matrix[unfavorable_states_in_seq, favorable_states_in_seq, drop = FALSE])
      
      total_transitions <- sum(transition_matrix)
      if (total_transitions > 0) {
        instability_factor <- (fav_to_unfav - unfav_to_fav) / total_transitions
        # Adjust complexity by instability (more instability = higher complexity)
        basic_complexity <- basic_complexity * (1 + abs(instability_factor))
      }
    }
  }
  
  return(basic_complexity)
}

#' Compute initial state influence measures
compute_initial_state_influence <- function(states) {
  if (length(states) <= 1) {
    return(list(persistence = 1, return_rate = 0, influence_decay = 0))
  }
  
  initial_state <- states[1]
  
  # 1. Initial state persistence - how long does initial state last?
  initial_spell_length <- 1
  for (i in 2:length(states)) {
    if (states[i] == initial_state) {
      initial_spell_length <- initial_spell_length + 1
    } else {
      break
    }
  }
  persistence <- initial_spell_length / length(states)
  
  # 2. Return rate to initial state
  initial_positions <- which(states == initial_state)
  if (length(initial_positions) <= 1) {
    return_rate <- 0
  } else {
    # How often does the system return to initial state after leaving it?
    returns <- length(initial_positions) - 1  # Subtract the initial occurrence
    opportunities <- length(states) - initial_spell_length
    return_rate <- ifelse(opportunities == 0, 0, returns / opportunities)
  }
  
  # 3. Influence decay - how does initial state influence decay over time?
  # Measure presence of initial state in different thirds of sequence
  n <- length(states)
  first_third <- states[1:ceiling(n/3)]
  last_third <- states[ceiling(2*n/3):n]
  
  initial_prop_early <- sum(first_third == initial_state) / length(first_third)
  initial_prop_late <- sum(last_third == initial_state) / length(last_third)
  
  influence_decay <- initial_prop_early - initial_prop_late
  
  return(list(
    persistence = persistence,
    return_rate = return_rate,
    influence_decay = influence_decay
  ))
}

#' Compute feedback loop measures
compute_feedback_loops <- function(states, transition_matrix) {
  if (length(states) <= 1) {
    return(list(self_loop_tendency = 0, return_to_previous_rate = 0, cyclic_feedback_strength = 0))
  }
  
  # 1. Self-loop tendency - staying in same state
  self_loops <- 0
  total_transitions <- 0
  for (i in 1:(length(states) - 1)) {
    if (!is.na(states[i]) && !is.na(states[i + 1])) {
      if (states[i] == states[i + 1]) {
        self_loops <- self_loops + 1
      }
      total_transitions <- total_transitions + 1
    }
  }
  self_loop_tendency <- ifelse(total_transitions == 0, 0, self_loops / total_transitions)
  
  # 2. Return to previous state rate (A -> B -> A pattern)
  return_to_previous <- 0
  opportunities <- 0
  if (length(states) >= 3) {
    for (i in 3:length(states)) {
      if (!is.na(states[i-2]) && !is.na(states[i-1]) && !is.na(states[i])) {
        if (states[i-2] != states[i-1]) {  # There was a transition
          opportunities <- opportunities + 1
          if (states[i] == states[i-2]) {  # Returned to previous state
            return_to_previous <- return_to_previous + 1
          }
        }
      }
    }
  }
  return_to_previous_rate <- ifelse(opportunities == 0, 0, return_to_previous / opportunities)
  
  # 3. Cyclic feedback strength - longer cycles
  cyclic_feedback_strength <- compute_cyclic_strength(states)
  
  return(list(
    self_loop_tendency = self_loop_tendency,
    return_to_previous_rate = return_to_previous_rate,
    cyclic_feedback_strength = cyclic_feedback_strength
  ))
}

#' Compute Markov chain parameters
compute_markov_parameters <- function(states, transition_matrix) {
  if (length(states) <= 2) {
    return(list(stationarity = 0, mixing_time = 0, entropy_rate = 0))
  }
  
  # Convert counts to probabilities
  transition_probs <- transition_matrix
  for (i in 1:nrow(transition_matrix)) {
    row_sum <- sum(transition_matrix[i, ])
    if (row_sum > 0) {
      transition_probs[i, ] <- transition_matrix[i, ] / row_sum
    }
  }
  
  # 1. Stationarity measure - how close to steady state
  # Use the variation in state frequencies over time
  n <- length(states)
  if (n >= 6) {
    first_half <- states[1:floor(n/2)]
    second_half <- states[(floor(n/2)+1):n]
    
    freq1 <- table(factor(first_half, levels = unique(states)))
    freq2 <- table(factor(second_half, levels = unique(states)))
    
    prop1 <- freq1 / sum(freq1)
    prop2 <- freq2 / sum(freq2)
    
    # Stationarity = 1 - total variation distance
    stationarity <- 1 - 0.5 * sum(abs(prop1 - prop2))
  } else {
    stationarity <- 0
  }
  
  # 2. Mixing time approximation - simpler approach
  # Based on how quickly state distributions stabilize
  mixing_time <- 0
  
  if (length(states) >= 10) {
    # Compare early vs late state distributions
    early_third <- states[1:floor(length(states)/3)]
    late_third <- states[ceiling(2*length(states)/3):length(states)]
    
    early_freq <- table(factor(early_third, levels = unique(states)))
    late_freq <- table(factor(late_third, levels = unique(states)))
    
    early_prop <- early_freq / sum(early_freq)
    late_prop <- late_freq / sum(late_freq)
    
    # Mixing time is inversely related to how similar early and late distributions are
    similarity <- 1 - 0.5 * sum(abs(early_prop - late_prop))
    mixing_time <- ifelse(similarity > 0.1, 1/similarity, 10)  # Capped at 10
  }
  
  # 3. Entropy rate - average entropy per step
  entropy_rate <- 0
  if (nrow(transition_probs) > 0) {
    for (i in 1:nrow(transition_probs)) {
      row_probs <- transition_probs[i, ]
      row_probs <- row_probs[row_probs > 0]  # Remove zeros
      if (length(row_probs) > 0) {
        row_entropy <- -sum(row_probs * log(row_probs))
        # Weight by stationary probability (approximate with observed frequency)
        state_freq <- sum(states == rownames(transition_probs)[i]) / length(states)
        entropy_rate <- entropy_rate + state_freq * row_entropy
      }
    }
  }
  
  return(list(
    stationarity = stationarity,
    mixing_time = mixing_time,
    entropy_rate = entropy_rate
  ))
}

#' Compute favorable state measures
compute_favorable_state_measures <- function(states, favorable_states) {
  # Convert to favorable/unfavorable
  is_favorable <- states %in% favorable_states
  
  # Proportion of favorable states
  proportion_favorable <- mean(is_favorable)
  
  # Run analysis for favorable and unfavorable runs
  favorable_runs <- compute_runs(as.numeric(is_favorable), 1)
  unfavorable_runs <- compute_runs(as.numeric(is_favorable), 0)
  
  # Stability measures
  favorable_state_stability <- compute_state_type_stability(states, favorable_states, TRUE)
  unfavorable_state_volatility <- compute_state_type_stability(states, favorable_states, FALSE)
  
  # Ratio measure
  favorable_to_unfavorable_ratio <- ifelse(unfavorable_runs$count == 0, 
                                          Inf, 
                                          favorable_runs$count / unfavorable_runs$count)
  
  return(list(
    proportion_favorable_states = proportion_favorable,
    favorable_state_runs = favorable_runs$count,
    unfavorable_state_runs = unfavorable_runs$count,
    longest_favorable_run = favorable_runs$longest,
    longest_unfavorable_run = unfavorable_runs$longest,
    favorable_state_stability = favorable_state_stability,
    unfavorable_state_volatility = unfavorable_state_volatility,
    favorable_to_unfavorable_ratio = favorable_to_unfavorable_ratio
  ))
}

#' Compute stability for favorable or unfavorable states
compute_state_type_stability <- function(states, favorable_states, is_favorable_type) {
  if (is_favorable_type) {
    target_states <- favorable_states
  } else {
    target_states <- setdiff(unique(states), favorable_states)
  }
  
  # Find spells of target state type
  is_target <- states %in% target_states
  if (sum(is_target) == 0) return(0)
  
  # Compute average spell length for target state type
  rle_result <- rle(is_target)
  target_spells <- rle_result$lengths[rle_result$values == TRUE]
  
  if (length(target_spells) == 0) return(0)
  
  stability <- mean(target_spells)
  return(stability)
}

# ==============================================================================
# EXISTING HELPER FUNCTIONS (keeping the ones that work well)
# ==============================================================================

#' Compute spells in a sequence
compute_spells <- function(states) {
  if (length(states) == 0) {
    return(list(states = character(0), durations = numeric(0), unique_states = character(0)))
  }
  
  # Run-length encoding
  rle_result <- rle(states)
  
  return(list(
    states = rle_result$values,
    durations = rle_result$lengths,
    unique_states = unique(states)
  ))
}

#' Compute transitions in a sequence
compute_transitions <- function(states) {
  if (length(states) <= 1) {
    return(list(count = 0, matrix = matrix()))
  }
  
  # Count transitions
  transition_count <- sum(states[-1] != states[-length(states)])
  
  # Create transition matrix
  unique_states <- unique(states)
  n_states <- length(unique_states)
  trans_matrix <- matrix(0, nrow = n_states, ncol = n_states)
  rownames(trans_matrix) <- unique_states
  colnames(trans_matrix) <- unique_states
  
  for (i in 1:(length(states) - 1)) {
    from_state <- states[i]
    to_state <- states[i + 1]
    trans_matrix[from_state, to_state] <- trans_matrix[from_state, to_state] + 1
  }
  
  return(list(count = transition_count, matrix = trans_matrix))
}

#' Compute comprehensive durations considering all possible states
compute_comprehensive_durations <- function(observed_states, observed_durations, all_possible_states) {
  # Create duration vector for all possible states
  comprehensive_durations <- numeric(length(all_possible_states))
  names(comprehensive_durations) <- all_possible_states
  
  # Fill in observed durations
  for (i in seq_along(observed_states)) {
    state <- observed_states[i]
    if (state %in% all_possible_states) {
      comprehensive_durations[state] <- comprehensive_durations[state] + observed_durations[i]
    }
  }
  
  return(comprehensive_durations)
}

#' Compute concentration index (similar to Gini coefficient)
compute_concentration_index <- function(frequencies) {
  if (length(frequencies) <= 1) return(0)
  
  frequencies <- sort(frequencies)
  n <- length(frequencies)
  index <- 2 * sum((1:n) * frequencies) / (n * sum(frequencies)) - (n + 1) / n
  
  return(index)
}

#' Compute subsequences
compute_subsequences <- function(states) {
  if (length(states) <= 1) {
    return(list(unique_subsequences = states, longest_repeat = 0))
  }
  
  # Generate all subsequences of length 2+
  subsequences <- character(0)
  for (len in 2:min(length(states), 10)) {  # Limit to avoid explosion
    for (start in 1:(length(states) - len + 1)) {
      subseq <- paste(states[start:(start + len - 1)], collapse = "-")
      subsequences <- c(subsequences, subseq)
    }
  }
  
  unique_subsequences <- unique(subsequences)
  
  # Find longest repeated subsequence
  longest_repeat <- 0
  for (subseq in unique_subsequences) {
    occurrences <- sum(grepl(subseq, subsequences, fixed = TRUE))
    if (occurrences > 1) {
      subseq_length <- length(strsplit(subseq, "-")[[1]])
      longest_repeat <- max(longest_repeat, subseq_length)
    }
  }
  
  return(list(unique_subsequences = unique_subsequences, longest_repeat = longest_repeat))
}

#' Compute cyclic pattern strength
compute_cyclic_strength <- function(states) {
  if (length(states) <= 2) return(0)
  
  # Look for repeating patterns
  max_strength <- 0
  max_period <- min(length(states) %/% 2, 10)  # Limit search
  
  for (period in 2:max_period) {
    matches <- 0
    comparisons <- 0
    
    for (i in 1:(length(states) - period)) {
      if (i + period <= length(states)) {
        if (states[i] == states[i + period]) {
          matches <- matches + 1
        }
        comparisons <- comparisons + 1
      }
    }
    
    if (comparisons > 0) {
      strength <- matches / comparisons
      max_strength <- max(max_strength, strength)
    }
  }
  
  return(max_strength)
}

#' Compute attractor state - the state the system is drawn to
#' @param states Vector of states
#' @param method Method for combining components: "weighted", "all", or custom vector c(freq_weight, persistence_weight, return_weight)
#' @param normalize Whether to normalize the final strength (default TRUE)
compute_attractor_state <- function(states, method = "weighted", normalize = TRUE) {
  if (length(states) <= 1) {
    return(list(state = states[1], strength = 1))
  }
  
  # Determine weights based on method
  if (is.character(method)) {
    if (method == "weighted") {
      weights <- c(freq = 1.0, persistence = 0.5, return = 0.3)
    } else if (method == "all") {
      weights <- c(freq = 1.0, persistence = 1.0, return = 1.0)
    } else {
      stop("Method must be 'weighted', 'all', or a numeric vector c(freq_weight, persistence_weight, return_weight)")
    }
  } else if (is.numeric(method) && length(method) == 3) {
    weights <- c(freq = method[1], persistence = method[2], return = method[3])
  } else {
    stop("Method must be 'weighted', 'all', or a numeric vector c(freq_weight, persistence_weight, return_weight)")
  }
  
  # Attractor state is determined by:
  # 1. Overall frequency (like modal state)
  # 2. Tendency to persist when reached
  # 3. Tendency to return to after leaving
  
  state_counts <- table(states)
  state_props <- state_counts / length(states)
  
  # For each state, compute its "attractor strength"
  attractor_strengths <- numeric(length(state_counts))
  names(attractor_strengths) <- names(state_counts)
  
  for (state in names(state_counts)) {
    # Base strength from frequency
    freq_strength <- state_props[state]
    
    # Persistence strength - how long it lasts when reached
    state_positions <- which(states == state)
    if (length(state_positions) > 0) {
      spells <- compute_state_spells(states, state)
      persistence_strength <- mean(spells) / length(states)
    } else {
      persistence_strength <- 0
    }
    
    # Return strength - how often system returns to this state
    if (length(state_positions) > 1) {
      gaps <- diff(state_positions) - 1  # Gaps between occurrences
      return_strength <- 1 / (1 + mean(gaps))  # Inverse of average gap
    } else {
      return_strength <- 0
    }
    
    # Combined attractor strength using specified weights
    attractor_strengths[state] <- weights["freq"] * freq_strength + 
                                 weights["persistence"] * persistence_strength + 
                                 weights["return"] * return_strength
  }
  
  # Attractor state is simply the most frequent state (modal state)
  attractor_state <- names(which.max(state_props))
  
  # Attractor strength is simply the proportion of that state
  attractor_strength <- as.numeric(state_props[attractor_state])
  
  return(list(state = attractor_state, strength = attractor_strength))
}

#' Compute emergent state - only if it becomes MORE dominant than the attractor
compute_emergent_state <- function(states) {
  if (length(states) <= 3) {
    return(list(state = NA_character_, persistence = 0, proportion = 0))
  }
  
  # First identify the attractor state (most frequent)
  state_props <- table(states) / length(states)
  attractor_state <- names(which.max(state_props))
  attractor_proportion <- max(state_props)
  
  # Get attractor state's maximum persistence
  attractor_spells <- compute_state_spells(states, attractor_state)
  attractor_max_persistence <- ifelse(length(attractor_spells) > 0, max(attractor_spells), 0)
  
  initial_state <- states[1]
  
  # Calculate initial state first spell (not max, just first spell)
  initial_first_spell <- 1
  for (i in 2:length(states)) {
    if (states[i] != initial_state) {
      break
    }
    initial_first_spell <- initial_first_spell + 1
  }
  
  # Find states that truly emerge beyond the attractor's dominance
  best_emergent <- list(state = NA_character_, persistence = 0, proportion = 0)
  
  # Look for any state that becomes MORE dominant than the attractor
  for (state in unique(states)) {
    # Find all spells of this state
    state_spells <- compute_state_spells(states, state)
    
    if (length(state_spells) > 0) {
      # Longest persistence of this state
      max_persistence <- max(state_spells)
      state_proportion <- sum(states == state) / length(states)
      
      # Consider as emergent ONLY if it truly emerges beyond the initial dominant pattern:
      # 1. It has significant persistence (at least 3 time points)
      # 2. AND either:
      #    a) It's different from attractor AND surpasses attractor in persistence
      #    b) It IS the attractor but wasn't initially dominant, showing true emergence
      
      is_emergent <- FALSE
      
      if (state != attractor_state && max_persistence >= 3 && 
          max_persistence > attractor_max_persistence) {
        # Different from attractor and MORE persistent - true emergence of new dominant state
        is_emergent <- TRUE
      } else if (state == attractor_state && state != initial_state && max_persistence >= 3) {
        # This state became the attractor despite not being initial - true emergence
        is_emergent <- TRUE
      } else if (state == attractor_state && state == initial_state && 
                 max_persistence > initial_first_spell * 2 && max_persistence >= 3) {
        # Initial state that becomes attractor through MUCH higher persistence (2x initial)
        is_emergent <- TRUE
      }
      
      if (is_emergent) {
        # Update best emergent state if this one is better
        if (max_persistence > best_emergent$persistence) {
          best_emergent <- list(
            state = state,
            persistence = max_persistence,
            proportion = state_proportion
          )
        }
      }
    }
  }
  
  return(best_emergent)
}

#' Compute spells for a specific state
compute_state_spells <- function(states, target_state) {
  if (length(states) == 0 || !target_state %in% states) {
    return(numeric(0))
  }
  
  # Use run-length encoding to find spells
  rle_result <- rle(states == target_state)
  target_spells <- rle_result$lengths[rle_result$values == TRUE]
  
  return(target_spells)
}

#' Compute integrative potential using the mathematical formula
#' I_integr(x) = (Σ is.pos(xi) * i^ω) / (Σ i^ω)
#' @param states Vector of states
#' @param favorable_states Vector of states considered favorable/positive
#' @param omega Power parameter for position weighting (default 1)
#' @param with_missing Boolean to include missing values (default FALSE)
compute_integrative_potential <- function(states, favorable_states = NULL, omega = 1, with_missing = FALSE) {
  if (length(states) <= 1) return(1)
  
  # Validate input
  if (!is.vector(states)) {
    stop("[!] states must be a vector of states.")
  }
  
  # Handle missing values
  if (with_missing) {
    states <- na.omit(states)
  }
  
  # Create binary indicator for positive/favorable states
  is_positive <- states %in% favorable_states
  
  # Position weights: i^ω for i = 1, 2, ..., n
  n <- length(states)
  position_weights <- (1:n)^omega
  sum_integ <- sum(position_weights, na.rm = TRUE)
  
  # Calculate the numerator: Σ is.pos(xi) * i^ω
  numerator <- sum(is_positive * position_weights, na.rm = TRUE)
  
  # Calculate the integrative potential
  integrative_potential <- numerator / sum_integ
  
  return(integrative_potential)
}

#' Compute complexity index - clear measure of sequence unpredictability and variability
#' Higher values = more complex/unpredictable sequences
#' Range: 0 to 1, where 0 = completely predictable, 1 = maximally complex
compute_complexity_index <- function(states, spell_count, transition_count) {
  if (length(states) <= 1) return(0)
  
  sequence_length <- length(states)
  n_unique <- length(unique(states))
  
  # 1. ENTROPY COMPONENT (40%): How evenly distributed are the states?
  # 0 = all same state, 1 = perfectly even distribution
  state_probs <- table(states) / length(states)
  entropy <- -sum(state_probs * log(state_probs))
  max_entropy <- log(n_unique)
  entropy_component <- ifelse(max_entropy > 0, entropy / max_entropy, 0)
  
  # 2. TRANSITION COMPONENT (40%): How often does the sequence change?
  # 0 = no transitions, 1 = changes every time point
  max_transitions <- sequence_length - 1
  transition_component <- ifelse(max_transitions > 0, transition_count / max_transitions, 0)
  
  # 3. VARIABILITY COMPONENT (20%): How variable are the spell durations?
  # 0 = all spells same length, higher = more variable spell lengths
  spells <- compute_spells(states)
  if (length(spells$durations) > 1) {
    cv_spells <- sd(spells$durations) / mean(spells$durations)
    variability_component <- min(cv_spells, 1)  # Cap at 1
  } else {
    variability_component <- 0
  }
  
  # FINAL COMPLEXITY INDEX: Weighted combination
  # 0.4 * entropy + 0.4 * transitions + 0.2 * variability
  complexity <- 0.4 * entropy_component + 
                0.4 * transition_component + 
                0.2 * variability_component
  
  return(complexity)
}

#' Compute recurrence rate - how often states return after absence
compute_recurrence_rate <- function(states) {
  if (length(states) <= 2) return(0)
  
  # Recurrence rate measures how often states return after being absent
  # Higher values indicate more recurring, cyclical patterns
  
  unique_states <- unique(states)
  if (length(unique_states) <= 1) return(0)
  
  total_recurrences <- 0
  total_possible_recurrences <- 0
  
  for (state in unique_states) {
    state_positions <- which(states == state)
    
    if (length(state_positions) >= 2) {
      # For each occurrence after the first, check if it's a recurrence
      for (i in 2:length(state_positions)) {
        current_pos <- state_positions[i]
        previous_pos <- state_positions[i-1]
        
        # It's a recurrence if there's at least one different state in between
        if (current_pos - previous_pos > 1) {
          # Check if there are actually different states in between
          between_states <- states[(previous_pos + 1):(current_pos - 1)]
          if (any(between_states != state)) {
            total_recurrences <- total_recurrences + 1
          }
        }
        
        # Each occurrence after the first is a potential recurrence
        total_possible_recurrences <- total_possible_recurrences + 1
      }
    }
  }
  
  # Calculate recurrence rate
  recurrence_rate <- ifelse(total_possible_recurrences > 0, 
                           total_recurrences / total_possible_recurrences, 0)
  
  return(recurrence_rate)
}

#' Compute runs of a specific value
compute_runs <- function(binary_seq, target_value) {
  if (length(binary_seq) == 0) {
    return(list(count = 0, longest = 0))
  }
  
  # Find runs of target value
  runs <- rle(binary_seq)
  target_runs <- runs$lengths[runs$values == target_value]
  
  if (length(target_runs) == 0) {
    return(list(count = 0, longest = 0))
  }
  
  return(list(
    count = length(target_runs),
    longest = max(target_runs)
  ))
}

#' Create empty indices list for sequences with no valid data - CORRECTED
create_empty_indices_list <- function(sequence_length, valid_observations) {
  return(list(
    sequence_length = sequence_length,
    valid_observations = valid_observations,
    valid_proportion = valid_observations / sequence_length,
    
    unique_states_visited = 0,
    mean_spell_duration = 0,
    longitudinal_entropy = 0,
    simpson_diversity = 0,
    
    self_loop_tendency = 0,
    transition_rate = 0,
    transition_complexity = 0,
    
    initial_state_persistence = 0,
    initial_state_influence_decay = 0,
    
    cyclic_feedback_strength = 0,
    
    first_state = NA_character_,
    last_state = NA_character_,
    attractor_state = NA_character_,
    attractor_strength = 0,
    
    emergent_state = NA_character_,
    emergent_state_persistence = 0,
    emergent_state_proportion = 0,
    
    integrative_potential = 0,
    complexity_index = 0
  ))
}

#' Compute group-level summaries
compute_group_summaries <- function(results, groups) {
  if (is.null(groups)) return(NULL)
  
  # Numeric columns for summarization
  numeric_cols <- sapply(results, is.numeric)
  numeric_data <- results[, numeric_cols, drop = FALSE]
  
  # Compute summaries by group
  group_summary <- aggregate(numeric_data, by = list(Group = groups), 
                           FUN = function(x) c(mean = mean(x, na.rm = TRUE),
                                              sd = sd(x, na.rm = TRUE),
                                              median = median(x, na.rm = TRUE),
                                              min = min(x, na.rm = TRUE),
                                              max = max(x, na.rm = TRUE)))
  
  return(group_summary)
}

# ==============================================================================
# SUMMARY AND VISUALIZATION FUNCTIONS
# ==============================================================================

#' Print summary of sequence indices - Complexity Systems Version
#'
#' @param indices_result Result from compute_sequence_indices()
#' @param top_n Number of top sequences to show for each measure
#' @export
print_indices_summary <- function(indices_result, top_n = 5) {
  if (is.list(indices_result) && "individual_indices" %in% names(indices_result)) {
    data <- indices_result$individual_indices
    params <- indices_result$parameters
  } else {
    data <- indices_result
    params <- NULL
  }
  
  cat("=== SEQUENCE INDICES SUMMARY - COMPLEXITY SYSTEMS VERSION ===\n\n")
  
  if (!is.null(params)) {
    cat("Parameters:\n")
    cat("  Number of sequences:", params$n_sequences, "\n")
    cat("  Number of possible states:", params$n_possible_states, "\n")
    cat("  Possible states:", paste(params$all_possible_states, collapse = ", "), "\n")
    if (!is.null(params$favorable_states)) {
      cat("  Favorable states:", paste(params$favorable_states, collapse = ", "), "\n")
    }
    cat("  Minimum length:", params$min_length, "\n\n")
  }
  
  cat("Overall Statistics (STREAMLINED):\n")
  cat("  Mean sequence length:", round(mean(data$sequence_length, na.rm = TRUE), 2), "\n")
  cat("  Mean valid observations:", round(mean(data$valid_observations, na.rm = TRUE), 2), "\n")
  cat("  Mean unique states:", round(mean(data$unique_states_visited, na.rm = TRUE), 2), "\n")
  cat("  Mean spell duration:", round(mean(data$mean_spell_duration, na.rm = TRUE), 2), "\n")
  cat("  Mean entropy:", round(mean(data$longitudinal_entropy, na.rm = TRUE), 3), "\n")
  
  cat("\nComplexity Systems Measures:\n")
  cat("  Mean change volatility:", round(mean(data$change_volatility, na.rm = TRUE), 3), "\n")
  cat("  Mean transition complexity:", round(mean(data$transition_complexity, na.rm = TRUE), 3), "\n")
  
  cat("\nInitial State Influence:\n")
  cat("  Mean initial state persistence:", round(mean(data$initial_state_persistence, na.rm = TRUE), 3), "\n")
  cat("  Mean initial state influence decay:", round(mean(data$initial_state_influence_decay, na.rm = TRUE), 3), "\n")
  
  cat("\nFeedback Loops:\n")
  cat("  Mean cyclic feedback strength:", round(mean(data$cyclic_feedback_strength, na.rm = TRUE), 3), "\n")
  
  cat("\nMarkov Chain Properties:\n")
  cat("  Mean stationarity:", round(mean(data$markov_stationarity, na.rm = TRUE), 3), "\n")
  cat("  Mean entropy rate:", round(mean(data$markov_entropy_rate, na.rm = TRUE), 3), "\n")
  
  # Show favorable state measures if available
  if ("proportion_favorable_states" %in% names(data)) {
    cat("\nFavorable State Analysis (STREAMLINED):\n")
    cat("  Mean proportion favorable:", round(mean(data$proportion_favorable_states, na.rm = TRUE), 3), "\n")
    cat("  Mean favorable state stability:", round(mean(data$favorable_state_stability, na.rm = TRUE), 3), "\n")
  }
  
  # Show top sequences for key measures (STREAMLINED)
  key_measures <- c("change_volatility", "initial_state_persistence", 
                   "markov_stationarity", "integrative_potential", "complexity_index")
  
  for (measure in key_measures) {
    if (measure %in% names(data)) {
      cat("\nTop", top_n, "sequences by", measure, ":\n")
      top_idx <- order(data[[measure]], decreasing = TRUE)[1:min(top_n, nrow(data))]
      top_data <- data[top_idx, c("sequence_id", measure)]
      print(top_data)
    }
  }
}

cat("Sequence indices toolkit - Complexity Systems Version loaded successfully.\n")
cat("Main function: compute_sequence_indices()\n") 
cat("Use print_indices_summary() to view results.\n")
cat("Key improvements: favorable_states, initial state influence, feedback loops, Markov parameters\n") 
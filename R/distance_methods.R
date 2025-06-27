# ==============================================================================
# DISTANCE COMPUTATION METHODS FOR SEQUENCE CLUSTERING
# ==============================================================================
# 
# This module provides various distance computation methods for categorical
# sequence data, including traditional methods and integration with stringdist.
# All methods are optimized for performance and handle missing values robustly.
#
# ==============================================================================

#' Compute Euclidean Distance for Sequences
#'
#' Converts categorical sequences to numeric encoding and computes Euclidean distance.
#' This is the fastest method and works well for exploratory analysis.
#'
#' @param data A data frame or matrix where rows are sequences and columns are time points
#' @return A distance object of class 'dist'
#' @keywords internal
compute_euclidean_distance <- function(data) {
  if (!is.data.frame(data) && !is.matrix(data)) {
    stop("Data must be a data frame or matrix")
  }
  
  # Convert to character matrix for consistent processing
  char_matrix <- as.matrix(data)
  mode(char_matrix) <- "character"
  
  # Get all unique states excluding NA
  all_states <- unique(as.vector(char_matrix))
  all_states <- all_states[!is.na(all_states)]
  
  if (length(all_states) == 0) {
    stop("No non-missing states found in data")
  }
  
  # Create mapping using vectorized operations
  state_mapping <- seq_along(all_states)
  names(state_mapping) <- all_states
  
  # Vectorized conversion to numeric
  numeric_matrix <- matrix(state_mapping[char_matrix], nrow = nrow(data), ncol = ncol(data))
  numeric_matrix[is.na(char_matrix)] <- 0  # Handle NAs
  
  stats::dist(numeric_matrix, method = "euclidean")
}

#' Compute Hamming Distance for Sequences
#'
#' Computes position-wise differences between sequences.
#' Counts the number of positions where sequences differ.
#'
#' @param data A data frame or matrix where rows are sequences and columns are time points
#' @return A distance object of class 'dist'
#' @keywords internal
compute_hamming_distance <- function(data) {
  if (!is.data.frame(data) && !is.matrix(data)) {
    stop("Data must be a data frame or matrix")
  }
  
  # Convert to character matrix
  char_matrix <- as.matrix(data)
  mode(char_matrix) <- "character"
  char_matrix[is.na(char_matrix)] <- "NA_VALUE"  # Handle NAs explicitly
  
  n_seq <- nrow(char_matrix)
  
  # Vectorized computation using outer function
  hamming_distances <- outer(seq_len(n_seq), seq_len(n_seq), Vectorize(function(i, j) {
    if (i >= j) return(0)
    sum(char_matrix[i, ] != char_matrix[j, ])
  }))
  
  # Make symmetric
  hamming_distances <- hamming_distances + t(hamming_distances)
  
  as.dist(hamming_distances)
}

#' Compute LCS-based Distance for Sequences
#'
#' Computes dissimilarity based on Longest Common Subsequence.
#' Measures similarity by the proportion of common elements.
#'
#' @param data A data frame or matrix where rows are sequences and columns are time points
#' @return A distance object of class 'dist'
#' @keywords internal
compute_lcs_distance <- function(data) {
  if (!is.data.frame(data) && !is.matrix(data)) {
    stop("Data must be a data frame or matrix")
  }
  
  # Convert to character matrix and clean
  char_matrix <- as.matrix(data)
  mode(char_matrix) <- "character"
  
  # Pre-process sequences: convert to lists and remove NAs
  sequences <- lapply(seq_len(nrow(char_matrix)), function(i) {
    row_data <- char_matrix[i, ]
    row_data[!is.na(row_data)]
  })
  
  # Get sequence lengths
  seq_lengths <- lengths(sequences)
  
  n_seq <- length(sequences)
  dist_matrix <- matrix(0, nrow = n_seq, ncol = n_seq)
  
  # Compute pairwise distances
  for (i in seq_len(n_seq)) {
    for (j in i:n_seq) {
      if (i == j) {
        dist_matrix[i, j] <- 0
      } else {
        len1 <- seq_lengths[i]
        len2 <- seq_lengths[j]
        
        if (len1 == 0 || len2 == 0) {
          similarity <- 0
        } else {
          # Fast set intersection
          common_count <- length(intersect(sequences[[i]], sequences[[j]]))
          similarity <- common_count / max(len1, len2)
        }
        
        dist_matrix[i, j] <- dist_matrix[j, i] <- 1 - similarity
      }
    }
  }
  
  as.dist(dist_matrix)
}

#' Compute Start Position Distance for Sequences
#'
#' Computes distance based on first occurrence positions of each state.
#' Useful for analyzing timing patterns in sequences.
#'
#' @param data A data frame or matrix where rows are sequences and columns are time points
#' @return A distance object of class 'dist'
#' @keywords internal
compute_start_position_distance <- function(data) {
  if (!is.data.frame(data) && !is.matrix(data)) {
    stop("Data must be a data frame or matrix")
  }
  
  # Convert to character matrix
  char_matrix <- as.matrix(data)
  mode(char_matrix) <- "character"
  
  all_states <- unique(as.vector(char_matrix))
  all_states <- all_states[!is.na(all_states)]
  
  if (length(all_states) == 0) {
    stop("No non-missing states found in data")
  }
  
  n_seq <- nrow(char_matrix)
  n_states <- length(all_states)
  position_matrix <- matrix(ncol(data) + 1, nrow = n_seq, ncol = n_states)
  colnames(position_matrix) <- all_states
  
  # Vectorized first position finding
  for (state in all_states) {
    first_positions <- apply(char_matrix == state, 1, function(x) {
      pos <- which(x)[1]
      if (is.na(pos)) ncol(data) + 1 else pos
    })
    position_matrix[, state] <- first_positions
  }
  
  stats::dist(position_matrix, method = "euclidean")
}

#' Compute Advanced Transition Distance for Sequences
#'
#' Computes distance based on comprehensive Markov transition analysis including
#' higher-order patterns, time-weighted transitions, and pattern complexity.
#'
#' @param data A data frame or matrix where rows are sequences and columns are time points
#' @param first_order_weight Weight for basic state transitions (default: 0.3)
#' @param second_order_weight Weight for higher-order patterns (default: 0.2)
#' @param time_weighted_weight Weight for temporal emphasis (default: 0.2)
#' @param complexity_weight Weight for pattern diversity (default: 0.15)
#' @param persistence_weight Weight for state stability (default: 0.15)
#' @return A distance object of class 'dist'
#' @keywords internal
compute_transition_distance <- function(data, 
                                       first_order_weight = 0.3,
                                       second_order_weight = 0.2,
                                       time_weighted_weight = 0.2,
                                       complexity_weight = 0.15,
                                       persistence_weight = 0.15) {
  
  if (!is.data.frame(data) && !is.matrix(data)) {
    stop("Data must be a data frame or matrix")
  }
  
  # Normalize weights
  weights <- c(first_order_weight, second_order_weight, time_weighted_weight,
               complexity_weight, persistence_weight)
  weights <- weights / sum(weights)
  
  # Convert to character matrix
  char_matrix <- as.matrix(data)
  mode(char_matrix) <- "character"
  
  all_states <- unique(as.vector(char_matrix))
  all_states <- all_states[!is.na(all_states)]
  n_states <- length(all_states)
  n_seq <- nrow(char_matrix)
  
  if (n_states < 2) {
    stop("Need at least 2 different states for transition analysis")
  }
  
  # Create comprehensive feature matrices
  n_transitions <- n_states * n_states
  n_bigrams <- n_states * n_states
  
  # Feature matrices
  transition_features <- matrix(0, nrow = n_seq, ncol = n_transitions)
  second_order_features <- matrix(0, nrow = n_seq, ncol = n_bigrams)
  time_weighted_features <- matrix(0, nrow = n_seq, ncol = n_transitions)
  complexity_features <- matrix(0, nrow = n_seq, ncol = 3)  # entropy, unique transitions, repetitions
  persistence_features <- matrix(0, nrow = n_seq, ncol = n_states)
  
  # Process each sequence
  for (i in seq_len(n_seq)) {
    seq_data <- char_matrix[i, ]
    valid_positions <- which(!is.na(seq_data))
    
    if (length(valid_positions) < 2) next
    
    valid_seq <- seq_data[valid_positions]
    
    # 1. First-order transitions
    for (t in seq_len(length(valid_seq) - 1)) {
      from_state <- valid_seq[t]
      to_state <- valid_seq[t + 1]
      from_idx <- which(all_states == from_state)
      to_idx <- which(all_states == to_state)
      trans_idx <- (from_idx - 1) * n_states + to_idx
      transition_features[i, trans_idx] <- transition_features[i, trans_idx] + 1
    }
    
    # 2. Second-order patterns (trigrams)
    if (length(valid_seq) >= 3) {
      for (t in seq_len(length(valid_seq) - 2)) {
        state1 <- valid_seq[t]
        state2 <- valid_seq[t + 1]
        state3 <- valid_seq[t + 2]
        
        # Use combined state as feature
        idx1 <- which(all_states == state1)
        idx2 <- which(all_states == state2)
        bigram_idx <- (idx1 - 1) * n_states + idx2
        second_order_features[i, bigram_idx] <- second_order_features[i, bigram_idx] + 1
      }
    }
    
    # 3. Time-weighted transitions (early vs late transitions)
    seq_length <- length(valid_seq)
    for (t in seq_len(length(valid_seq) - 1)) {
      from_state <- valid_seq[t]
      to_state <- valid_seq[t + 1]
      from_idx <- which(all_states == from_state)
      to_idx <- which(all_states == to_state)
      trans_idx <- (from_idx - 1) * n_states + to_idx
      
      # Weight: higher for later transitions
      time_weight <- t / seq_length
      time_weighted_features[i, trans_idx] <- time_weighted_features[i, trans_idx] + time_weight
    }
    
    # 4. Complexity measures
    unique_states <- length(unique(valid_seq))
    state_counts <- table(valid_seq)
    entropy <- -sum((state_counts/sum(state_counts)) * log(state_counts/sum(state_counts) + 1e-10))
    
    # Count unique transitions
    unique_transitions <- length(unique(paste(valid_seq[-length(valid_seq)], 
                                             valid_seq[-1], sep = "->")))
    
    # Count repetitions (same state consecutively)
    repetitions <- sum(valid_seq[-1] == valid_seq[-length(valid_seq)])
    
    complexity_features[i, ] <- c(entropy, unique_transitions, repetitions)
    
    # 5. Persistence features (time spent in each state)
    for (state in all_states) {
      state_idx <- which(all_states == state)
      persistence_features[i, state_idx] <- sum(valid_seq == state) / length(valid_seq)
    }
  }
  
  # Normalize features
  transition_features <- transition_features / rowSums(transition_features + 1e-10)
  second_order_features <- second_order_features / rowSums(second_order_features + 1e-10)
  time_weighted_features <- time_weighted_features / rowSums(time_weighted_features + 1e-10)
  
  # Normalize complexity features to [0, 1]
  for (j in seq_len(ncol(complexity_features))) {
    range_val <- max(complexity_features[, j]) - min(complexity_features[, j])
    if (range_val > 0) {
      complexity_features[, j] <- (complexity_features[, j] - min(complexity_features[, j])) / range_val
    }
  }
  
  # Combine all features with weights
  combined_features <- cbind(
    weights[1] * transition_features,
    weights[2] * second_order_features,
    weights[3] * time_weighted_features,
    weights[4] * complexity_features,
    weights[5] * persistence_features
  )
  
  # Replace any NaN with 0
  combined_features[is.nan(combined_features)] <- 0
  
  stats::dist(combined_features, method = "euclidean")
}

#' Compute Optimal Matching Distance for Sequences
#'
#' Implements optimal matching algorithm with dynamic programming.
#' Allows for substitution and indel costs.
#'
#' @param data A data frame or matrix where rows are sequences and columns are time points
#' @param substitution_cost Cost for substituting one state for another (default: 1)
#' @param indel_cost Cost for insertions and deletions (default: 1)
#' @return A distance object of class 'dist'
#' @keywords internal
compute_optimal_matching_distance <- function(data, substitution_cost = 1, indel_cost = 1) {
  if (!is.data.frame(data) && !is.matrix(data)) {
    stop("Data must be a data frame or matrix")
  }
  
  # Convert to character matrix
  char_matrix <- as.matrix(data)
  mode(char_matrix) <- "character"
  
  # Remove NAs and create clean sequences
  sequences <- lapply(seq_len(nrow(char_matrix)), function(i) {
    row_data <- char_matrix[i, ]
    row_data[!is.na(row_data)]
  })
  
  n_seq <- length(sequences)
  dist_matrix <- matrix(0, nrow = n_seq, ncol = n_seq)
  
  # Optimal matching algorithm
  .optimal_matching <- function(seq1, seq2) {
    n1 <- length(seq1)
    n2 <- length(seq2)
    
    if (n1 == 0) return(n2 * indel_cost)
    if (n2 == 0) return(n1 * indel_cost)
    
    # Dynamic programming matrix
    dp <- matrix(0, nrow = n1 + 1, ncol = n2 + 1)
    
    # Initialize first row and column
    dp[1, ] <- seq(0, n2) * indel_cost
    dp[, 1] <- seq(0, n1) * indel_cost
    
    # Fill the matrix
    for (i in 2:(n1 + 1)) {
      for (j in 2:(n2 + 1)) {
        match_cost <- ifelse(seq1[i-1] == seq2[j-1], 0, substitution_cost)
        dp[i, j] <- min(
          dp[i-1, j] + indel_cost,      # deletion
          dp[i, j-1] + indel_cost,      # insertion
          dp[i-1, j-1] + match_cost     # substitution
        )
      }
    }
    
    return(dp[n1 + 1, n2 + 1])
  }
  
  # Compute pairwise distances
  for (i in seq_len(n_seq)) {
    for (j in i:n_seq) {
      if (i == j) {
        dist_matrix[i, j] <- 0
      } else {
        dist_matrix[i, j] <- dist_matrix[j, i] <- .optimal_matching(sequences[[i]], sequences[[j]])
      }
    }
  }
  
  as.dist(dist_matrix)
}

# ==============================================================================
# STRINGDIST-BASED DISTANCE METHODS
# ==============================================================================

#' Convert sequences to strings for stringdist methods
#' @param data A data frame or matrix where rows are sequences and columns are time points
#' @return Character vector of concatenated sequences
#' @keywords internal
.sequences_to_strings <- function(data) {
  char_matrix <- as.matrix(data)
  mode(char_matrix) <- "character"
  
  apply(char_matrix, 1, function(row) {
    # Remove NAs and concatenate
    clean_row <- row[!is.na(row)]
    if (length(clean_row) == 0) return("")
    paste(clean_row, collapse = "")
  })
}

#' Create pairwise distance matrix using stringdist methods
#' @param sequences Character vector of sequences
#' @param method StringDist method name
#' @param ... Additional parameters for stringdist
#' @return A distance object of class 'dist'
#' @keywords internal
.create_stringdist_matrix <- function(sequences, method, ...) {
  n_seq <- length(sequences)
  dist_matrix <- matrix(0, nrow = n_seq, ncol = n_seq)
  
  # Compute pairwise distances
  for (i in seq_len(n_seq)) {
    for (j in i:n_seq) {
      if (i == j) {
        dist_matrix[i, j] <- 0
      } else {
        dist_value <- stringdist::stringdist(sequences[i], sequences[j], method = method, ...)
        dist_matrix[i, j] <- dist_matrix[j, i] <- dist_value
      }
    }
  }
  
  as.dist(dist_matrix)
}

#' Compute OSA (Optimal String Alignment) Distance
#' @param data A data frame or matrix where rows are sequences and columns are time points
#' @return A distance object of class 'dist'
#' @keywords internal
compute_osa_distance <- function(data) {
  sequences <- .sequences_to_strings(data)
  .create_stringdist_matrix(sequences, "osa")
}

#' Compute Levenshtein Distance
#' @param data A data frame or matrix where rows are sequences and columns are time points
#' @return A distance object of class 'dist'
#' @keywords internal
compute_lv_distance <- function(data) {
  sequences <- .sequences_to_strings(data)
  .create_stringdist_matrix(sequences, "lv")
}

#' Compute Damerau-Levenshtein Distance
#' @param data A data frame or matrix where rows are sequences and columns are time points
#' @return A distance object of class 'dist'
#' @keywords internal
compute_dl_distance <- function(data) {
  sequences <- .sequences_to_strings(data)
  .create_stringdist_matrix(sequences, "dl")
}

#' Compute Q-gram Distance
#' @param data A data frame or matrix where rows are sequences and columns are time points
#' @param q Q-gram size (default: 2)
#' @return A distance object of class 'dist'
#' @keywords internal
compute_qgram_distance <- function(data, q = 2) {
  sequences <- .sequences_to_strings(data)
  .create_stringdist_matrix(sequences, "qgram", q = q)
}

#' Compute Jaro Distance
#' @param data A data frame or matrix where rows are sequences and columns are time points
#' @return A distance object of class 'dist'
#' @keywords internal
compute_jaro_distance <- function(data) {
  sequences <- .sequences_to_strings(data)
  .create_stringdist_matrix(sequences, "jaro")
}

#' Compute Jaro-Winkler Distance
#' @param data A data frame or matrix where rows are sequences and columns are time points
#' @param p Prefix scaling factor (default: 0.1)
#' @return A distance object of class 'dist'
#' @keywords internal
compute_jw_distance <- function(data, p = 0.1) {
  sequences <- .sequences_to_strings(data)
  .create_stringdist_matrix(sequences, "jw", p = p)
}

#' Compute Cosine Distance on Q-grams
#' @param data A data frame or matrix where rows are sequences and columns are time points
#' @param q Q-gram size (default: 2)
#' @return A distance object of class 'dist'
#' @keywords internal
compute_cosine_distance <- function(data, q = 2) {
  sequences <- .sequences_to_strings(data)
  .create_stringdist_matrix(sequences, "cosine", q = q)
}

#' Compute Jaccard Distance on Q-grams
#' @param data A data frame or matrix where rows are sequences and columns are time points
#' @param q Q-gram size (default: 2)
#' @return A distance object of class 'dist'
#' @keywords internal
compute_jaccard_distance <- function(data, q = 2) {
  sequences <- .sequences_to_strings(data)
  .create_stringdist_matrix(sequences, "jaccard", q = q)
}

#' Compute Soundex Distance
#' @param data A data frame or matrix where rows are sequences and columns are time points
#' @return A distance object of class 'dist'
#' @keywords internal
compute_soundex_distance <- function(data) {
  sequences <- .sequences_to_strings(data)
  .create_stringdist_matrix(sequences, "soundex")
}

#' Compute LCS Distance via stringdist
#' @param data A data frame or matrix where rows are sequences and columns are time points
#' @return A distance object of class 'dist'
#' @keywords internal
compute_lcs_stringdist_distance <- function(data) {
  sequences <- .sequences_to_strings(data)
  .create_stringdist_matrix(sequences, "lcs")
} 
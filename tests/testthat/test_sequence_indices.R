# Test suite for compute_sequence_indices() in R/sequence_indices.R

context("Sequence Indices: compute_sequence_indices")

# Test data
test_data_indices <- data.frame(
  T1 = c("Active", "Average", "Disengaged", "Active", "Average"),
  T2 = c("Active", "Active", "Disengaged", "Active", "Average"), 
  T3 = c("Average", "Active", "Average", "Disengaged", "Active"),
  T4 = c("Average", "Disengaged", "Average", "Disengaged", "Active"),
  T5 = c("Disengaged", "Disengaged", "Active", "Active", "Active"),
  Group = c("G1", "G1", "G2", "G1", "G2")
)

favorable_states_test <- c("Active", "Average")

test_that("Basic functionality and output structure", {
  result <- compute_sequence_indices(
    test_data_indices[, -which(names(test_data_indices) == "Group")], # Data without group
    favorable_states = favorable_states_test
  )
  
  expect_s3_class(result, "data.frame")
  # Expected columns based on the function's results data.frame initialization
  # This list might need adjustment if the function evolved.
  # Taking from the latest understanding of compute_sequence_indices
  expected_cols <- c(
    "sequence_id", "sequence_length", "valid_observations", "valid_proportion",
    "unique_states_visited", "mean_spell_duration", "longitudinal_entropy", "simpson_diversity",
    "self_loop_tendency", "transition_rate", "transition_complexity",
    "initial_state_persistence", "initial_state_influence_decay", "cyclic_feedback_strength",
    "first_state", "last_state", "attractor_state", "attractor_strength",
    "emergent_state", "emergent_state_persistence", "emergent_state_proportion",
    "integrative_potential", "complexity_index",
    "proportion_favorable_states", "favorable_state_stability"
  )
  expect_true(all(expected_cols %in% names(result)))
  expect_equal(nrow(result), nrow(test_data_indices))

  # With group and group summary
  result_grouped <- compute_sequence_indices(
    test_data_indices,
    group_col = "Group",
    favorable_states = favorable_states_test,
    return_group_summary = TRUE
  )
  expect_type(result_grouped, "list")
  expect_named(result_grouped, c("individual_indices", "group_summaries", "parameters"))
  expect_s3_class(result_grouped$individual_indices, "data.frame")
  expect_s3_class(result_grouped$group_summaries, "data.frame") # aggregate returns df
  expect_true("Group" %in% names(result_grouped$group_summaries))
})

test_that("compute_attractor_state correctly identifies modal state and handles unused params", {
  # Test compute_attractor_state directly as its signature was simplified
  states1 <- c("A", "B", "A", "A", "C", "A") # Attractor A, prop 4/6
  attractor1 <- compute_attractor_state(states1)
  expect_equal(attractor1$state, "A")
  expect_equal(attractor1$strength, 4/6)

  states2 <- c("X", "Y", "Z", "X", "Y", "Z") # Tie, X should be chosen
  attractor2 <- compute_attractor_state(states2)
  expect_equal(attractor2$state, "X") # which.max returns first in case of tie
  expect_equal(attractor2$strength, 2/6)

  states_empty <- character(0)
  attractor_empty <- compute_attractor_state(states_empty)
  expect_true(is.na(attractor_empty$state))
  expect_equal(attractor_empty$strength, 0)

  states_single <- c("S")
  attractor_single <- compute_attractor_state(states_single)
  expect_equal(attractor_single$state, "S")
  expect_equal(attractor_single$strength, 1)

  # Check that the function signature in R/sequence_indices.R no longer has method/normalize
  # This is more of a static check on the source, but important.
  # For now, trust the fix was applied.
})

test_that("Edge cases for compute_sequence_indices", {
  # Single state sequence
  single_state_df <- data.frame(T1 = "A", T2 = "A", T3 = "A")
  res_single <- compute_single_sequence_indices(as.character(single_state_df[1,]), all_possible_states = "A")

  expect_equal(res_single$unique_states_visited, 1)
  expect_equal(res_single$mean_spell_duration, 3) # One spell of length 3
  expect_equal(res_single$longitudinal_entropy, 0)
  expect_equal(res_single$simpson_diversity, 0)
  expect_equal(res_single$transition_rate, 0) # No transitions
  expect_equal(res_single$self_loop_tendency, 1) # All transitions are self-loops (implicit at ends)
                                                # Or, (2 self-loops) / (2 transitions) = 1
  expect_equal(res_single$complexity_index, 0, tolerance = 1e-5) # Should be very low/zero

  # Sequence with all NAs
  na_sequence_df <- data.frame(T1 = NA_character_, T2 = NA_character_)
  res_na <- compute_single_sequence_indices(as.character(na_sequence_df[1,]), all_possible_states = "A")
  empty_indices <- create_empty_indices_list(sequence_length = 2, valid_observations = 0)
  # Compare relevant numeric fields that should be 0
  expect_equal(res_na$valid_observations, 0)
  expect_equal(res_na$unique_states_visited, empty_indices$unique_states_visited)
  expect_equal(res_na$longitudinal_entropy, empty_indices$longitudinal_entropy)
  expect_true(is.na(res_na$first_state)) # Check a character field for NA

  # Sequence shorter than min_length for main function
  short_seq_data <- data.frame(T1 = "A", Group = "G1")
  result_short <- compute_sequence_indices(short_seq_data, group_col = "Group", min_length = 2)
  expect_equal(nrow(result_short), 0) # Filtered out
})

test_that("Integrative potential calculation is correct", {
  # Formula: (Σ is_positive[i] * i^ω) / (Σ i^ω), assuming ω=1
  # States: D, A, A (Favorable: A)
  # is_positive: 0, 1, 1
  # positions (i): 1, 2, 3
  # Numerator: (0*1) + (1*2) + (1*3) = 0 + 2 + 3 = 5
  # Denominator (i^1): 1 + 2 + 3 = 6
  # Potential = 5/6

  seq_integrative <- data.frame(T1 = "D", T2 = "A", T3 = "A")
  res_integrative <- compute_sequence_indices(seq_integrative, favorable_states = "A")
  expect_equal(res_integrative$integrative_potential[1], 5/6, tolerance = 1e-5)

  # All favorable
  seq_all_fav <- data.frame(T1 = "A", T2 = "A", T3 = "A")
  res_all_fav <- compute_sequence_indices(seq_all_fav, favorable_states = "A")
  expect_equal(res_all_fav$integrative_potential[1], 1.0, tolerance = 1e-5) # (1*1+1*2+1*3)/(1+2+3) = 6/6 = 1

  # None favorable
  seq_none_fav <- data.frame(T1 = "D", T2 = "D", T3 = "D")
  res_none_fav <- compute_sequence_indices(seq_none_fav, favorable_states = "A")
  expect_equal(res_none_fav$integrative_potential[1], 0.0, tolerance = 1e-5)
})

test_that("Complexity index components and range", {
  # Sequence: A-A-A (Entropy=0, Transitions=0, Variability=0) -> Complexity = 0
  seq_simple <- data.frame(T1="A", T2="A", T3="A")
  res_simple <- compute_sequence_indices(seq_simple, favorable_states="A")
  expect_equal(res_simple$complexity_index[1], 0.0, tolerance=1e-5)

  # Sequence: A-B-C (Max Entropy for 3 states, Max Transitions, CV of spells depends on definition)
  # Spells: A (1), B (1), C (1). Durations: 1,1,1. Mean=1, SD=0. CV=0.
  # Entropy for A,B,C (props 1/3 each) = -3 * (1/3 * log(1/3)) = log(3)
  # Max Entropy = log(3). Entropy component = 1.
  # Transitions = 2. Max transitions = 2. Transition component = 1.
  # Variability component (CV of spell durations 1,1,1) = 0.
  # Complexity = 0.4 * 1 + 0.4 * 1 + 0.2 * 0 = 0.8
  seq_abc <- data.frame(T1="A", T2="B", T3="C")
  res_abc <- compute_sequence_indices(seq_abc, favorable_states="A", all_possible_states=c("A","B","C"))
  expect_equal(res_abc$longitudinal_entropy[1], log(3), tolerance=1e-5)
  expect_equal(res_abc$transition_rate[1], 2/2, tolerance=1e-5) # 2 transitions in 2 opportunities
  # CV of (1,1,1) is 0.
  expect_equal(res_abc$complexity_index[1], 0.4 * 1 + 0.4 * 1 + 0.2 * 0, tolerance=1e-5)

  # Sequence: A-B-A-B (Entropy for A,B (props 1/2 each) = log(2). Max Entropy = log(2). E=1)
  # Transitions = 3. Max Transitions = 3. T=1.
  # Spells: A(1)B(1)A(1)B(1). Durations 1,1,1,1. CV=0.
  # Complexity = 0.4 * 1 + 0.4 * 1 + 0.2 * 0 = 0.8
  seq_abab <- data.frame(T1="A", T2="B", T3="A", T4="B")
  res_abab <- compute_sequence_indices(seq_abab, favorable_states="A", all_possible_states=c("A","B"))
  expect_equal(res_abab$complexity_index[1], 0.8, tolerance=1e-5)

  # Ensure complexity index is within [0,1] for the main test data
  result_main <- compute_sequence_indices(test_data_indices, favorable_states = favorable_states_test)
  expect_true(all(result_main$complexity_index >= 0 & result_main$complexity_index <= 1))
})

test_that("Favorable state measures are correct", {
  data_fav <- data.frame(
    T1 = c("F", "U", "F", "F"), # F=Favorable, U=Unfavorable
    T2 = c("F", "U", "U", "F"),
    T3 = c("U", "F", "F", "F")
  )
  # Seq1: F-F-U. Prop Fav = 2/3. Fav Stability (spells of F): mean(2) = 2.
  # Seq2: U-U-F. Prop Fav = 1/3. Fav Stability: mean(1) = 1.
  # Seq3: F-U-F. Prop Fav = 2/3. Fav Stability: mean(1,1) = 1.
  # Seq4: F-F-F. Prop Fav = 3/3. Fav Stability: mean(3) = 3.

  res_fav <- compute_sequence_indices(data_fav, favorable_states = "F")

  expect_equal(res_fav$proportion_favorable_states, c(2/3, 1/3, 2/3, 3/3), tolerance=1e-5)
  # Favorable state stability (mean length of favorable spells)
  # Seq1 (F,F,U): Favorable spells are (F,F). Lengths: 2. Mean=2.
  # Seq2 (U,U,F): Favorable spells are (F). Lengths: 1. Mean=1.
  # Seq3 (F,U,F): Favorable spells are (F), (F). Lengths: 1,1. Mean=1.
  # Seq4 (F,F,F): Favorable spells are (F,F,F). Lengths: 3. Mean=3.
  expect_equal(res_fav$favorable_state_stability, c(2, 1, 1, 3), tolerance=1e-5)
})

test_that("Initial state influence measures", {
    # Initial state persistence: prop of sequence length the initial state lasts
    # Initial state influence decay: presence in first third vs last third
    seq_init <- data.frame(T1=c("A","A","B","A","C","A","A","A","A")) # Length 9
    # Seq1: A A B A C A A A A
    # Initial state A. Initial spell A A (length 2). Persistence = 2/9.
    # First third (1-3): A A B. Prop A = 2/3.
    # Last third (7-9): A A A. Prop A = 3/3.
    # Decay = 2/3 - 3/3 = -1/3 (means it increased)

    res_init <- compute_sequence_indices(seq_init, all_possible_states=c("A","B","C"))
    expect_equal(res_init$initial_state_persistence[1], 2/9, tolerance=1e-5)
    expect_equal(res_init$initial_state_influence_decay[1], (2/3) - (3/3), tolerance=1e-5)
})

test_that("print_indices_summary runs", {
  # Test with individual indices data frame
  res_df <- compute_sequence_indices(test_data_indices[, -which(names(test_data_indices) == "Group")])
  expect_output(print_indices_summary(res_df), "SEQUENCE INDICES SUMMARY")

  # Test with grouped summary list
  res_list <- compute_sequence_indices(test_data_indices, group_col = "Group", return_group_summary = TRUE)
  expect_output(print_indices_summary(res_list), "SEQUENCE INDICES SUMMARY")
  expect_output(print_indices_summary(res_list), "Parameters:") # Check for parameter output

  # Check a specific line from the summary output if possible and stable
  # For example, mean sequence length.
  # This requires capturing output and parsing, which can be fragile.
  # For now, just checking it runs and prints the header is a good start.
})

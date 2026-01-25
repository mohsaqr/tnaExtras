# ==============================================================================
# PATHWAY ANALYSIS FUNCTIONS
# ==============================================================================
#
# Functions for analyzing sequence pathways and their relationship to outcomes.
# Supports group_tna objects for group comparison analysis.
#
# Functions:
# - pathway_outcome(): Main analysis function for pathway-outcome relationships
# - extract_sequences_from_wide(): Extract learning sequences from wide-format data
# - extract_learning_sequences(): Extract sequences from long-format data
# - classify_by_last_state(): Classify sequences by final state
#
# ==============================================================================

# Ensure data.table knows this package is aware of its semantics
.datatable.aware <- TRUE

#' @importFrom data.table := .SD .BY .N .I
NULL

#' Classify Wide Data by Last Valid State
#'
#' This function looks at a set of columns (wide format), ignores NAs and empty
#' strings, and returns the last valid state for each row.
#'
#' @param data A data.frame or data.table.
#' @param cols A character vector of column names or integer indices representing
#'   the sequence columns.
#'
#' @return A character vector of classifications (the last non-NA, non-empty
#'   value for each row).
#'
#' @examples
#' df <- data.frame(
#'   T1 = c("Wrong", "Step", "Wrong"),
#'   T2 = c("Step", "Right", "Step"),
#'   T3 = c("Right", NA, "Step"),
#'   T4 = c(NA, NA, "Right")
#' )
#' classify_by_last_state(df, c("T1", "T2", "T3", "T4"))
#'
#' @export
classify_by_last_state <- function(data, cols) {
  # Handle data.table column subsetting

  if (inherits(data, "data.table")) {
    state_subset <- data[, cols, with = FALSE]
  } else {
    state_subset <- data[, cols, drop = FALSE]
  }

  wide_matrix <- as.matrix(state_subset)

  # Apply row-wise logic: find last non-NA, non-empty value

  group <- apply(wide_matrix, 1, function(row) {
    valid_values <- row[!is.na(row) & row != ""]
    if (length(valid_values) == 0) {
      return(NA_character_)
    }
    return(valid_values[length(valid_values)])
  })

  return(as.character(group))
}


#' Extract Learning Sequences from Long-Format Data
#'
#' Extracts learning sequences from long-format event data where each row is an
#' event. Sequences start with a specified start state and end with one of the
#' specified end states.
#'
#' @param data A data.frame or data.table containing event logs.
#' @param state_col Character. The name of the column containing event types
#'   (e.g., "Event_Type").
#' @param group_cols Character vector. Column names to group by when extracting
#'   sequences (e.g., c("student_id", "skill_id")).
#' @param start_state Character. The state that initiates a sequence
#'   (default is "Wrong").
#' @param end_states Character vector. States that terminate a sequence
#'   (default is c("Correct", "Wrong")).
#'
#' @return A data.table containing the extracted sequences with added columns:
#'   \item{sequence_id}{Unique identifier for each sequence}
#'   \item{sequence_end_type}{The state that ended this sequence}
#'
#' @details
#' The function processes event logs and extracts sequences that:
#' \itemize{
#'   \item Start when the specified start_state is encountered
#'   \item End when any of the end_states is encountered
#'   \item Include all events between start and end
#' }
#'
#' This is useful for analyzing learning trajectories where you want to
#' understand what happens between an initial error and the eventual outcome.
#'
#' @examples
#' # Create sample long-format data
#' events <- data.frame(
#'   student_id = c(1, 1, 1, 1, 1, 2, 2, 2),
#'   Event_Type = c("Wrong", "Hint", "Step", "Correct", "Wrong", "Wrong", "Step", "Correct"),
#'   timestamp = 1:8
#' )
#' sequences <- extract_learning_sequences(
#'   events,
#'   state_col = "Event_Type",
#'   group_cols = "student_id",
#'   start_state = "Wrong",
#'   end_states = c("Correct", "Wrong")
#' )
#'
#' @export
extract_learning_sequences <- function(data,
                                       state_col = "Event_Type",
                                       group_cols = c("studentDiagnosticsId", "SkillID"),
                                       start_state = "Wrong",
                                       end_states = c("Correct", "Wrong")) {

  # Handle group_tna input
  if (is_group_tna(data)) {
    converted <- convert_group_tna(data)
    data <- converted$data
    # Add Group to group_cols if not already present
    if (!converted$group_col %in% group_cols) {
      group_cols <- c(group_cols, converted$group_col)
    }
  }

  dt <- data.table::as.data.table(data)

  # Internal helper to process a single group
  process_group <- function(sub_dt, group_vals) {
    states <- sub_dt[[state_col]]

    # Identify outcome boundaries
    is_result <- states %in% end_states
    is_start <- states == start_state

    res_idx <- which(is_result)
    start_idx_list <- which(is_start)

    sequences <- list()

    if (length(start_idx_list) == 0) {
      return(NULL)
    }

    for (s_idx in start_idx_list) {
      # Find first result index that is > s_idx
      e_idx <- res_idx[res_idx > s_idx][1]

      if (!is.na(e_idx)) {
        seg <- sub_dt[s_idx:e_idx]

        # Create a unique ID using grouping values and start position
        group_id_str <- paste(group_vals, collapse = "_")
        seg[, sequence_id := paste0(group_id_str, "_seq_", s_idx)]
        seg[, sequence_end_type := states[e_idx]]

        sequences[[length(sequences) + 1]] <- data.table::copy(seg)
      }
    }
    return(data.table::rbindlist(sequences))
  }

  # Apply the helper grouped by the specified columns
  result <- dt[, process_group(.SD, .BY), by = group_cols]

  return(result)
}


#' Extract Learning Sequences from Wide-Format Data
#'
#' Converts wide-format sequence data (one row per sequence) into individual
#' learning sequences, where each sequence starts with "Wrong" and ends with
#' either "Right" or "Wrong".
#'
#' @param data A data.frame or data.table in wide format with columns T1, T2, etc.
#' @param id_col Character. Name of the ID column. If NULL, auto-detects as the
#'   first non-T column.
#' @param t_col_pattern Character. Regex pattern to identify time columns
#'   (default is "^T[0-9]+$").
#' @param start_state Character. State that initiates a new sequence
#'   (default is "Wrong").
#' @param end_states Character vector. States that terminate a sequence
#'   (default is c("Right", "Wrong")).
#'
#' @return A data.table with extracted sequences in wide format, including:
#'   \item{id_col}{Original ID column}
#'   \item{seq_id}{Sequence number within each ID}
#'   \item{T1, T2, ...}{Renumbered time steps for each sequence}
#'
#' @details
#' The function:
#' \enumerate{
#'   \item Melts wide data to long format
#'   \item Identifies sequence boundaries based on start/end states
#'   \item Handles overlapping sequences (when one ends with start_state)
#'   \item Casts back to wide format with renumbered time steps
#' }
#'
#' This is useful for pre-processing data before pathway analysis when you have
#' wide-format data with multiple learning attempts per row.
#'
#' @examples
#' df <- data.frame(
#'   student_id = 1:2,
#'   T1 = c("Wrong", "Wrong"),
#'   T2 = c("Step", "Step"),
#'   T3 = c("Right", "Wrong"),
#'   T4 = c(NA, "Step"),
#'   T5 = c(NA, "Right")
#' )
#' extract_sequences_from_wide(df)
#'
#' @export
extract_sequences_from_wide <- function(data,
                                        id_col = NULL,
                                        t_col_pattern = "^T[0-9]+$",
                                        start_state = "Wrong",
                                        end_states = c("Right", "Wrong")) {

  # Handle group_tna input
  group_col <- NULL
  if (is_group_tna(data)) {
    converted <- convert_group_tna(data)
    data <- converted$data
    group_col <- converted$group_col
  }

  dt_wide <- data.table::as.data.table(data)

  # Identify T columns automatically
  t_cols <- grep(t_col_pattern, names(dt_wide), value = TRUE)

  # Sort T columns numerically to ensure correct order
  t_nums <- as.integer(gsub("^T", "", t_cols))
  t_cols <- t_cols[order(t_nums)]

  # Identify ID column if not specified
  if (is.null(id_col)) {
    non_t_cols <- setdiff(names(dt_wide), t_cols)
    # Exclude group column from auto-detection
    if (!is.null(group_col)) {
      non_t_cols <- setdiff(non_t_cols, group_col)
    }
    id_col <- non_t_cols[1]
  }

  # Build id.vars for melt
  id_vars <- id_col
  if (!is.null(group_col)) {
    id_vars <- c(id_vars, group_col)
  }

  # Melt to long format - PRESERVE ORIGINAL ORDER
  dt_long <- data.table::melt(dt_wide,
                              id.vars = id_vars,
                              measure.vars = t_cols,
                              variable.name = "time_step",
                              value.name = "Event",
                              na.rm = TRUE)

  # Extract numeric time for proper ordering
  dt_long[, time_num := as.integer(gsub("^T", "", time_step))]
  data.table::setorderv(dt_long, c(id_col, "time_num"))

  # Identify sequence boundaries
  dt_long[, is_start := Event == start_state]
  dt_long[, is_end := Event %in% end_states]

  # Assign sequence numbers based on starts
  dt_long[, seq_num := 0]
  dt_long[is_start == TRUE, seq_num := 1:.N, by = id_col]

  # Forward fill sequence numbers
  dt_long[, seq_id := cummax(seq_num), by = id_col]

  # Keep only rows that belong to a sequence (seq_id > 0)
  dt_seq <- dt_long[seq_id > 0]

  # For each sequence, keep events until (and including) the first end event
  dt_seq[, cumul_ends := cumsum(is_end), by = c(id_col, "seq_id")]
  dt_seq <- dt_seq[cumul_ends <= 1]

  # Handle overlap: if sequence ends on start_state, duplicate for next sequence
  dt_seq[, ends_on_start := is_end & Event == start_state]

  # Identify which sequences need their last event duplicated
  overlap_info <- dt_seq[, .(
    last_time = max(time_num),
    ends_start = any(ends_on_start)
  ), by = c(id_col, "seq_id")]

  overlap_info <- overlap_info[ends_start == TRUE]

  if (nrow(overlap_info) > 0) {
    # For each overlapping sequence, duplicate the last event into next sequence
    overlap_rows <- dt_seq[overlap_info, on = c(id_col, "seq_id"), nomatch = 0]
    overlap_rows <- overlap_rows[, .SD[.N], by = c(id_col, "seq_id")]
    overlap_rows[, seq_id := seq_id + 1]

    # Add duplicated rows
    dt_seq <- data.table::rbindlist(list(dt_seq, overlap_rows), fill = TRUE)
  }

  # Re-sort by ORIGINAL time order
  data.table::setorderv(dt_seq, c(id_col, "seq_id", "time_num"))

  # Create position within sequence
  dt_seq[, step_in_seq := 1:.N, by = c(id_col, "seq_id")]
  dt_seq[, time_step_new := paste0("T", step_in_seq)]

  # Build cast formula
  if (!is.null(group_col)) {
    cast_formula <- stats::as.formula(paste(id_col, "+", group_col, "+ seq_id ~ time_step_new"))
  } else {
    cast_formula <- stats::as.formula(paste(id_col, "+ seq_id ~ time_step_new"))
  }

  # Cast back to wide format
  result <- data.table::dcast(dt_seq, cast_formula, value.var = "Event")

  # Force numerical sorting of T columns
  t_result_cols <- grep("^T[0-9]+$", names(result), value = TRUE)
  t_result_nums <- as.integer(gsub("^T", "", t_result_cols))
  t_result_cols_sorted <- t_result_cols[order(t_result_nums)]

  other_cols <- setdiff(names(result), t_result_cols)
  data.table::setcolorder(result, c(other_cols, t_result_cols_sorted))

  return(result)
}


#' Analyze Pathway-Outcome Relationships
#'
#' Analyzes sequence pathways and their relationship to outcomes, with optional
#' group comparison. Provides statistical analysis of action presence, exact
#' sequences, and logistic regression modeling.
#'
#' @param data A data.frame, data.table, or group_tna object. If group_tna,
#'   groups are automatically extracted for comparison.
#' @param seq_cols Character vector. Names of sequence columns. If NULL,
#'   auto-detects columns matching "^T\\d+" pattern.
#' @param id_col Character. Name of the ID column. If NULL, uses first column.
#' @param outcome_col Character. Name of outcome column. If NULL, uses last
#'   non-NA value in each sequence as the outcome.
#' @param group_col Character. Name of group column for comparison analysis.
#'   If using group_tna input, this is auto-populated.
#' @param ignore_initial_values Character vector. States to ignore if they
#'   appear as the first value (default is c("Wrong")).
#' @param min_count Integer. Minimum count threshold for including patterns
#'   in analysis (default is 100).
#' @param min_seq_length Integer. Minimum sequence length to include
#'   (default is 1).
#' @param max_depth Integer. Maximum sequence depth to analyze (default is 7).
#' @param top_n_patterns Integer. Number of top patterns to include in
#'   detailed analysis and plots (default is 20).
#' @param aggregate_rare Logical. Whether to aggregate rare patterns into
#'   "Other" category in GLM (default is TRUE).
#' @param run_glm Logical. Whether to run logistic regression models
#'   (default is TRUE).
#' @param success_class Character. The outcome value considered as success
#'   (default is "Right").
#' @param verbose Logical. Whether to print progress and results
#'   (default is TRUE).
#' @param plot Logical. Whether to generate plots (default is TRUE).
#'
#' @return A list containing:
#'   \item{tables}{List of result tables:
#'     \itemize{
#'       \item{contains}{Action presence analysis with chi-square tests}
#'       \item{exact}{Exact sequence analysis with chi-square tests}
#'     }
#'   }
#'   \item{plots}{List of generated plots}
#'   \item{models}{List of fitted GLM models (if run_glm=TRUE):
#'     \itemize{
#'       \item{freq}{Model with action frequency predictors}
#'       \item{has}{Model with action presence (binary) predictors}
#'       \item{pattern}{Model with top pattern predictors}
#'     }
#'   }
#'   \item{status}{Named vector of analysis status for each component}
#'   \item{data}{Processed data used in analysis}
#'
#' @details
#' The function performs three main analyses:
#'
#' \strong{1. Action Presence Analysis ("contains")}
#' Tests whether the presence of each action in a sequence is associated with
#' success. Uses chi-square tests with FDR correction.
#'
#' \strong{2. Exact Sequence Analysis ("exact")}
#' Analyzes complete sequences (pathways) and their success rates. Identifies
#' which exact pathways lead to better outcomes.
#'
#' \strong{3. GLM Analysis (optional)}
#' Fits logistic regression models:
#' \itemize{
#'   \item Model 1: Action counts as predictors
#'   \item Model 2: Action presence (binary) as predictors
#'   \item Model 3: Top patterns as categorical predictors
#' }
#'
#' \strong{Group Comparison}
#' When a group_tna object or group_col is provided, the analysis includes:
#' \itemize{
#'   \item Group variable in GLM models as predictor
#'   \item Separate analysis tables by group
#'   \item Visual comparison across groups
#' }
#'
#' @examples
#' # Create sample data
#' df <- data.frame(
#'   id = 1:100,
#'   T1 = sample(c("Wrong", "Step"), 100, replace = TRUE),
#'   T2 = sample(c("Step", "Hint"), 100, replace = TRUE),
#'   T3 = sample(c("Step", "Right", "Wrong"), 100, replace = TRUE),
#'   T4 = sample(c("Right", "Wrong", NA), 100, replace = TRUE)
#' )
#'
#' # Basic analysis
#' result <- pathway_outcome(df, success_class = "Right", min_count = 5,
#'                           verbose = FALSE, plot = FALSE)
#'
#' @export
pathway_outcome <- function(data,
                            seq_cols = NULL,
                            id_col = NULL,
                            outcome_col = NULL,
                            group_col = NULL,
                            ignore_initial_values = c("Wrong"),
                            min_count = 100,
                            min_seq_length = 1,
                            max_depth = 7,
                            top_n_patterns = 20,
                            aggregate_rare = TRUE,
                            run_glm = TRUE,
                            success_class = "Right",
                            verbose = TRUE,
                            plot = TRUE) {

  status_report <- c(Preprocessing = "Pending", Contains = "Pending",
                     Exact = "Pending", GLM = "Skipped")
  log_msg <- function(...) if (verbose) cat(...)

  # === A. PREPROCESSING ===
  log_msg("\n=== PREPROCESSING ===\n")

  # Handle group_tna input
  from_group_tna <- FALSE
  if (is_group_tna(data)) {
    converted <- convert_group_tna(data)
    data <- converted$data
    group_col <- converted$group_col
    from_group_tna <- TRUE
    log_msg("Detected group_tna object with groups:",
            paste(converted$group_info$levels, collapse = ", "), "\n")
  }

  # Convert to data.table
  dt <- data.table::as.data.table(data)

  # 1. Auto-detection of columns
  # First detect seq_cols to avoid picking a sequence column as ID
  if (is.null(seq_cols)) {
    seq_cols <- grep("^T\\d+", names(dt), value = TRUE)
  }

  if (is.null(id_col)) {
    # For group_tna or if no obvious ID, create a synthetic row ID
    if (from_group_tna || length(seq_cols) > 0 && names(dt)[1] %in% seq_cols) {
      dt[, ..row_id := .I]
      id_col <- "..row_id"
      log_msg("Created synthetic ID column\n")
    } else {
      id_col <- names(dt)[1]
      log_msg("Auto-detected ID:", id_col, "\n")
    }
  }

  if (length(seq_cols) == 0) {
    exclude_cols <- c(id_col, outcome_col, group_col)
    seq_cols <- setdiff(names(dt), exclude_cols)
  }
  log_msg("Auto-detected sequence columns:", length(seq_cols), "columns\n")

  # 2. Add row identifier
  dt[, ..rid := .I]

  # 3. Pivot to long format using melt
  cols_to_keep <- c("..rid", id_col, seq_cols)
  if (!is.null(outcome_col)) cols_to_keep <- c(cols_to_keep, outcome_col)
  if (!is.null(group_col)) cols_to_keep <- c(cols_to_keep, group_col)

  dt_subset <- dt[, cols_to_keep, with = FALSE]

  id_vars <- c("..rid", id_col)
  if (!is.null(outcome_col)) id_vars <- c(id_vars, outcome_col)
  if (!is.null(group_col)) id_vars <- c(id_vars, group_col)

  long_data <- data.table::melt(dt_subset,
                                id.vars = id_vars,
                                measure.vars = seq_cols,
                                variable.name = "time_step",
                                value.name = "act",
                                na.rm = TRUE)

  # Remove empty values
  long_data <- long_data[act != ""]

  # 4. Clean Initial Values
  if (!is.null(ignore_initial_values) && length(ignore_initial_values) > 0) {
    # Add step index within each sequence
    long_data[, step_idx := seq_len(.N), by = ..rid]
    # Remove ignored initial values
    long_data <- long_data[!(step_idx == 1 & act %in% ignore_initial_values)]
    long_data[, step_idx := NULL]
  }

  # 5. Outcome Construction
  if (!is.null(outcome_col)) {
    # Use provided outcome column - it should already be in long_data
    data.table::setnames(long_data, outcome_col, "outcome", skip_absent = TRUE)
  } else {
    # Use last value as outcome
    long_data[, outcome := act[.N], by = ..rid]
  }

  long_data <- long_data[!is.na(outcome)]

  # Validate Success Class
  if (!any(long_data$outcome == success_class)) {
    warning(paste("Success class", success_class, "not found. GLM may fail."))
    status_report["Preprocessing"] <- "Warning: No Success Class"
  } else {
    status_report["Preprocessing"] <- "Success"
  }

  # 6. Collapse to pattern level
  # Build by columns
  by_cols <- c(id_col, "..rid", "outcome")
  if (!is.null(group_col)) by_cols <- c(by_cols, group_col)

  pattern_df <- long_data[, {
    acts <- act
    n_acts <- length(acts)
    pattern_full <- paste(acts, collapse = "->")

    # Body: if outcome from last act, exclude it; otherwise include all
    if (is.null(outcome_col)) {
      body_raw <- if (n_acts > 1) paste(head(acts, -1), collapse = "->") else ""
    } else {
      body_raw <- paste(acts, collapse = "->")
    }

    list(
      pattern_full = pattern_full,
      body_raw = body_raw,
      len = n_acts
    )
  }, by = by_cols]

  # Clean body
  if (!is.null(ignore_initial_values) && length(ignore_initial_values) > 0) {
    ignore_pattern <- paste0("->(", paste(ignore_initial_values, collapse = "|"), ")$")
    pattern_df[, body := gsub(ignore_pattern, "", body_raw)]
  } else {
    pattern_df[, body := body_raw]
  }

  pattern_df[body == "", body := "(Immediate Outcome)"]
  pattern_df[, depth := lengths(strsplit(body, "->", fixed = TRUE))]
  pattern_df <- pattern_df[len > 0]

  # Calculate global success rate
  succ_rate <- mean(pattern_df$outcome == success_class, na.rm = TRUE)
  log_msg(sprintf("Analyzed %d sequences. Global Success Rate (%s): %.2f%%\n",
                  nrow(pattern_df), success_class, succ_rate * 100))

  # === B. TABLES & PLOTS ===
  plots <- list()
  contains_tbl <- NULL
  seq_df <- NULL

  # 1. CONTAINS ANALYSIS (Action Presence)
  log_msg("\n--- TABLE 1: ACTION PRESENCE ---\n")
  valid_actions <- setdiff(unique(long_data$act), c(success_class, ignore_initial_values))

  if (length(valid_actions) > 0) {
    presence_data <- long_data[act %in% valid_actions]

    # Get distinct (..rid, act, outcome) combinations
    if (!is.null(group_col)) {
      presence_distinct <- unique(presence_data[, c("..rid", "act", "outcome", group_col), with = FALSE])
    } else {
      presence_distinct <- unique(presence_data[, .(..rid, act, outcome)])
    }

    # Summarize by action
    presence_df <- presence_distinct[, .(
      Count = .N,
      Success = sum(outcome == success_class, na.rm = TRUE),
      Failure = sum(outcome != success_class, na.rm = TRUE)
    ), by = act]

    presence_df[, Rate := ifelse(Count > 0, Success / Count, 0)]
    presence_df <- presence_df[Count >= min_count]

    if (nrow(presence_df) > 0) {
      # Chi-square tests
      presence_df[, c("ChiSq", "P_Value") := {
        res <- tryCatch({
          test <- stats::chisq.test(c(Success, Failure), p = c(succ_rate, 1 - succ_rate))
          list(as.numeric(test$statistic), test$p.value)
        }, error = function(e) list(NA_real_, NA_real_))
        list(res[[1]], res[[2]])
      }, by = act]

      # Filter and adjust p-values
      presence_df <- presence_df[!is.na(P_Value)]
      if (nrow(presence_df) > 0) {
        presence_df[, Adj_P_Value := stats::p.adjust(P_Value, method = "fdr")]
        data.table::setorder(presence_df, -Rate)
        data.table::setnames(presence_df, "act", "Action")

        # Remove test_res if exists and select final columns
        contains_tbl <- presence_df[, .(Action, Count, Success, Failure, Rate,
                                        ChiSq, P_Value, Adj_P_Value)]

        if (verbose) {
          print_df <- as.data.frame(contains_tbl)
          print_df[sapply(print_df, is.numeric)] <-
            lapply(print_df[sapply(print_df, is.numeric)], round, 4)
          print(print_df)
        }

        # Plot
        if (plot && nrow(contains_tbl) > 0) {
          plots$contains <- .plot_action_presence(contains_tbl, succ_rate)
        }

        status_report["Contains"] <- "Success"
      } else {
        status_report["Contains"] <- "Empty"
      }
    } else {
      status_report["Contains"] <- "Empty"
    }
  } else {
    status_report["Contains"] <- "Empty"
  }

  # 2. EXACT SEQUENCE ANALYSIS
  log_msg("\n--- TABLE 2: EXACT SEQUENCES ---\n")

  seq_stats <- pattern_df[depth >= min_seq_length & depth <= max_depth,
                          .(Count = .N,
                            Success = sum(outcome == success_class, na.rm = TRUE),
                            Failure = sum(outcome != success_class, na.rm = TRUE)),
                          by = .(body, depth)]

  seq_stats <- seq_stats[Count >= min_count]
  seq_stats[, Rate := ifelse(Count > 0, Success / Count, 0)]

  if (nrow(seq_stats) > 0) {
    # Chi-square tests
    seq_stats[, c("ChiSq", "P_Value") := {
      res <- tryCatch({
        test <- stats::chisq.test(c(Success, Failure), p = c(succ_rate, 1 - succ_rate))
        list(as.numeric(test$statistic), test$p.value)
      }, error = function(e) list(NA_real_, NA_real_))
      list(res[[1]], res[[2]])
    }, by = body]

    # Filter and adjust
    seq_stats <- seq_stats[!is.na(P_Value)]
    if (nrow(seq_stats) > 0) {
      seq_stats[, Adj_P_Value := stats::p.adjust(P_Value, method = "fdr")]
      data.table::setnames(seq_stats, c("body", "depth"), c("Sequence", "Length"))
      seq_df <- seq_stats[, .(Sequence, Length, Count, Success, Failure, Rate,
                              ChiSq, P_Value, Adj_P_Value)]
      data.table::setorder(seq_df, -Count)

      if (verbose) {
        print_df <- as.data.frame(head(seq_df, 20))
        print_df[sapply(print_df, is.numeric)] <-
          lapply(print_df[sapply(print_df, is.numeric)], round, 4)
        print(print_df)
      }

      # Plots
      if (plot) {
        # Top overall sequences
        top_overall <- head(seq_df, top_n_patterns)
        if (nrow(top_overall) > 0) {
          plots$overall_top <- .plot_sequences(top_overall,
                                               paste("Top", top_n_patterns, "Most Frequent Sequences"))
        }

        # By length
        lens <- sort(unique(seq_df$Length))
        for (l in lens) {
          p_data <- head(seq_df[Length == l], 10)
          if (nrow(p_data) > 0) {
            plots[[paste0("len_", l)]] <- .plot_sequences(p_data, paste("Length", l, "Outcomes"))
          }
        }
      }

      status_report["Exact"] <- "Success"
    } else {
      status_report["Exact"] <- "Empty"
    }
  } else {
    status_report["Exact"] <- "Empty"
  }

  # === C. GLM ANALYSIS ===
  glm_res <- list()

  if (run_glm && !is.null(contains_tbl) && nrow(contains_tbl) > 0) {
    status_report["GLM"] <- "Running"
    log_msg("\n--- GLM ANALYSIS ---\n")

    glm_df <- data.table::copy(pattern_df)
    glm_df[, y := as.integer(outcome == success_class)]

    # Get action counts per sequence
    act_counts <- long_data[act %in% contains_tbl$Action,
                            .N,
                            by = .(..rid, act)]
    act_counts_wide <- data.table::dcast(act_counts, ..rid ~ act, value.var = "N", fill = 0)

    # Merge with glm_df
    glm_df <- merge(glm_df, act_counts_wide, by = "..rid", all.x = TRUE)

    # Fill NAs with 0 for action columns
    pred_cols <- intersect(names(glm_df), contains_tbl$Action)
    for (col in pred_cols) {
      glm_df[is.na(get(col)), (col) := 0]
    }

    # Sanitize column names for formula (handle special chars)
    safe_names <- make.names(pred_cols)
    if (!identical(safe_names, pred_cols)) {
      for (i in seq_along(pred_cols)) {
        if (pred_cols[i] != safe_names[i]) {
          data.table::setnames(glm_df, pred_cols[i], safe_names[i])
        }
      }
      pred_cols <- safe_names
    }

    # Model 1: Action frequency
    if (length(pred_cols) > 0) {
      f1_terms <- pred_cols
      if (!is.null(group_col)) f1_terms <- c(f1_terms, group_col)
      f1 <- stats::as.formula(paste("y ~", paste(f1_terms, collapse = " + ")))

      glm_res$freq <- tryCatch({
        m <- stats::glm(f1, data = glm_df, family = stats::binomial)
        if (verbose) .print_glm_results(m, "Model 1: Action Frequency")
        m
      }, error = function(e) {
        log_msg("  Error in Model 1:", e$message, "\n")
        NULL
      })
    }

    # Model 2: Action presence (binary)
    if (length(pred_cols) > 0) {
      has_cols <- paste0("has_", pred_cols)
      for (i in seq_along(pred_cols)) {
        glm_df[, (has_cols[i]) := as.integer(get(pred_cols[i]) > 0)]
      }

      f2_terms <- has_cols
      if (!is.null(group_col)) f2_terms <- c(f2_terms, group_col)
      f2 <- stats::as.formula(paste("y ~", paste(f2_terms, collapse = " + ")))

      glm_res$has <- tryCatch({
        m <- stats::glm(f2, data = glm_df, family = stats::binomial)
        if (verbose) .print_glm_results(m, "Model 2: Action Presence")
        m
      }, error = function(e) {
        log_msg("  Error in Model 2:", e$message, "\n")
        NULL
      })
    }

    # Model 3: Top patterns
    top_pats <- pattern_df[, .N, by = body][order(-N)][1:min(top_n_patterns, .N), body]

    if (length(top_pats) > 0) {
      if (aggregate_rare) {
        glm_df[, pat_group := ifelse(body %in% top_pats, body, "Other")]
        ref_lvl <- if ("Other" %in% glm_df$pat_group) "Other" else
          names(sort(table(glm_df$pat_group), decreasing = TRUE))[1]
      } else {
        glm_df[, pat_group := ifelse(body %in% top_pats, body, "Baseline")]
        ref_lvl <- "Baseline"
      }
      glm_df[, pat_group := factor(pat_group)]
      glm_df[, pat_group := stats::relevel(pat_group, ref = ref_lvl)]

      if (data.table::uniqueN(glm_df$pat_group) > 1) {
        f3_terms <- "pat_group"
        if (!is.null(group_col)) f3_terms <- c(f3_terms, group_col)
        f3 <- stats::as.formula(paste("y ~", paste(f3_terms, collapse = " + ")))

        glm_res$pattern <- tryCatch({
          m <- stats::glm(f3, data = glm_df, family = stats::binomial)
          if (verbose) .print_glm_results(m, "Model 3: Top Patterns")
          m
        }, error = function(e) {
          log_msg("  Error in Model 3:", e$message, "\n")
          NULL
        })
      }
    }

    status_report["GLM"] <- "Success"
  }

  return(list(
    tables = list(contains = contains_tbl, exact = seq_df),
    plots = plots,
    models = glm_res,
    status = status_report,
    data = list(pattern_df = pattern_df, long_data = long_data)
  ))
}


# ==============================================================================
# INTERNAL HELPER FUNCTIONS
# ==============================================================================

#' Print GLM results (internal)
#' @param model GLM model object
#' @param title Title for output
#' @noRd
.print_glm_results <- function(model, title) {
  cat("\n", title, "\n")

  coef_summary <- summary(model)$coefficients
  if (is.null(coef_summary) || nrow(coef_summary) == 0) {
    cat("  No coefficients to display\n")
    return(invisible(NULL))
  }

  # Create results data.frame
  res <- data.frame(
    Term = rownames(coef_summary),
    Estimate = coef_summary[, "Estimate"],
    Std_Error = coef_summary[, "Std. Error"],
    Statistic = coef_summary[, "z value"],
    P_Value = coef_summary[, "Pr(>|z|)"],
    stringsAsFactors = FALSE
  )

  # Add odds ratio
  res$OR <- exp(res$Estimate)

  # Remove intercept and clean names
  res <- res[res$Term != "(Intercept)", ]
  res$Term <- gsub("`", "", res$Term)

  if (nrow(res) > 0) {
    # Round and print
    res[sapply(res, is.numeric)] <- lapply(res[sapply(res, is.numeric)], round, 3)
    res <- res[order(-res$OR), ]
    print(res, row.names = FALSE)
  }

  # Print fit statistics
  cat("\n  AIC:", round(stats::AIC(model), 2), "\n")

  invisible(NULL)
}


#' Plot action presence results (internal)
#' @param tbl Data table with action presence results
#' @param succ_rate Global success rate for reference line
#' @noRd
.plot_action_presence <- function(tbl, succ_rate) {
  # Order by Rate
  tbl <- as.data.frame(tbl)
  tbl <- tbl[order(tbl$Rate), ]

  # Set up colors
  colors <- ifelse(tbl$Rate >= succ_rate, "forestgreen", "grey60")

  # Create plot
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par))

  graphics::par(mar = c(5, 10, 4, 2) + 0.1)

  bp <- graphics::barplot(
    tbl$Rate,
    horiz = TRUE,
    names.arg = tbl$Action,
    las = 1,
    col = colors,
    border = NA,
    xlim = c(0, max(tbl$Rate, succ_rate) * 1.1),
    main = "Success Rate by Action Presence",
    xlab = "Success Rate"
  )

  # Add reference line
  graphics::abline(v = succ_rate, lty = 2, col = "red", lwd = 2)

  # Add rate labels
  graphics::text(tbl$Rate, bp, labels = sprintf("%.1f%%", tbl$Rate * 100),
                 pos = 4, cex = 0.8)

  invisible(NULL)
}


#' Plot sequence results (internal)
#' @param tbl Data table with sequence results
#' @param title Plot title
#' @noRd
.plot_sequences <- function(tbl, title) {
  tbl <- as.data.frame(tbl)

  # Order by Count (ascending for bottom-to-top display)
  tbl <- tbl[order(tbl$Count), ]

  # Calculate percentages
  tbl$Pct_Success <- tbl$Success / tbl$Count
  tbl$Pct_Failure <- tbl$Failure / tbl$Count

  # Set up plot
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par))

  # Calculate left margin based on longest sequence name
  max_chars <- max(nchar(as.character(tbl$Sequence)))
  left_mar <- min(max(max_chars * 0.5, 8), 20)

  graphics::par(mar = c(5, left_mar, 4, 2) + 0.1)

  n <- nrow(tbl)

  # Create stacked barplot
  bar_matrix <- rbind(tbl$Failure, tbl$Success)
  rownames(bar_matrix) <- c("Failure", "Success")

  bp <- graphics::barplot(
    bar_matrix,
    horiz = TRUE,
    names.arg = tbl$Sequence,
    las = 1,
    col = c("#b22222", "#228b22"),
    border = NA,
    main = title,
    xlab = "Count",
    cex.names = 0.7
  )

  # Add percentage labels in center of each segment
  for (i in 1:n) {
    # Failure label (left side)
    if (tbl$Failure[i] > 0) {
      x_fail <- tbl$Failure[i] / 2
      graphics::text(x_fail, bp[i],
                     labels = sprintf("%d%%", round(tbl$Pct_Failure[i] * 100)),
                     col = "white", cex = 0.7, font = 2)
    }
    # Success label (right side)
    if (tbl$Success[i] > 0) {
      x_succ <- tbl$Failure[i] + tbl$Success[i] / 2
      graphics::text(x_succ, bp[i],
                     labels = sprintf("%d%%", round(tbl$Pct_Success[i] * 100)),
                     col = "white", cex = 0.7, font = 2)
    }
  }

  # Add legend
  graphics::legend("bottomright",
                   legend = c("Success", "Failure"),
                   fill = c("#228b22", "#b22222"),
                   border = NA,
                   bty = "n")

  invisible(NULL)
}

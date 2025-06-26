# ==============================================================================
# GROUP_TNA OBJECT SUPPORT
# ==============================================================================
# 
# Support functions for integrating with group_tna objects from the tna package.
# These functions allow tnaExtras to work seamlessly with tna package outputs.
#
# Functions:
# - is_group_tna(): Check if object is a group_tna
# - convert_group_tna(): Convert group_tna to tnaExtras format
# - extract_group_tna_info(): Extract metadata from group_tna objects
#
# ==============================================================================

#' Null coalescing operator
#' @param x First value
#' @param y Second value (used if x is NULL)
#' @return x if not NULL, otherwise y
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

#' Check if object is a group_tna object
#'
#' @param x Object to check
#' @return Logical indicating if object is group_tna
is_group_tna <- function(x) {
  inherits(x, "group_tna")
}

#' Extract metadata from group_tna object
#'
#' @param group_tna_obj A group_tna object
#' @return List with metadata information
extract_group_tna_info <- function(group_tna_obj) {
  info <- list(
    label = attr(group_tna_obj, "label"),
    levels = attr(group_tna_obj, "levels"),
    na_rm = attr(group_tna_obj, "na.rm"),
    n_groups = length(group_tna_obj)
  )
  
  # Extract state labels from first group (assuming consistent across groups)
  if (length(group_tna_obj) > 0 && !is.null(group_tna_obj[[1]]$labels)) {
    info$state_labels <- group_tna_obj[[1]]$labels
  }
  
  return(info)
}

#' Convert group_tna Object to tnaExtras Format
#'
#' Converts a group_tna object (list of tna models) to the format expected by tnaExtras functions.
#' Combines data from all groups and converts numeric codes back to text labels.
#'
#' @param group_tna_obj A group_tna object from the tna package
#' @return List with data, groups, group_col, and group_info
#' @export
convert_group_tna <- function(group_tna_obj) {
  if (!is_group_tna(group_tna_obj)) {
    stop("Input must be a group_tna object")
  }
  
  # Extract group names
  group_names <- names(group_tna_obj)
  if (is.null(group_names) || length(group_names) == 0) {
    stop("group_tna object must have named groups")
  }
  
  # Get labels from first group (should be same across all groups)
  first_group <- group_tna_obj[[1]]
  if (!"labels" %in% names(first_group)) {
    stop("tna objects must contain 'labels' element")
  }
  labels <- first_group$labels
  
  # Initialize lists to store all data
  all_data_frames <- list()
  all_groups <- character()
  
  # Process each group
  for (group_name in group_names) {
    tna_obj <- group_tna_obj[[group_name]]
    
    # Validate tna object structure
    if (!"data" %in% names(tna_obj)) {
      stop("tna object for group '", group_name, "' must contain 'data' element")
    }
    
    # Get numeric data
    numeric_data <- tna_obj$data
    n_sequences <- nrow(numeric_data)
    
    # Convert numeric codes to text labels
    text_data <- numeric_data  # Start with the same structure
    for (col_idx in 1:ncol(numeric_data)) {
      col_values <- numeric_data[, col_idx]
      # Convert numeric codes to text labels
      text_values <- character(length(col_values))
      valid_indices <- !is.na(col_values) & col_values >= 1 & col_values <= length(labels)
      text_values[valid_indices] <- labels[col_values[valid_indices]]
      text_values[!valid_indices] <- NA_character_
      text_data[, col_idx] <- text_values
    }
    
    # Convert to data.frame
    text_data <- as.data.frame(text_data, stringsAsFactors = FALSE)
    # Ensure column names are preserved
    colnames(text_data) <- colnames(numeric_data)
    
    # Store the data frame and group info
    all_data_frames[[group_name]] <- text_data
    all_groups <- c(all_groups, rep(group_name, n_sequences))
  }
  
  # Combine all data frames using rbind
  combined_data <- do.call(rbind, all_data_frames)
  
  # Reset row names
  rownames(combined_data) <- NULL
  
  # Add group column
  combined_data$Group <- all_groups
  
  # Create metadata
  group_info <- list(
    label = "Group_TNA",
    levels = group_names,
    original_labels = labels,
    n_groups = length(group_names)
  )
  
  return(list(
    data = combined_data,
    groups = all_groups,
    group_col = "Group",
    group_info = group_info
  ))
}

#' Create a mock group_tna object for testing
#'
#' @param groups List of group names
#' @param n_sequences Number of sequences per group
#' @param n_timepoints Number of time points
#' @param states Possible states
#' @return Mock group_tna object
create_mock_group_tna <- function(groups = c("GroupA", "GroupB", "GroupC"), 
                                  n_sequences = 5, 
                                  n_timepoints = 4,
                                  states = c("Active", "Average", "Disengaged")) {
  
  # Create mock data for each group
  group_data <- list()
  
  for (i in seq_along(groups)) {
    # Create sequence data for this group
    sequences <- matrix(
      sample(states, n_sequences * n_timepoints, replace = TRUE),
      nrow = n_sequences,
      ncol = n_timepoints
    )
    colnames(sequences) <- paste0("T", 1:n_timepoints)
    
    # Create group data structure
    group_data[[i]] <- list(
      data = as.data.frame(sequences),
      labels = states,
      n_sequences = n_sequences
    )
  }
  
  names(group_data) <- groups
  
  # Add attributes to mimic group_tna structure
  attr(group_data, "class") <- c("group_tna", "list")
  attr(group_data, "label") <- "Treatment"
  attr(group_data, "levels") <- groups
  attr(group_data, "na.rm") <- FALSE
  
  return(group_data)
}

#' Print method for mock group_tna objects
#'
#' @param x Mock group_tna object
#' @param ... Additional arguments
print.group_tna <- function(x, ...) {
  cat("Mock group_tna object\n")
  cat("Groups:", length(x), "\n")
  cat("Group names:", paste(attr(x, "levels"), collapse = ", "), "\n")
  cat("Label:", attr(x, "label"), "\n")
  
  for (i in seq_along(x)) {
    group_name <- attr(x, "levels")[i]
    if (is.list(x[[i]]) && "data" %in% names(x[[i]])) {
      n_seq <- nrow(x[[i]]$data)
      n_time <- ncol(x[[i]]$data)
      cat("  ", group_name, ":", n_seq, "sequences x", n_time, "timepoints\n")
    }
  }
} 
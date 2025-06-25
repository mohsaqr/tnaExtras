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

#' Convert group_tna object to tnaExtras-compatible format
#'
#' @param group_tna_obj A group_tna object from the tna package
#' @return List with $data (data.frame) and $group_info (metadata)
convert_group_tna <- function(group_tna_obj) {
  if (!is_group_tna(group_tna_obj)) {
    stop("Object is not a group_tna object")
  }
  
  # Extract metadata
  group_info <- extract_group_tna_info(group_tna_obj)
  
  # Try to use combine_data if available (from tna package)
  if (exists("combine_data") && is.function(get("combine_data"))) {
    tryCatch({
      combined_data <- get("combine_data")(group_tna_obj)
      
      # The combined data should have a .group column
      if (".group" %in% names(combined_data)) {
        return(list(
          data = combined_data,
          group_info = group_info,
          group_col = ".group"
        ))
      }
    }, error = function(e) {
      warning("Failed to use combine_data(), falling back to manual conversion")
    })
  }
  
  # Manual conversion if combine_data is not available
  combined_data <- convert_group_tna_manual(group_tna_obj, group_info)
  
  return(list(
    data = combined_data,
    group_info = group_info,
    group_col = ".group"
  ))
}

#' Manual conversion of group_tna to data.frame
#'
#' @param group_tna_obj A group_tna object
#' @param group_info Metadata extracted from the object
#' @return Data.frame with .group column
convert_group_tna_manual <- function(group_tna_obj, group_info) {
  # Initialize result list
  all_data <- list()
  
  # Get group names
  group_names <- group_info$levels
  if (is.null(group_names)) {
    group_names <- names(group_tna_obj)
  }
  if (is.null(group_names)) {
    group_names <- paste0("Group_", seq_along(group_tna_obj))
  }
  
  # Process each group
  for (i in seq_along(group_tna_obj)) {
    group_data <- group_tna_obj[[i]]
    group_name <- group_names[i]
    
    # Handle different possible structures
    if (is.data.frame(group_data)) {
      # If it's already a data.frame, use it directly
      group_df <- group_data
    } else if (is.matrix(group_data)) {
      # If it's a matrix, convert to data.frame
      group_df <- as.data.frame(group_data)
    } else if (is.list(group_data) && "data" %in% names(group_data)) {
      # If it's a list with a data component
      group_df <- as.data.frame(group_data$data)
    } else if (is.list(group_data) && length(group_data) > 0) {
      # Try to convert list to data.frame
      tryCatch({
        group_df <- as.data.frame(group_data)
      }, error = function(e) {
        # If conversion fails, try to extract sequence-like data
        sequence_cols <- sapply(group_data, function(x) is.vector(x) && length(x) > 1)
        if (any(sequence_cols)) {
          group_df <- as.data.frame(group_data[sequence_cols])
        } else {
          warning("Could not convert group ", i, " to data.frame")
          group_df <- data.frame()
        }
      })
    } else {
      warning("Unknown group_tna structure for group ", i)
      group_df <- data.frame()
    }
    
    # Add group column
    if (nrow(group_df) > 0) {
      group_df$.group <- group_name
      all_data[[i]] <- group_df
    }
  }
  
  # Combine all groups
  if (length(all_data) > 0) {
    combined_data <- do.call(rbind, all_data)
    rownames(combined_data) <- NULL
  } else {
    combined_data <- data.frame(.group = character(0))
  }
  
  return(combined_data)
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
# ==============================================================================
# ASSOCIATION RULES UTILITIES
# ==============================================================================

#' Prepare Transaction Data for Association Rule Mining
#'
#' Converts various input formats to standardized transaction format.
#' Supports list, matrix, and data.frame inputs.
#'
#' @param data Transaction data in various formats
#' @param item_cols Column names/indices for items (for data.frame input)
#' @param transaction_col Column name/index for transaction ID (for data.frame input)
#' @return List with standardized transactions and metadata
#' @export
#' @examples
#' # List format - learning activities
#' learning_sessions <- list(c("plan", "discuss"), c("discuss", "execute"), c("plan", "execute"))
#' prepared <- prepare_transactions(learning_sessions)
#'
#' # Data frame format - student actions
#' student_actions <- data.frame(
#'   student_id = c(1, 1, 2, 2, 3, 3),
#'   action = c("plan", "discuss", "discuss", "execute", "plan", "execute")
#' )
#' prepared <- prepare_transactions(student_actions, transaction_col = "student_id", item_cols = "action")
prepare_transactions <- function(data, item_cols = NULL, transaction_col = NULL) {
  if (is.list(data) && !is.data.frame(data)) {
    # List format: each element is a transaction
    transactions <- lapply(data, function(x) {
      if (is.character(x) || is.factor(x)) {
        as.character(unique(x))
      } else {
        stop("List elements must be character or factor vectors")
      }
    })
  } else if (is.matrix(data)) {
    # Matrix format: binary or transaction matrix
    if (all(data %in% c(0, 1, TRUE, FALSE))) {
      # Binary matrix
      transactions <- apply(data, 1, function(row) {
        colnames(data)[which(row > 0)]
      })
      transactions <- transactions[sapply(transactions, length) > 0]
    } else {
      stop("Matrix must be binary (0/1 or TRUE/FALSE)")
    }
  } else if (is.data.frame(data)) {
    # Data frame format
    if (is.null(transaction_col) && is.null(item_cols)) {
      # Auto-detect format
      if (ncol(data) == 2) {
        # Assume first column is transaction ID, second is item (long format)
        transaction_col <- 1
        item_cols <- 2
      } else {
        # Check if all columns are character/factor (wide format)
        char_cols <- sapply(data, function(x) is.character(x) || is.factor(x))

        if (all(char_cols)) {
          # Wide format: each row is a transaction, each column contains items
          # This is the format used by tna::group_regulation
          transactions <- apply(data, 1, function(row) {
            items <- as.character(row)
            items <- items[!is.na(items)] # Remove NA values
            unique(items) # Get unique items
          })
          # Filter out empty transactions
          transactions <- transactions[sapply(transactions, length) > 0]

          # Skip the long-format processing below
          transaction_col <- NULL
          item_cols <- NULL
        } else if (sum(char_cols) >= 2) {
          # Long format: find transaction ID and item columns
          transaction_col <- which(char_cols)[1]
          item_cols <- which(char_cols)[2]
        } else {
          stop("Unable to auto-detect transaction and item columns. Please specify explicitly.")
        }
      }
    }

    # Process long format if columns were specified or detected
    if (!is.null(transaction_col) && !is.null(item_cols)) {
      # Convert to transactions (long format)
      trans_df <- data[, c(transaction_col, item_cols), drop = FALSE]
      colnames(trans_df) <- c("transaction", "item")

      # Group by transaction
      transaction_list <- split(as.character(trans_df$item), trans_df$transaction)
      transactions <- lapply(transaction_list, unique)
    }
  } else {
    stop("Data must be a list, matrix, or data.frame")
  }

  # Validate and clean transactions
  transactions <- transactions[sapply(transactions, length) > 0]

  if (length(transactions) == 0) {
    stop("No valid transactions found")
  }

  # Get all unique items
  all_items <- sort(unique(unlist(transactions)))

  return(list(
    transactions = transactions,
    n_transactions = length(transactions),
    n_items = length(all_items),
    items = all_items,
    avg_transaction_length = mean(sapply(transactions, length))
  ))
}

#' Create Binary Transaction Matrix
#'
#' Converts transaction list to binary matrix format for efficient processing.
#'
#' @param transactions List of transactions
#' @param items Optional vector of item names (uses all unique items if NULL)
#' @return Binary matrix where rows are transactions and columns are items
#' @export
create_transaction_matrix <- function(transactions, items = NULL) {
  if (is.null(items)) {
    items <- sort(unique(unlist(transactions)))
  }

  n_trans <- length(transactions)
  n_items <- length(items)

  # Create binary matrix
  trans_matrix <- matrix(FALSE, nrow = n_trans, ncol = n_items)
  colnames(trans_matrix) <- items
  rownames(trans_matrix) <- paste0("T", seq_len(n_trans))

  # Fill matrix
  for (i in seq_len(n_trans)) {
    trans_items <- transactions[[i]]
    valid_items <- intersect(trans_items, items)
    if (length(valid_items) > 0) {
      trans_matrix[i, valid_items] <- TRUE
    }
  }

  return(trans_matrix)
}

#' Validate Association Rule Parameters
#'
#' Checks that support, confidence, and lift parameters are valid.
#'
#' @param min_support Minimum support threshold
#' @param min_confidence Minimum confidence threshold
#' @param min_lift Minimum lift threshold
#' @export
validate_parameters <- function(min_support, min_confidence, min_lift) {
  if (!is.numeric(min_support) || length(min_support) != 1 ||
    min_support < 0 || min_support > 1) {
    stop("min_support must be a single number between 0 and 1")
  }

  if (!is.numeric(min_confidence) || length(min_confidence) != 1 ||
    min_confidence < 0 || min_confidence > 1) {
    stop("min_confidence must be a single number between 0 and 1")
  }

  if (!is.numeric(min_lift) || length(min_lift) != 1 || min_lift < 0) {
    stop("min_lift must be a single non-negative number")
  }

  if (min_support == 0) {
    warning("min_support = 0 may result in very large number of patterns")
  }

  if (min_confidence == 0) {
    warning("min_confidence = 0 may result in meaningless rules")
  }
}

#' Calculate Association Rule Metrics
#'
#' Computes support, confidence, lift, and conviction for an association rule.
#'
#' @param antecedent Vector of antecedent items
#' @param consequent Vector of consequent items
#' @param trans_matrix Binary transaction matrix
#' @return List of calculated metrics
#' @export
calculate_rule_metrics <- function(antecedent, consequent, trans_matrix) {
  # Count transactions containing antecedent
  antecedent_mask <- apply(trans_matrix[, antecedent, drop = FALSE], 1, all)
  antecedent_count <- sum(antecedent_mask)

  # Count transactions containing consequent
  consequent_mask <- apply(trans_matrix[, consequent, drop = FALSE], 1, all)
  consequent_count <- sum(consequent_mask)

  # Count transactions containing both antecedent and consequent
  both_mask <- antecedent_mask & consequent_mask
  both_count <- sum(both_mask)

  n_transactions <- nrow(trans_matrix)

  # Calculate metrics
  support <- both_count / n_transactions
  confidence <- if (antecedent_count > 0) both_count / antecedent_count else 0

  # Lift = P(A ∩ B) / (P(A) * P(B))
  lift <- if (antecedent_count > 0 && consequent_count > 0) {
    (both_count / n_transactions) / ((antecedent_count / n_transactions) * (consequent_count / n_transactions))
  } else {
    0
  }

  # Conviction = (1 - P(B)) / (1 - confidence)
  conviction <- if (confidence < 1 && consequent_count < n_transactions) {
    (1 - consequent_count / n_transactions) / (1 - confidence)
  } else {
    Inf
  }

  return(list(
    support = support,
    confidence = confidence,
    lift = lift,
    conviction = conviction,
    count = both_count,
    antecedent_count = antecedent_count,
    consequent_count = consequent_count
  ))
}

#' Calculate Support for Itemset
#'
#' Computes the support (frequency) of an itemset in the transaction data.
#'
#' @param itemset Vector of items
#' @param trans_matrix Binary transaction matrix
#' @return Support value (proportion of transactions containing the itemset)
#' @export
calculate_itemset_support <- function(itemset, trans_matrix) {
  if (length(itemset) == 0) {
    return(0)
  }

  if (length(itemset) == 1) {
    return(sum(trans_matrix[, itemset]) / nrow(trans_matrix))
  }

  # For multiple items, check transactions containing all items
  itemset_mask <- apply(trans_matrix[, itemset, drop = FALSE], 1, all)
  return(sum(itemset_mask) / nrow(trans_matrix))
}

#' Filter Association Rules
#'
#' Filters rules based on quality metrics and other criteria.
#'
#' @param rules Data frame of association rules
#' @param min_support Minimum support threshold
#' @param min_confidence Minimum confidence threshold
#' @param min_lift Minimum lift threshold
#' @param max_length Maximum rule length (antecedent + consequent)
#' @param min_count Minimum absolute count
#' @return Filtered rules data frame
#' @export
filter_association_rules <- function(rules, min_support = 0, min_confidence = 0,
                                     min_lift = 0, max_length = Inf, min_count = 0) {
  if (nrow(rules) == 0) {
    return(rules)
  }

  # Apply filters
  mask <- rep(TRUE, nrow(rules))

  if (min_support > 0) {
    mask <- mask & (rules$support >= min_support)
  }

  if (min_confidence > 0) {
    mask <- mask & (rules$confidence >= min_confidence)
  }

  if (min_lift > 0) {
    mask <- mask & (rules$lift >= min_lift)
  }

  if (min_count > 0) {
    mask <- mask & (rules$count >= min_count)
  }

  if (is.finite(max_length)) {
    antecedent_lengths <- sapply(strsplit(rules$antecedent, ","), length)
    consequent_lengths <- sapply(strsplit(rules$consequent, ","), length)
    total_lengths <- antecedent_lengths + consequent_lengths
    mask <- mask & (total_lengths <= max_length)
  }

  filtered_rules <- rules[mask, , drop = FALSE]
  rownames(filtered_rules) <- NULL

  return(filtered_rules)
}

#' Rank Association Rules
#'
#' Sorts association rules by specified metric.
#'
#' @param rules Data frame of association rules
#' @param by Metric to sort by ("support", "confidence", "lift", "conviction")
#' @param decreasing Whether to sort in decreasing order
#' @return Sorted rules data frame
#' @export
rank_association_rules <- function(rules, by = "lift", decreasing = TRUE) {
  if (nrow(rules) == 0) {
    return(rules)
  }

  if (!by %in% names(rules)) {
    stop("Metric '", by, "' not found in rules. Available: ", paste(names(rules), collapse = ", "))
  }

  sorted_rules <- rules[order(rules[[by]], decreasing = decreasing), , drop = FALSE]
  rownames(sorted_rules) <- NULL

  return(sorted_rules)
}

#' Find Redundant Rules
#'
#' Identifies rules that are redundant based on more general rules with higher confidence.
#'
#' @param rules Data frame of association rules
#' @param confidence_threshold Minimum confidence difference to consider redundancy
#' @return Logical vector indicating which rules are redundant
#' @export
find_redundant_rules <- function(rules, confidence_threshold = 0.01) {
  if (nrow(rules) == 0) {
    return(logical(0))
  }

  redundant <- rep(FALSE, nrow(rules))

  for (i in 1:nrow(rules)) {
    rule_i <- rules[i, ]
    antecedent_i <- strsplit(rule_i$antecedent, ",")[[1]]
    consequent_i <- strsplit(rule_i$consequent, ",")[[1]]

    for (j in 1:nrow(rules)) {
      if (i != j && !redundant[i]) {
        rule_j <- rules[j, ]
        antecedent_j <- strsplit(rule_j$antecedent, ",")[[1]]
        consequent_j <- strsplit(rule_j$consequent, ",")[[1]]

        # Check if rule j is more general than rule i (subset of antecedent, same consequent)
        if (setequal(consequent_i, consequent_j) &&
          all(antecedent_j %in% antecedent_i) &&
          length(antecedent_j) < length(antecedent_i) &&
          rule_j$confidence >= rule_i$confidence - confidence_threshold) {
          redundant[i] <- TRUE
          break
        }
      }
    }
  }

  return(redundant)
}

#' Extract Rules by Item
#'
#' Extracts rules containing specific items in antecedent or consequent.
#'
#' @param rules Data frame of association rules
#' @param items Vector of items to search for
#' @param side Where to search ("both", "antecedent", "consequent")
#' @param match_type How to match ("any", "all", "exact")
#' @return Filtered rules containing the specified items
#' @export
extract_rules_by_item <- function(rules, items, side = "both", match_type = "any") {
  if (nrow(rules) == 0) {
    return(rules)
  }

  if (!side %in% c("both", "antecedent", "consequent")) {
    stop("side must be 'both', 'antecedent', or 'consequent'")
  }

  if (!match_type %in% c("any", "all", "exact")) {
    stop("match_type must be 'any', 'all', or 'exact'")
  }

  mask <- rep(FALSE, nrow(rules))

  for (i in 1:nrow(rules)) {
    rule <- rules[i, ]
    antecedent_items <- strsplit(rule$antecedent, ",")[[1]]
    consequent_items <- strsplit(rule$consequent, ",")[[1]]

    # Check antecedent
    antecedent_match <- FALSE
    if (side %in% c("both", "antecedent")) {
      if (match_type == "any") {
        antecedent_match <- any(items %in% antecedent_items)
      } else if (match_type == "all") {
        antecedent_match <- all(items %in% antecedent_items)
      } else if (match_type == "exact") {
        antecedent_match <- setequal(items, antecedent_items)
      }
    }

    # Check consequent
    consequent_match <- FALSE
    if (side %in% c("both", "consequent")) {
      if (match_type == "any") {
        consequent_match <- any(items %in% consequent_items)
      } else if (match_type == "all") {
        consequent_match <- all(items %in% consequent_items)
      } else if (match_type == "exact") {
        consequent_match <- setequal(items, consequent_items)
      }
    }

    # Combine matches based on side
    if (side == "both") {
      mask[i] <- antecedent_match || consequent_match
    } else if (side == "antecedent") {
      mask[i] <- antecedent_match
    } else if (side == "consequent") {
      mask[i] <- consequent_match
    }
  }

  filtered_rules <- rules[mask, , drop = FALSE]
  rownames(filtered_rules) <- NULL

  return(filtered_rules)
}

#' Convert Rules to Transactions
#'
#' Converts association rules back to transaction format for further analysis.
#'
#' @param rules Data frame of association rules
#' @param include_metrics Whether to include rule metrics as items
#' @return List of transactions representing the rules
#' @export
rules_to_transactions <- function(rules, include_metrics = FALSE) {
  if (nrow(rules) == 0) {
    return(list())
  }

  transactions <- list()

  for (i in 1:nrow(rules)) {
    rule <- rules[i, ]
    antecedent_items <- strsplit(rule$antecedent, ",")[[1]]
    consequent_items <- strsplit(rule$consequent, ",")[[1]]

    # Create transaction from rule items
    rule_items <- c(antecedent_items, consequent_items)

    if (include_metrics) {
      # Add metrics as categorical items
      support_cat <- paste0("support_", cut(rule$support,
        breaks = c(0, 0.1, 0.3, 0.5, 1),
        labels = c("low", "medium", "high", "very_high")
      ))
      confidence_cat <- paste0("confidence_", cut(rule$confidence,
        breaks = c(0, 0.5, 0.8, 0.95, 1),
        labels = c("low", "medium", "high", "very_high")
      ))
      lift_cat <- paste0("lift_", cut(rule$lift,
        breaks = c(0, 1, 2, 5, Inf),
        labels = c("negative", "neutral", "positive", "strong")
      ))

      rule_items <- c(rule_items, support_cat, confidence_cat, lift_cat)
    }

    transactions[[i]] <- unique(rule_items)
  }

  return(transactions)
}

#' Calculate Rule Overlap
#'
#' Calculates the overlap between rules based on shared items.
#'
#' @param rules Data frame of association rules
#' @param method Overlap method ("jaccard", "overlap", "dice")
#' @return Matrix of pairwise rule similarities
#' @export
calculate_rule_overlap <- function(rules, method = "jaccard") {
  if (nrow(rules) == 0) {
    return(matrix(numeric(0), nrow = 0, ncol = 0))
  }

  if (!method %in% c("jaccard", "overlap", "dice")) {
    stop("method must be 'jaccard', 'overlap', or 'dice'")
  }

  n_rules <- nrow(rules)
  similarity_matrix <- matrix(0, nrow = n_rules, ncol = n_rules)

  # Extract all items for each rule
  rule_items <- list()
  for (i in 1:n_rules) {
    antecedent <- strsplit(rules$antecedent[i], ",")[[1]]
    consequent <- strsplit(rules$consequent[i], ",")[[1]]
    rule_items[[i]] <- unique(c(antecedent, consequent))
  }

  # Calculate pairwise similarities
  for (i in 1:n_rules) {
    for (j in i:n_rules) {
      items_i <- rule_items[[i]]
      items_j <- rule_items[[j]]

      intersection <- length(intersect(items_i, items_j))
      union_size <- length(union(items_i, items_j))

      if (method == "jaccard") {
        similarity <- if (union_size > 0) intersection / union_size else 0
      } else if (method == "overlap") {
        similarity <- if (min(length(items_i), length(items_j)) > 0) {
          intersection / min(length(items_i), length(items_j))
        } else {
          0
        }
      } else if (method == "dice") {
        similarity <- if (length(items_i) + length(items_j) > 0) {
          2 * intersection / (length(items_i) + length(items_j))
        } else {
          0
        }
      }

      similarity_matrix[i, j] <- similarity_matrix[j, i] <- similarity
    }
  }

  rownames(similarity_matrix) <- colnames(similarity_matrix) <- paste0("Rule_", 1:n_rules)

  return(similarity_matrix)
}

#' Export Association Rules
#'
#' Exports association rules to various formats.
#'
#' @param rules Association rules object or data frame
#' @param file Output file path
#' @param format Output format ("csv", "json", "txt")
#' @param include_summary Whether to include summary information
#' @export
export_association_rules <- function(rules, file, format = "csv", include_summary = TRUE) {
  # Handle association_rules object
  if (inherits(rules, "association_rules")) {
    rules_df <- rules$rules
    summary_info <- rules$summary
    parameters <- rules$parameters
  } else if (is.data.frame(rules)) {
    rules_df <- rules
    summary_info <- NULL
    parameters <- NULL
  } else {
    stop("rules must be an association_rules object or data frame")
  }

  if (nrow(rules_df) == 0) {
    warning("No rules to export")
    return(invisible(NULL))
  }

  # Export based on format
  if (format == "csv") {
    write.csv(rules_df, file, row.names = FALSE)

    if (include_summary && !is.null(summary_info)) {
      # Write summary to separate file
      summary_file <- sub("\\.csv$", "_summary.txt", file)
      cat("Association Rules Summary\n", file = summary_file)
      cat("========================\n\n", file = summary_file, append = TRUE)

      if (!is.null(parameters)) {
        cat("Parameters:\n", file = summary_file, append = TRUE)
        for (param in names(parameters)) {
          cat(paste0("  ", param, ": ", parameters[[param]], "\n"),
            file = summary_file, append = TRUE
          )
        }
        cat("\n", file = summary_file, append = TRUE)
      }

      cat("Summary:\n", file = summary_file, append = TRUE)
      for (stat in names(summary_info)) {
        cat(paste0("  ", stat, ": ", summary_info[[stat]], "\n"),
          file = summary_file, append = TRUE
        )
      }
    }
  } else if (format == "json") {
    # Simple JSON export (no external dependencies)
    json_content <- paste0('{"rules": [')

    rule_jsons <- character(nrow(rules_df))
    for (i in 1:nrow(rules_df)) {
      rule <- rules_df[i, ]
      rule_json <- paste0(
        '{"antecedent": "', rule$antecedent, '", ',
        '"consequent": "', rule$consequent, '", ',
        '"support": ', rule$support, ", ",
        '"confidence": ', rule$confidence, ", ",
        '"lift": ', rule$lift, ", ",
        '"count": ', rule$count, "}"
      )
      rule_jsons[i] <- rule_json
    }

    json_content <- paste0(json_content, paste(rule_jsons, collapse = ", "), "]")

    if (include_summary && !is.null(summary_info)) {
      json_content <- paste0(json_content, ', "summary": {')
      summary_parts <- paste0('"', names(summary_info), '": ', summary_info, collapse = ", ")
      json_content <- paste0(json_content, summary_parts, "}")
    }

    json_content <- paste0(json_content, "}")

    writeLines(json_content, file)
  } else if (format == "txt") {
    # Human-readable text format
    cat("Association Rules\n", file = file)
    cat("================\n\n", file = file, append = TRUE)

    if (include_summary && !is.null(summary_info)) {
      cat("Summary:\n", file = file, append = TRUE)
      for (stat in names(summary_info)) {
        cat(paste0("  ", stat, ": ", summary_info[[stat]], "\n"),
          file = file, append = TRUE
        )
      }
      cat("\n", file = file, append = TRUE)
    }

    cat("Rules:\n", file = file, append = TRUE)
    cat("------\n", file = file, append = TRUE)

    for (i in 1:nrow(rules_df)) {
      rule <- rules_df[i, ]
      cat(sprintf("%d. %s => %s\n", i, rule$antecedent, rule$consequent),
        file = file, append = TRUE
      )
      cat(
        sprintf(
          "   Support: %.3f, Confidence: %.3f, Lift: %.3f\n\n",
          rule$support, rule$confidence, rule$lift
        ),
        file = file, append = TRUE
      )
    }
  } else {
    stop("format must be 'csv', 'json', or 'txt'")
  }

  cat("Rules exported to:", file, "\n")
  return(invisible(NULL))
}

#' Get Rules Network Data
#'
#' Extracts network data (nodes and edges) from association rules with flexible
#' aggregation options for handling duplicate edges.
#'
#' @param rules Association rules object or data frame
#' @param top_n Number of top rules to include (default: Inf for all)
#' @param weight_by Metric to use for edge weights. Options:
#'   \itemize{
#'     \item \code{"lift"} - Use lift values (default)
#'     \item \code{"confidence"} - Use confidence values
#'     \item \code{"support"} - Use support values
#'   }
#' @param agg_fun Aggregation function for duplicate edges. Options:
#'   \itemize{
#'     \item \code{"mean"} - Average values (default)
#'     \item \code{"max"} - Maximum value
#'     \item \code{"min"} - Minimum value
#'     \item \code{"sum"} - Sum of values
#'     \item \code{"median"} - Median value
#'   }
#'
#' @return List containing 'nodes' (data frame) and 'edges' (data frame)
#'
#' @details
#' When multiple rules produce the same antecedent→consequent item pair,
#' the edges are aggregated using the specified aggregation function.
#' The weight is determined by the \code{weight_by} parameter, while
#' lift, confidence, and support are all aggregated using \code{agg_fun}.
#'
#' @examples
#' \dontrun{
#' # Get network with lift-based weights, averaged
#' network <- get_rules_network(rules, top_n = 20)
#'
#' # Use confidence as weights, take maximum
#' network <- get_rules_network(rules,
#'   top_n = 20,
#'   weight_by = "confidence",
#'   agg_fun = "max"
#' )
#'
#' # Use support as weights, sum duplicates
#' network <- get_rules_network(rules,
#'   top_n = 20,
#'   weight_by = "support",
#'   agg_fun = "sum"
#' )
#' }
#'
#' @export
get_rules_network <- function(rules, top_n = Inf,
                              weight_by = c("lift", "confidence", "support"),
                              agg_fun = c("mean", "max", "min", "sum", "median")) {
  # Match arguments
  weight_by <- match.arg(weight_by)
  agg_fun <- match.arg(agg_fun)

  # Handle association_rules object
  if (inherits(rules, "association_rules")) {
    rules_df <- rules$rules
  } else if (is.data.frame(rules)) {
    rules_df <- rules
  } else {
    stop("rules must be an association_rules object or data frame")
  }

  if (nrow(rules_df) == 0) {
    return(list(
      nodes = data.frame(id = character(0), label = character(0), stringsAsFactors = FALSE),
      edges = data.frame(
        from = character(0), to = character(0),
        weight = numeric(0), lift = numeric(0),
        confidence = numeric(0), support = numeric(0),
        stringsAsFactors = FALSE
      )
    ))
  }

  # Select top rules
  if (is.finite(top_n) && top_n < nrow(rules_df)) {
    rules_df <- head(rules_df[order(rules_df$lift, decreasing = TRUE), ], top_n)
  }

  # Get aggregation function
  agg_func <- switch(agg_fun,
    "mean" = mean,
    "max" = max,
    "min" = min,
    "sum" = sum,
    "median" = median
  )

  # Initialize edges data frame
  # Pre-allocate lists for better performance
  edge_list <- vector("list", nrow(rules_df))

  for (i in seq_len(nrow(rules_df))) {
    antecedent_items <- strsplit(rules_df$antecedent[i], ",")[[1]]
    consequent_items <- strsplit(rules_df$consequent[i], ",")[[1]]

    # Create edges from each antecedent to each consequent
    # Expand grid of all combinations
    combinations <- expand.grid(from = antecedent_items, to = consequent_items, stringsAsFactors = FALSE)

    # Add metrics
    combinations$lift <- rules_df$lift[i]
    combinations$confidence <- rules_df$confidence[i]
    combinations$support <- rules_df$support[i]

    edge_list[[i]] <- combinations
  }

  # Combine all edges
  edges <- do.call(rbind, edge_list)

  # Aggregate duplicate edges (same from-to pair)
  # Group by from-to and aggregate metrics
  edge_key <- paste(edges$from, edges$to, sep = "->")
  unique_keys <- unique(edge_key)

  aggregated_edges <- do.call(rbind, lapply(unique_keys, function(key) {
    idx <- which(edge_key == key)
    if (length(idx) == 1) {
      # Single edge: add weight based on weight_by parameter
      result <- edges[idx, ]
      result$weight <- result[[weight_by]]
      result
    } else {
      # Multiple edges with same from-to: aggregate
      data.frame(
        from = edges$from[idx[1]],
        to = edges$to[idx[1]],
        weight = agg_func(edges[[weight_by]][idx]),
        lift = agg_func(edges$lift[idx]),
        confidence = agg_func(edges$confidence[idx]),
        support = agg_func(edges$support[idx]),
        stringsAsFactors = FALSE
      )
    }
  }))

  # Extract unique nodes
  all_items <- unique(c(aggregated_edges$from, aggregated_edges$to))
  nodes <- data.frame(
    id = all_items,
    label = all_items,
    stringsAsFactors = FALSE
  )

  return(list(nodes = nodes, edges = aggregated_edges))
}

# ==============================================================================
# ASSOCIATION RULES ANALYSIS
# ==============================================================================

#' Association Rules Analysis using Apriori Algorithm
#'
#' Discovers association rules from transaction data using the classic Apriori algorithm.
#' This algorithm uses a level-wise approach to find frequent itemsets and generate rules.
#'
#' @param transactions Transaction data in various formats (list, matrix, or data.frame)
#' @param min_support Minimum support threshold (default: 0.1)
#' @param min_confidence Minimum confidence threshold (default: 0.8)
#' @param min_lift Minimum lift threshold (default: 1.0)
#' @param max_length Maximum length of itemsets to consider (default: 10)
#' @param verbose Whether to show progress messages (default: TRUE)
#' @return Object of class "association_rules" containing discovered rules
#' @export
#' @examples
#' # Example with learning activity data
#' transactions <- list(
#'   c("plan", "discuss", "execute", "reflect"),
#'   c("plan", "research", "analyze"),
#'   c("discuss", "execute", "collaborate", "reflect"),
#'   c("plan", "discuss", "execute", "evaluate"),
#'   c("research", "analyze", "collaborate")
#' )
#' 
#' rules <- apriori_rules(transactions, min_support = 0.2, min_confidence = 0.6)
#' print(rules)
#' summary(rules)
apriori_rules <- function(transactions, min_support = 0.1, min_confidence = 0.8, 
                         min_lift = 1.0, max_length = 10, verbose = TRUE) {
  
  # Input validation
  if (missing(transactions)) {
    stop("transactions parameter is required")
  }
  
  validate_parameters(min_support, min_confidence, min_lift)
  
  if (verbose) cat("Preparing transaction data...\n")
  
  # Prepare transactions
  trans_data <- prepare_transactions(transactions)
  trans_matrix <- create_transaction_matrix(trans_data$transactions)
  
  if (verbose) {
    cat("Transaction summary:\n")
    cat("  - Number of transactions:", nrow(trans_matrix), "\n")
    cat("  - Number of unique items:", ncol(trans_matrix), "\n")
    cat("  - Average transaction length:", round(mean(rowSums(trans_matrix)), 2), "\n")
  }
  
  if (verbose) cat("Finding frequent itemsets using Apriori...\n")
  
  # Find frequent itemsets
  frequent_itemsets <- find_frequent_itemsets_apriori(trans_matrix, min_support, max_length, verbose)
  
  if (length(frequent_itemsets) == 0) {
    warning("No frequent itemsets found with given support threshold")
    return(create_empty_rules_object("apriori"))
  }
  
  if (verbose) cat("Generating association rules...\n")
  
  # Generate rules
  rules <- generate_association_rules(frequent_itemsets, trans_matrix, min_confidence, min_lift)
  
  # Create result object
  result <- list(
    rules = rules,
    algorithm = "apriori",
    parameters = list(
      min_support = min_support,
      min_confidence = min_confidence,
      min_lift = min_lift,
      max_length = max_length
    ),
    summary = list(
      n_transactions = nrow(trans_matrix),
      n_items = ncol(trans_matrix),
      n_frequent_itemsets = length(frequent_itemsets),
      n_rules = nrow(rules)
    ),
    item_names = colnames(trans_matrix)
  )
  
  class(result) <- "association_rules"
  
  if (verbose) {
    cat("Association rule mining completed!\n")
    cat("  - Frequent itemsets found:", length(frequent_itemsets), "\n")
    cat("  - Rules generated:", nrow(rules), "\n")
  }
  
  return(result)
}

#' Association Rules Analysis using FP-Growth Algorithm
#'
#' Discovers association rules from transaction data using the FP-Growth algorithm.
#' This algorithm builds an FP-tree for efficient frequent pattern mining.
#'
#' @param transactions Transaction data in various formats (list, matrix, or data.frame)
#' @param min_support Minimum support threshold (default: 0.1)
#' @param min_confidence Minimum confidence threshold (default: 0.8)
#' @param min_lift Minimum lift threshold (default: 1.0)
#' @param verbose Whether to show progress messages (default: TRUE)
#' @return Object of class "association_rules" containing discovered rules
#' @export
#' @examples
#' # Example with learning activity data
#' transactions <- list(
#'   c("plan", "discuss", "execute", "reflect"),
#'   c("plan", "research", "analyze"),
#'   c("discuss", "execute", "collaborate", "reflect"),
#'   c("plan", "discuss", "execute", "evaluate"),
#'   c("research", "analyze", "collaborate")
#' )
#' 
#' rules <- fp_growth_rules(transactions, min_support = 0.2, min_confidence = 0.6)
#' print(rules)
#' summary(rules)
fp_growth_rules <- function(transactions, min_support = 0.1, min_confidence = 0.8, 
                           min_lift = 1.0, verbose = TRUE) {
  
  # Input validation
  if (missing(transactions)) {
    stop("transactions parameter is required")
  }
  
  validate_parameters(min_support, min_confidence, min_lift)
  
  if (verbose) cat("Preparing transaction data...\n")
  
  # Prepare transactions
  trans_data <- prepare_transactions(transactions)
  trans_matrix <- create_transaction_matrix(trans_data$transactions)
  
  if (verbose) {
    cat("Transaction summary:\n")
    cat("  - Number of transactions:", nrow(trans_matrix), "\n")
    cat("  - Number of unique items:", ncol(trans_matrix), "\n")
    cat("  - Average transaction length:", round(mean(rowSums(trans_matrix)), 2), "\n")
  }
  
  if (verbose) cat("Building FP-Tree...\n")
  
  # Build FP-tree and mine patterns
  fp_tree <- build_fp_tree(trans_data$transactions, min_support, verbose)
  
  if (verbose) cat("Mining frequent patterns...\n")
  
  frequent_patterns <- mine_fp_tree(fp_tree, min_support, verbose)
  
  if (length(frequent_patterns) == 0) {
    warning("No frequent patterns found with given support threshold")
    return(create_empty_rules_object("fp_growth"))
  }
  
  if (verbose) cat("Generating association rules...\n")
  
  # Generate rules from patterns
  rules <- generate_rules_from_patterns(frequent_patterns, trans_matrix, min_confidence, min_lift)
  
  # Create result object
  result <- list(
    rules = rules,
    algorithm = "fp_growth",
    parameters = list(
      min_support = min_support,
      min_confidence = min_confidence,
      min_lift = min_lift
    ),
    summary = list(
      n_transactions = nrow(trans_matrix),
      n_items = ncol(trans_matrix),
      n_frequent_patterns = length(frequent_patterns),
      n_rules = nrow(rules)
    ),
    item_names = colnames(trans_matrix)
  )
  
  class(result) <- "association_rules"
  
  if (verbose) {
    cat("Association rule mining completed!\n")
    cat("  - Frequent patterns found:", length(frequent_patterns), "\n")
    cat("  - Rules generated:", nrow(rules), "\n")
  }
  
  return(result)
}

#' Compare Association Rule Mining Algorithms
#'
#' Compares results from different association rule mining algorithms.
#'
#' @param rules_list List of association rules objects from different algorithms
#' @param metrics Vector of metrics to compare (default: c("support", "confidence", "lift"))
#' @return Comparison summary
#' @export
compare_rule_algorithms <- function(rules_list, metrics = c("support", "confidence", "lift")) {
  
  if (!is.list(rules_list) || length(rules_list) < 2) {
    stop("rules_list must be a list with at least 2 association rules objects")
  }
  
  # Validate all objects are association_rules
  algorithms <- sapply(rules_list, function(x) {
    if (!inherits(x, "association_rules")) {
      stop("All objects must be of class 'association_rules'")
    }
    x$algorithm
  })
  
  if (is.null(names(rules_list))) {
    names(rules_list) <- algorithms
  }
  
  cat("Association Rule Algorithm Comparison\n")
  cat("=====================================\n\n")
  
  # Summary comparison
  summary_df <- data.frame(
    Algorithm = algorithms,
    N_Rules = sapply(rules_list, function(x) x$summary$n_rules),
    N_Frequent_Items = sapply(rules_list, function(x) {
      if ("n_frequent_itemsets" %in% names(x$summary)) {
        x$summary$n_frequent_itemsets
      } else {
        x$summary$n_frequent_patterns
      }
    }),
    Min_Support = sapply(rules_list, function(x) x$parameters$min_support),
    Min_Confidence = sapply(rules_list, function(x) x$parameters$min_confidence),
    Min_Lift = sapply(rules_list, function(x) x$parameters$min_lift),
    stringsAsFactors = FALSE
  )
  
  print(summary_df)
  cat("\n")
  
  # Metric comparison for rules
  for (metric in metrics) {
    cat(paste("Top 10 Rules by", toupper(metric), ":\n"))
    cat(paste(rep("-", 30), collapse = ""), "\n")
    
    for (name in names(rules_list)) {
      rules_obj <- rules_list[[name]]
      if (nrow(rules_obj$rules) > 0 && metric %in% names(rules_obj$rules)) {
        top_rules <- head(rules_obj$rules[order(rules_obj$rules[[metric]], decreasing = TRUE), ], 5)
        cat(paste(name, "algorithm:\n"))
        for (i in 1:nrow(top_rules)) {
          cat(sprintf("  %s => %s (%.3f)\n", 
                     top_rules$antecedent[i], 
                     top_rules$consequent[i], 
                     top_rules[[metric]][i]))
        }
        cat("\n")
      }
    }
  }
  
  return(invisible(summary_df))
}

# ==============================================================================
# APRIORI ALGORITHM CORE FUNCTIONS
# ==============================================================================

#' Find Frequent Itemsets using Apriori Algorithm
#'
#' @param trans_matrix Binary transaction matrix
#' @param min_support Minimum support threshold
#' @param max_length Maximum itemset length
#' @param verbose Whether to show progress
#' @return List of frequent itemsets
find_frequent_itemsets_apriori <- function(trans_matrix, min_support, max_length, verbose = TRUE) {
  
  n_transactions <- nrow(trans_matrix)
  min_count <- ceiling(min_support * n_transactions)
  
  # Find frequent 1-itemsets
  item_counts <- colSums(trans_matrix)
  frequent_1 <- names(item_counts[item_counts >= min_count])
  
  if (length(frequent_1) == 0) {
    return(list())
  }
  
  if (verbose) cat("  - Level 1: Found", length(frequent_1), "frequent items\n")
  
  # Initialize result with 1-itemsets
  all_frequent <- list()
  for (item in frequent_1) {
    all_frequent[[item]] <- list(
      itemset = item,
      support = item_counts[item] / n_transactions,
      count = item_counts[item]
    )
  }
  
  current_frequent <- frequent_1
  level <- 2
  
  # Generate higher-level itemsets
  while (length(current_frequent) > 1 && level <= max_length) {
    
    # Generate candidates
    candidates <- generate_candidates_apriori(current_frequent, level)
    
    if (length(candidates) == 0) break
    
    # Count candidate support
    frequent_candidates <- character(0)
    
    for (candidate in candidates) {
      items <- strsplit(candidate, ",")[[1]]
      
      # Count transactions containing all items
      count <- sum(apply(trans_matrix[, items, drop = FALSE], 1, all))
      
      if (count >= min_count) {
        frequent_candidates <- c(frequent_candidates, candidate)
        all_frequent[[candidate]] <- list(
          itemset = items,
          support = count / n_transactions,
          count = count
        )
      }
    }
    
    if (verbose && length(frequent_candidates) > 0) {
      cat("  - Level", level, ": Found", length(frequent_candidates), "frequent itemsets\n")
    }
    
    current_frequent <- frequent_candidates
    level <- level + 1
  }
  
  return(all_frequent)
}

#' Generate Candidate Itemsets for Apriori
#'
#' @param frequent_itemsets Character vector of frequent itemsets from previous level
#' @param level Current level (itemset size)
#' @return Character vector of candidate itemsets
generate_candidates_apriori <- function(frequent_itemsets, level) {
  
  if (level == 2) {
    # Generate 2-itemsets from 1-itemsets
    candidates <- character(0)
    for (i in 1:(length(frequent_itemsets) - 1)) {
      for (j in (i + 1):length(frequent_itemsets)) {
        candidate <- paste(sort(c(frequent_itemsets[i], frequent_itemsets[j])), collapse = ",")
        candidates <- c(candidates, candidate)
      }
    }
    return(candidates)
  }
  
  # For level > 2, join frequent itemsets that differ by only one item
  candidates <- character(0)
  
  for (i in 1:(length(frequent_itemsets) - 1)) {
    items_i <- sort(strsplit(frequent_itemsets[i], ",")[[1]])
    
    for (j in (i + 1):length(frequent_itemsets)) {
      items_j <- sort(strsplit(frequent_itemsets[j], ",")[[1]])
      
      # Check if they can be joined (differ by only the last item)
      if (length(items_i) == length(items_j) && 
          all(items_i[1:(length(items_i) - 1)] == items_j[1:(length(items_j) - 1)])) {
        
        candidate_items <- sort(unique(c(items_i, items_j)))
        if (length(candidate_items) == level) {
          candidate <- paste(candidate_items, collapse = ",")
          candidates <- c(candidates, candidate)
        }
      }
    }
  }
  
  return(unique(candidates))
}

# ==============================================================================
# FP-GROWTH ALGORITHM CORE FUNCTIONS  
# ==============================================================================

#' Build FP-Tree Structure
#'
#' @param transactions List of transactions
#' @param min_support Minimum support threshold
#' @param verbose Whether to show progress
#' @return FP-tree structure
build_fp_tree <- function(transactions, min_support, verbose = TRUE) {
  
  n_transactions <- length(transactions)
  min_count <- ceiling(min_support * n_transactions)
  
  # Count item frequencies
  item_counts <- table(unlist(transactions))
  frequent_items <- names(item_counts[item_counts >= min_count])
  
  if (length(frequent_items) == 0) {
    return(list(root = list(), header_table = list(), item_counts = item_counts))
  }
  
  # Sort items by frequency (descending)
  frequent_items <- frequent_items[order(item_counts[frequent_items], decreasing = TRUE)]
  
  if (verbose) cat("  - Found", length(frequent_items), "frequent items\n")
  
  # Initialize FP-tree
  fp_tree <- list(
    root = list(item = "root", count = 0, children = list(), parent = NULL),
    header_table = list(),
    item_counts = item_counts[frequent_items]
  )
  
  # Build header table
  for (item in frequent_items) {
    fp_tree$header_table[[item]] <- list(count = item_counts[item], nodes = list())
  }
  
  # Insert transactions into FP-tree
  for (transaction in transactions) {
    # Filter and sort items by frequency
    filtered_items <- intersect(frequent_items, transaction)
    if (length(filtered_items) > 0) {
      insert_transaction_fp(fp_tree, filtered_items)
    }
  }
  
  return(fp_tree)
}

#' Insert Transaction into FP-Tree
#'
#' @param fp_tree FP-tree structure
#' @param items Sorted items to insert
insert_transaction_fp <- function(fp_tree, items) {
  current_node <- fp_tree$root
  
  for (item in items) {
    # Check if item exists in children
    if (item %in% names(current_node$children)) {
      current_node <- current_node$children[[item]]
      current_node$count <- current_node$count + 1
    } else {
      # Create new node
      new_node <- list(
        item = item,
        count = 1,
        children = list(),
        parent = current_node
      )
      current_node$children[[item]] <- new_node
      
      # Add to header table
      fp_tree$header_table[[item]]$nodes <- c(fp_tree$header_table[[item]]$nodes, list(new_node))
      
      current_node <- new_node
    }
  }
}

#' Mine FP-Tree for Frequent Patterns
#'
#' @param fp_tree FP-tree structure
#' @param min_support Minimum support threshold
#' @param verbose Whether to show progress
#' @return List of frequent patterns
mine_fp_tree <- function(fp_tree, min_support, verbose = TRUE) {
  
  frequent_patterns <- list()
  items <- names(fp_tree$header_table)
  
  # Mine patterns for each item (bottom-up)
  for (i in length(items):1) {
    item <- items[i]
    
    # Generate conditional pattern base
    conditional_patterns <- generate_conditional_patterns(fp_tree, item)
    
    if (length(conditional_patterns) > 0) {
      # Build conditional FP-tree
      conditional_tree <- build_conditional_fp_tree(conditional_patterns, min_support)
      
      # Add single item pattern
      support <- fp_tree$header_table[[item]]$count / sum(fp_tree$item_counts)
      frequent_patterns[[item]] <- list(
        pattern = item,
        support = support,
        count = fp_tree$header_table[[item]]$count
      )
      
      # Recursively mine conditional tree
      if (length(conditional_tree$header_table) > 0) {
        conditional_patterns_result <- mine_fp_tree(conditional_tree, min_support, FALSE)
        
        # Add item to each conditional pattern
        for (pattern_name in names(conditional_patterns_result)) {
          pattern <- conditional_patterns_result[[pattern_name]]
          new_pattern <- paste(c(pattern$pattern, item), collapse = ",")
          frequent_patterns[[new_pattern]] <- list(
            pattern = c(pattern$pattern, item),
            support = pattern$support,
            count = pattern$count
          )
        }
      }
    }
  }
  
  if (verbose && length(frequent_patterns) > 0) {
    cat("  - Found", length(frequent_patterns), "frequent patterns\n")
  }
  
  return(frequent_patterns)
}

#' Generate Conditional Pattern Base
#'
#' @param fp_tree FP-tree structure
#' @param item Target item
#' @return List of conditional patterns
generate_conditional_patterns <- function(fp_tree, item) {
  
  conditional_patterns <- list()
  nodes <- fp_tree$header_table[[item]]$nodes
  
  for (node in nodes) {
    if (node$count > 0) {
      # Trace path to root
      path <- character(0)
      current <- node$parent
      
      while (!is.null(current) && current$item != "root") {
        path <- c(current$item, path)
        current <- current$parent
      }
      
      if (length(path) > 0) {
        conditional_patterns[[length(conditional_patterns) + 1]] <- list(
          pattern = path,
          count = node$count
        )
      }
    }
  }
  
  return(conditional_patterns)
}

#' Build Conditional FP-Tree
#'
#' @param conditional_patterns List of conditional patterns
#' @param min_support Minimum support threshold
#' @return Conditional FP-tree
build_conditional_fp_tree <- function(conditional_patterns, min_support) {
  
  if (length(conditional_patterns) == 0) {
    return(list(root = list(), header_table = list(), item_counts = numeric(0)))
  }
  
  # Count item frequencies in conditional patterns
  item_counts <- table(unlist(lapply(conditional_patterns, function(x) {
    rep(x$pattern, x$count)
  })))
  
  total_count <- sum(sapply(conditional_patterns, function(x) x$count))
  min_count <- ceiling(min_support * total_count)
  
  frequent_items <- names(item_counts[item_counts >= min_count])
  
  if (length(frequent_items) == 0) {
    return(list(root = list(), header_table = list(), item_counts = numeric(0)))
  }
  
  # Sort by frequency
  frequent_items <- frequent_items[order(item_counts[frequent_items], decreasing = TRUE)]
  
  # Build conditional FP-tree
  conditional_tree <- list(
    root = list(item = "root", count = 0, children = list(), parent = NULL),
    header_table = list(),
    item_counts = item_counts[frequent_items]
  )
  
  # Initialize header table
  for (item in frequent_items) {
    conditional_tree$header_table[[item]] <- list(count = item_counts[item], nodes = list())
  }
  
  # Insert conditional patterns
  for (pattern_info in conditional_patterns) {
    filtered_items <- intersect(frequent_items, pattern_info$pattern)
    if (length(filtered_items) > 0) {
      for (i in 1:pattern_info$count) {
        insert_transaction_fp(conditional_tree, filtered_items)
      }
    }
  }
  
  return(conditional_tree)
}

# ==============================================================================
# RULE GENERATION FUNCTIONS
# ==============================================================================

#' Generate Association Rules from Frequent Itemsets
#'
#' @param frequent_itemsets List of frequent itemsets
#' @param trans_matrix Transaction matrix
#' @param min_confidence Minimum confidence threshold
#' @param min_lift Minimum lift threshold
#' @return Data frame of association rules
generate_association_rules <- function(frequent_itemsets, trans_matrix, min_confidence, min_lift) {
  
  rules <- data.frame(
    antecedent = character(0),
    consequent = character(0),
    support = numeric(0),
    confidence = numeric(0),
    lift = numeric(0),
    conviction = numeric(0),
    count = integer(0),
    stringsAsFactors = FALSE
  )
  
  # Generate rules from itemsets with 2+ items
  for (itemset_name in names(frequent_itemsets)) {
    itemset_info <- frequent_itemsets[[itemset_name]]
    items <- itemset_info$itemset
    
    if (length(items) >= 2) {
      # Generate all possible antecedent-consequent combinations
      for (i in 1:(length(items) - 1)) {
        antecedent_combinations <- combn(items, i, simplify = FALSE)
        
        for (antecedent in antecedent_combinations) {
          consequent <- setdiff(items, antecedent)
          
          # Calculate metrics
          metrics <- calculate_rule_metrics(antecedent, consequent, trans_matrix)
          
          # Check thresholds
          if (metrics$confidence >= min_confidence && metrics$lift >= min_lift) {
            rules <- rbind(rules, data.frame(
              antecedent = paste(sort(antecedent), collapse = ","),
              consequent = paste(sort(consequent), collapse = ","),
              support = metrics$support,
              confidence = metrics$confidence,
              lift = metrics$lift,
              conviction = metrics$conviction,
              count = metrics$count,
              stringsAsFactors = FALSE
            ))
          }
        }
      }
    }
  }
  
  # Sort by lift (descending)
  if (nrow(rules) > 0) {
    rules <- rules[order(rules$lift, decreasing = TRUE), ]
    rownames(rules) <- NULL
  }
  
  return(rules)
}

#' Generate Rules from FP-Growth Patterns
#'
#' @param frequent_patterns List of frequent patterns
#' @param trans_matrix Transaction matrix
#' @param min_confidence Minimum confidence threshold
#' @param min_lift Minimum lift threshold
#' @return Data frame of association rules
generate_rules_from_patterns <- function(frequent_patterns, trans_matrix, min_confidence, min_lift) {
  
  # Convert patterns to itemsets format for rule generation
  frequent_itemsets <- list()
  
  for (pattern_name in names(frequent_patterns)) {
    pattern_info <- frequent_patterns[[pattern_name]]
    if (is.character(pattern_info$pattern)) {
      items <- pattern_info$pattern
    } else {
      items <- pattern_name
    }
    
    frequent_itemsets[[pattern_name]] <- list(
      itemset = if (length(items) == 1) items else strsplit(items, ",")[[1]],
      support = pattern_info$support,
      count = pattern_info$count
    )
  }
  
  # Use existing rule generation function
  return(generate_association_rules(frequent_itemsets, trans_matrix, min_confidence, min_lift))
}

# ==============================================================================
# UTILITY FUNCTIONS
# ==============================================================================

#' Create Empty Rules Object
#'
#' @param algorithm Algorithm name
#' @return Empty association rules object
create_empty_rules_object <- function(algorithm) {
  
  empty_rules <- data.frame(
    antecedent = character(0),
    consequent = character(0),
    support = numeric(0),
    confidence = numeric(0),
    lift = numeric(0),
    conviction = numeric(0),
    count = integer(0),
    stringsAsFactors = FALSE
  )
  
  result <- list(
    rules = empty_rules,
    algorithm = algorithm,
    parameters = list(),
    summary = list(
      n_transactions = 0,
      n_items = 0,
      n_frequent_itemsets = 0,
      n_rules = 0
    ),
    item_names = character(0)
  )
  
  class(result) <- "association_rules"
  return(result)
}

# ==============================================================================
# S3 METHODS
# ==============================================================================

#' Print Method for Association Rules
#'
#' @param x Association rules object
#' @param ... Additional arguments
#' @export
print.association_rules <- function(x, ...) {
  cat("Association Rules (", x$algorithm, " algorithm)\n", sep = "")
  cat("==========================================\n\n")
  
  cat("Parameters:\n")
  cat("  Min Support:", x$parameters$min_support, "\n")
  cat("  Min Confidence:", x$parameters$min_confidence, "\n")
  cat("  Min Lift:", x$parameters$min_lift, "\n\n")
  
  cat("Summary:\n")
  cat("  Transactions:", x$summary$n_transactions, "\n")
  cat("  Items:", x$summary$n_items, "\n")
  if ("n_frequent_itemsets" %in% names(x$summary)) {
    cat("  Frequent Itemsets:", x$summary$n_frequent_itemsets, "\n")
  } else {
    cat("  Frequent Patterns:", x$summary$n_frequent_patterns, "\n")
  }
  cat("  Rules Generated:", x$summary$n_rules, "\n\n")
  
  if (nrow(x$rules) > 0) {
    cat("Top 10 Rules by Lift:\n")
    top_rules <- head(x$rules[order(x$rules$lift, decreasing = TRUE), ], 10)
    for (i in 1:nrow(top_rules)) {
      cat(sprintf("%2d. %s => %s\n", i, top_rules$antecedent[i], top_rules$consequent[i]))
      cat(sprintf("    Support: %.3f, Confidence: %.3f, Lift: %.3f\n", 
                 top_rules$support[i], top_rules$confidence[i], top_rules$lift[i]))
    }
  } else {
    cat("No rules found with the given thresholds.\n")
  }
}

#' Summary Method for Association Rules
#'
#' @param object Association rules object
#' @param ... Additional arguments
#' @export
summary.association_rules <- function(object, ...) {
  
  cat("Association Rules Summary\n")
  cat("========================\n\n")
  
  cat("Algorithm:", object$algorithm, "\n")
  cat("Transactions:", object$summary$n_transactions, "\n")
  cat("Items:", object$summary$n_items, "\n")
  cat("Rules:", object$summary$n_rules, "\n\n")
  
  if (nrow(object$rules) > 0) {
    cat("Rule Quality Metrics:\n")
    cat("---------------------\n")
    
    metrics <- c("support", "confidence", "lift", "conviction")
    for (metric in metrics) {
      if (metric %in% names(object$rules)) {
        values <- object$rules[[metric]]
        cat(sprintf("%-10s: Min=%.3f, Mean=%.3f, Max=%.3f\n", 
                   toupper(metric), min(values), mean(values), max(values)))
      }
    }
    
    cat("\nRule Length Distribution:\n")
    antecedent_lengths <- sapply(strsplit(object$rules$antecedent, ","), length)
    consequent_lengths <- sapply(strsplit(object$rules$consequent, ","), length)
    total_lengths <- antecedent_lengths + consequent_lengths
    
    length_table <- table(total_lengths)
    for (len in names(length_table)) {
      cat(sprintf("  %s items: %d rules\n", len, length_table[len]))
    }
  }
}

#' Head Method for Association Rules
#'
#' @param x Association rules object
#' @param n Number of rules to show
#' @param ... Additional arguments
#' @export
head.association_rules <- function(x, n = 6, ...) {
  if (nrow(x$rules) == 0) {
    cat("No rules to display.\n")
    return(invisible(NULL))
  }
  
  top_rules <- head(x$rules, n)
  
  cat("Top", min(n, nrow(x$rules)), "Association Rules:\n")
  cat(rep("=", 50), "\n", sep = "")
  
  for (i in 1:nrow(top_rules)) {
    cat(sprintf("%d. %s => %s\n", i, top_rules$antecedent[i], top_rules$consequent[i]))
    cat(sprintf("   Sup: %.3f | Conf: %.3f | Lift: %.3f | Conv: %.3f\n\n",
               top_rules$support[i], top_rules$confidence[i], 
               top_rules$lift[i], top_rules$conviction[i]))
  }
  
  return(invisible(top_rules))
} 
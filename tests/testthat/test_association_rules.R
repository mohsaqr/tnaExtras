# =============================================================================
# TESTS FOR ASSOCIATION RULES FUNCTIONALITY
# =============================================================================

# Test data - learning activities
test_transactions_list <- list(
  c("plan", "discuss", "execute"),
  c("plan", "reflect"),
  c("discuss", "execute", "analyze"),
  c("plan", "discuss", "reflect"),
  c("plan", "execute", "analyze"),
  c("discuss", "reflect", "analyze")
)

test_transactions_df <- data.frame(
  transaction = c(1, 1, 1, 2, 2, 3, 3, 3, 4, 4, 4),
  item = c("A", "B", "C", "A", "D", "B", "C", "E", "A", "B", "D")
)

test_binary_matrix <- matrix(c(
  1, 1, 1, 0, 0,
  1, 0, 0, 1, 0,
  0, 1, 1, 0, 1,
  1, 1, 0, 1, 0,
  1, 0, 1, 0, 1,
  0, 1, 0, 1, 1
), nrow = 6, byrow = TRUE)
colnames(test_binary_matrix) <- c("A", "B", "C", "D", "E")

# =============================================================================
# DATA PREPARATION TESTS
# =============================================================================

test_that("prepare_transactions works with list input", {
  
  result <- prepare_transactions(test_transactions_list)
  
  expect_type(result, "list")
  expect_named(result, c("transactions", "n_transactions", "n_items", "items", "avg_transaction_length"))
  expect_equal(result$n_transactions, 6)
  expect_equal(result$n_items, 5)
  expect_setequal(result$items, c("A", "B", "C", "D", "E"))
  expect_gt(result$avg_transaction_length, 0)
})

test_that("prepare_transactions works with data.frame input", {
  
  result <- prepare_transactions(test_transactions_df, 
                               transaction_col = "transaction", 
                               item_cols = "item")
  
  expect_type(result, "list")
  expect_equal(result$n_transactions, 4)
  expect_equal(result$n_items, 5)
  expect_setequal(result$items, c("A", "B", "C", "D", "E"))
})

test_that("prepare_transactions auto-detects data.frame format", {
  
  result <- prepare_transactions(test_transactions_df)
  
  expect_type(result, "list")
  expect_equal(result$n_transactions, 4)
  expect_equal(result$n_items, 5)
})

test_that("prepare_transactions works with matrix input", {
  
  result <- prepare_transactions(test_binary_matrix)
  
  expect_type(result, "list")
  expect_equal(result$n_transactions, 6)
  expect_equal(result$n_items, 5)
  expect_setequal(result$items, c("A", "B", "C", "D", "E"))
})

test_that("prepare_transactions handles errors correctly", {
  
  expect_error(prepare_transactions("invalid"), "Data must be a list, matrix, or data.frame")
  expect_error(prepare_transactions(list()), "No valid transactions found")
  
  # Test with non-binary matrix
  non_binary_matrix <- matrix(1:10, nrow = 2)
  expect_error(prepare_transactions(non_binary_matrix), "Matrix must be binary")
})

test_that("create_transaction_matrix works correctly", {
  
  transactions <- list(c("A", "B"), c("B", "C"), c("A", "C"))
  result <- create_transaction_matrix(transactions)
  
  expect_true(is.matrix(result))
  expect_equal(nrow(result), 3)
  expect_equal(ncol(result), 3)
  expect_setequal(colnames(result), c("A", "B", "C"))
  expect_true(all(result %in% c(TRUE, FALSE)))
})

# =============================================================================
# PARAMETER VALIDATION TESTS
# =============================================================================

test_that("validate_parameters works correctly", {
  
  expect_silent(validate_parameters(0.1, 0.5, 1.0))
  expect_silent(validate_parameters(0, 0, 0))
  expect_silent(validate_parameters(1, 1, 10))
  
  expect_error(validate_parameters(-0.1, 0.5, 1.0), "min_support must be a single number between 0 and 1")
  expect_error(validate_parameters(1.1, 0.5, 1.0), "min_support must be a single number between 0 and 1")
  expect_error(validate_parameters(0.5, -0.1, 1.0), "min_confidence must be a single number between 0 and 1")
  expect_error(validate_parameters(0.5, 1.1, 1.0), "min_confidence must be a single number between 0 and 1")
  expect_error(validate_parameters(0.5, 0.5, -1.0), "min_lift must be a single non-negative number")
  
  expect_warning(validate_parameters(0, 0.5, 1.0), "min_support = 0 may result in very large number of patterns")
  expect_warning(validate_parameters(0.5, 0, 1.0), "min_confidence = 0 may result in meaningless rules")
})

# =============================================================================
# ASSOCIATION RULE MINING TESTS
# =============================================================================

test_that("apriori_rules works with basic input", {
  
  result <- apriori_rules(test_transactions_list, 
                         min_support = 0.3, 
                         min_confidence = 0.5, 
                         verbose = FALSE)
  
  expect_s3_class(result, "association_rules")
  expect_named(result, c("rules", "algorithm", "parameters", "summary", "item_names"))
  expect_equal(result$algorithm, "apriori")
  expect_true(is.data.frame(result$rules))
  expect_named(result$rules, c("antecedent", "consequent", "support", "confidence", "lift", "conviction", "count"))
})

test_that("fp_growth_rules works with basic input", {
  
  result <- fp_growth_rules(test_transactions_list, 
                           min_support = 0.3, 
                           min_confidence = 0.5, 
                           verbose = FALSE)
  
  expect_s3_class(result, "association_rules")
  expect_named(result, c("rules", "algorithm", "parameters", "summary", "item_names"))
  expect_equal(result$algorithm, "fp_growth")
  expect_true(is.data.frame(result$rules))
  expect_named(result$rules, c("antecedent", "consequent", "support", "confidence", "lift", "conviction", "count"))
})

test_that("association rule algorithms handle edge cases", {
  
  # High support threshold (no rules expected)
  result_high_support <- apriori_rules(test_transactions_list, 
                                      min_support = 0.9, 
                                      verbose = FALSE)
  expect_equal(nrow(result_high_support$rules), 0)
  
  # High confidence threshold
  result_high_confidence <- fp_growth_rules(test_transactions_list, 
                                           min_confidence = 0.99, 
                                           verbose = FALSE)
  expect_gte(nrow(result_high_confidence$rules), 0)
  
  # Single transaction
  single_transaction <- list(c("A", "B", "C"))
  result_single <- apriori_rules(single_transaction, verbose = FALSE)
  expect_gte(nrow(result_single$rules), 0)
})

test_that("association rule algorithms produce valid metrics", {
  
  result <- apriori_rules(test_transactions_list, 
                         min_support = 0.2, 
                         min_confidence = 0.5, 
                         verbose = FALSE)
  
  if (nrow(result$rules) > 0) {
    # Support should be between 0 and 1
    expect_true(all(result$rules$support >= 0 & result$rules$support <= 1))
    
    # Confidence should be between 0 and 1
    expect_true(all(result$rules$confidence >= 0 & result$rules$confidence <= 1))
    
    # Lift should be positive
    expect_true(all(result$rules$lift >= 0))
    
    # Count should be positive integers
    expect_true(all(result$rules$count >= 0))
    expect_true(all(result$rules$count == round(result$rules$count)))
  }
})

# =============================================================================
# RULE METRICS CALCULATION TESTS
# =============================================================================

test_that("calculate_rule_metrics works correctly", {
  
  trans_matrix <- create_transaction_matrix(test_transactions_list)
  
  # Test simple rule: A => B
  metrics <- calculate_rule_metrics(c("A"), c("B"), trans_matrix)
  
  expect_named(metrics, c("support", "confidence", "lift", "conviction", "count", 
                         "antecedent_count", "consequent_count"))
  expect_true(all(sapply(metrics, is.numeric)))
  expect_gte(metrics$support, 0)
  expect_lte(metrics$support, 1)
  expect_gte(metrics$confidence, 0)
  expect_lte(metrics$confidence, 1)
  expect_gte(metrics$lift, 0)
})

test_that("calculate_itemset_support works correctly", {
  
  trans_matrix <- create_transaction_matrix(test_transactions_list)
  
  # Single item
  support_A <- calculate_itemset_support("A", trans_matrix)
  expect_gte(support_A, 0)
  expect_lte(support_A, 1)
  
  # Multiple items
  support_AB <- calculate_itemset_support(c("A", "B"), trans_matrix)
  expect_gte(support_AB, 0)
  expect_lte(support_AB, 1)
  expect_lte(support_AB, support_A)  # AB support <= A support
  
  # Empty itemset
  support_empty <- calculate_itemset_support(character(0), trans_matrix)
  expect_equal(support_empty, 0)
})

# =============================================================================
# RULE FILTERING AND RANKING TESTS
# =============================================================================

test_that("filter_association_rules works correctly", {
  
  # Create some test rules
  test_rules <- data.frame(
    antecedent = c("A", "B", "C"),
    consequent = c("B", "C", "A"),
    support = c(0.5, 0.3, 0.8),
    confidence = c(0.7, 0.9, 0.6),
    lift = c(1.2, 2.0, 1.5),
    conviction = c(2.0, 3.0, 1.8),
    count = c(5, 3, 8),
    stringsAsFactors = FALSE
  )
  
  # Filter by support
  filtered_support <- filter_association_rules(test_rules, min_support = 0.4)
  expect_lte(nrow(filtered_support), nrow(test_rules))
  expect_true(all(filtered_support$support >= 0.4))
  
  # Filter by confidence
  filtered_confidence <- filter_association_rules(test_rules, min_confidence = 0.8)
  expect_true(all(filtered_confidence$confidence >= 0.8))
  
  # Filter by lift
  filtered_lift <- filter_association_rules(test_rules, min_lift = 1.5)
  expect_true(all(filtered_lift$lift >= 1.5))
  
  # Filter by count
  filtered_count <- filter_association_rules(test_rules, min_count = 4)
  expect_true(all(filtered_count$count >= 4))
  
  # Multiple filters
  filtered_multiple <- filter_association_rules(test_rules, 
                                               min_support = 0.4, 
                                               min_confidence = 0.6)
  expect_true(all(filtered_multiple$support >= 0.4))
  expect_true(all(filtered_multiple$confidence >= 0.6))
})

test_that("rank_association_rules works correctly", {
  
  test_rules <- data.frame(
    antecedent = c("A", "B", "C"),
    consequent = c("B", "C", "A"),
    support = c(0.3, 0.5, 0.8),
    confidence = c(0.7, 0.9, 0.6),
    lift = c(1.2, 2.0, 1.5),
    stringsAsFactors = FALSE
  )
  
  # Rank by support (descending)
  ranked_support <- rank_association_rules(test_rules, by = "support", decreasing = TRUE)
  expect_equal(ranked_support$support, sort(test_rules$support, decreasing = TRUE))
  
  # Rank by confidence (ascending)
  ranked_confidence <- rank_association_rules(test_rules, by = "confidence", decreasing = FALSE)
  expect_equal(ranked_confidence$confidence, sort(test_rules$confidence, decreasing = FALSE))
  
  # Invalid metric
  expect_error(rank_association_rules(test_rules, by = "invalid"), 
               "Metric 'invalid' not found in rules")
})

test_that("extract_rules_by_item works correctly", {
  
  test_rules <- data.frame(
    antecedent = c("A", "B,C", "D"),
    consequent = c("B", "A", "E"),
    support = c(0.5, 0.3, 0.8),
    confidence = c(0.7, 0.9, 0.6),
    lift = c(1.2, 2.0, 1.5),
    stringsAsFactors = FALSE
  )
  
  # Extract rules with item A
  rules_with_A <- extract_rules_by_item(test_rules, "A")
  expect_true(nrow(rules_with_A) >= 1)
  
  # Extract rules with A in antecedent only
  rules_A_antecedent <- extract_rules_by_item(test_rules, "A", side = "antecedent")
  expect_true(all(grepl("A", rules_A_antecedent$antecedent)))
  
  # Extract rules with A in consequent only
  rules_A_consequent <- extract_rules_by_item(test_rules, "A", side = "consequent")
  expect_true(all(grepl("A", rules_A_consequent$consequent)))
  
  # Extract with non-existent item
  rules_empty <- extract_rules_by_item(test_rules, "Z")
  expect_equal(nrow(rules_empty), 0)
})

# =============================================================================
# UTILITY FUNCTION TESTS
# =============================================================================

test_that("find_redundant_rules works correctly", {
  
  test_rules <- data.frame(
    antecedent = c("A", "A,B", "C"),
    consequent = c("B", "C", "D"),
    support = c(0.5, 0.3, 0.8),
    confidence = c(0.7, 0.8, 0.6),
    lift = c(1.2, 2.0, 1.5),
    stringsAsFactors = FALSE
  )
  
  redundant <- find_redundant_rules(test_rules)
  expect_type(redundant, "logical")
  expect_equal(length(redundant), nrow(test_rules))
})

test_that("calculate_rule_overlap works correctly", {
  
  test_rules <- data.frame(
    antecedent = c("A", "B", "C"),
    consequent = c("B", "C", "A"),
    support = c(0.5, 0.3, 0.8),
    stringsAsFactors = FALSE
  )
  
  # Jaccard similarity
  overlap_jaccard <- calculate_rule_overlap(test_rules, method = "jaccard")
  expect_true(is.matrix(overlap_jaccard))
  expect_equal(nrow(overlap_jaccard), nrow(test_rules))
  expect_equal(ncol(overlap_jaccard), nrow(test_rules))
  expect_true(all(overlap_jaccard >= 0 & overlap_jaccard <= 1))
  expect_true(all(diag(overlap_jaccard) == 1))  # Self-similarity should be 1
  
  # Dice similarity
  overlap_dice <- calculate_rule_overlap(test_rules, method = "dice")
  expect_true(is.matrix(overlap_dice))
  expect_true(all(overlap_dice >= 0 & overlap_dice <= 1))
  
  # Overlap method
  overlap_overlap <- calculate_rule_overlap(test_rules, method = "overlap")
  expect_true(is.matrix(overlap_overlap))
  expect_true(all(overlap_overlap >= 0 & overlap_overlap <= 1))
  
  # Invalid method
  expect_error(calculate_rule_overlap(test_rules, method = "invalid"), 
               "method must be 'jaccard', 'overlap', or 'dice'")
})

test_that("rules_to_transactions works correctly", {
  
  test_rules <- data.frame(
    antecedent = c("A", "B"),
    consequent = c("B", "C"),
    support = c(0.5, 0.3),
    confidence = c(0.7, 0.9),
    lift = c(1.2, 2.0),
    stringsAsFactors = FALSE
  )
  
  # Without metrics
  transactions <- rules_to_transactions(test_rules, include_metrics = FALSE)
  expect_type(transactions, "list")
  expect_equal(length(transactions), nrow(test_rules))
  expect_true(all(sapply(transactions, is.character)))
  
  # With metrics
  transactions_with_metrics <- rules_to_transactions(test_rules, include_metrics = TRUE)
  expect_true(length(transactions_with_metrics[[1]]) > length(transactions[[1]]))
})

# =============================================================================
# ALGORITHM COMPARISON TESTS
# =============================================================================

test_that("compare_rule_algorithms works correctly", {
  
  apriori_result <- apriori_rules(test_transactions_list, 
                                 min_support = 0.2, 
                                 verbose = FALSE)
  fp_growth_result <- fp_growth_rules(test_transactions_list, 
                                     min_support = 0.2, 
                                     verbose = FALSE)
  
  # Test comparison
  expect_output(comparison <- compare_rule_algorithms(
    list(Apriori = apriori_result, FP_Growth = fp_growth_result)
  ))
  
  expect_type(comparison, "data.frame")
  expect_true("Algorithm" %in% names(comparison))
  expect_true("N_Rules" %in% names(comparison))
  
  # Test with insufficient inputs
  expect_error(compare_rule_algorithms(list(apriori_result)), 
               "rules_list must be a list with at least 2")
  
  # Test with invalid objects
  expect_error(compare_rule_algorithms(list("invalid", "objects")), 
               "All objects must be of class 'association_rules'")
})

# =============================================================================
# S3 METHODS TESTS
# =============================================================================

test_that("print.association_rules works correctly", {
  
  result <- apriori_rules(test_transactions_list, verbose = FALSE)
  
  expect_output(print(result), "Association Rules")
  expect_output(print(result), "Parameters:")
  expect_output(print(result), "Summary:")
})

test_that("summary.association_rules works correctly", {
  
  result <- apriori_rules(test_transactions_list, verbose = FALSE)
  
  expect_output(summary(result), "Association Rules Summary")
  expect_output(summary(result), "Algorithm:")
  expect_output(summary(result), "Rules:")
})

test_that("head.association_rules works correctly", {
  
  result <- apriori_rules(test_transactions_list, verbose = FALSE)
  
  if (nrow(result$rules) > 0) {
    expect_output(head_result <- head(result, 3), "Association Rules")
    expect_true(is.data.frame(head_result))
    expect_lte(nrow(head_result), 3)
  }
})

# =============================================================================
# ERROR HANDLING TESTS
# =============================================================================

test_that("association rule functions handle errors correctly", {
  
  # Missing required parameter
  expect_error(apriori_rules(), "transactions parameter is required")
  expect_error(fp_growth_rules(), "transactions parameter is required")
  
  # Invalid parameter values
  expect_error(apriori_rules(test_transactions_list, min_support = -0.1), 
               "min_support must be a single number between 0 and 1")
  expect_error(fp_growth_rules(test_transactions_list, min_confidence = 1.5), 
               "min_confidence must be a single number between 0 and 1")
  
  # Empty input
  expect_warning(apriori_rules(list(), verbose = FALSE), 
                 "No frequent itemsets found")
})

# =============================================================================
# EXPORT FUNCTION TESTS
# =============================================================================

test_that("export_association_rules works correctly", {
  
  result <- apriori_rules(test_transactions_list, verbose = FALSE)
  
  # Test with association_rules object
  temp_file <- tempfile(fileext = ".csv")
  
  if (nrow(result$rules) > 0) {
    expect_silent(export_association_rules(result, temp_file, format = "csv"))
    expect_true(file.exists(temp_file))
    unlink(temp_file)
  } else {
    expect_warning(export_association_rules(result, temp_file), "No rules to export")
  }
  
  # Test with data frame
  if (nrow(result$rules) > 0) {
    temp_file2 <- tempfile(fileext = ".json")
    expect_silent(export_association_rules(result$rules, temp_file2, format = "json"))
    expect_true(file.exists(temp_file2))
    unlink(temp_file2)
  }
  
  # Test invalid format
  expect_error(export_association_rules(result, "test.xyz", format = "xyz"), 
               "format must be 'csv', 'json', or 'txt'")
  
  # Test invalid input
  expect_error(export_association_rules("invalid", "test.csv"), 
               "rules must be an association_rules object or data frame")
})

# =============================================================================
# EDGE CASE AND PERFORMANCE TESTS
# =============================================================================

test_that("association rules handle edge cases", {
  
  # Empty transactions
  empty_transactions <- list()
  expect_error(prepare_transactions(empty_transactions), "No valid transactions found")
  
  # Single item transactions
  single_item_transactions <- list(c("A"), c("B"), c("A"))
  result_single <- apriori_rules(single_item_transactions, verbose = FALSE)
  expect_s3_class(result_single, "association_rules")
  
  # Identical transactions
  identical_transactions <- list(c("A", "B"), c("A", "B"), c("A", "B"))
  result_identical <- fp_growth_rules(identical_transactions, verbose = FALSE)
  expect_s3_class(result_identical, "association_rules")
  
  # Very long transactions
  long_transaction <- list(c(letters[1:20]))
  result_long <- apriori_rules(long_transaction, verbose = FALSE)
  expect_s3_class(result_long, "association_rules")
})

test_that("association rules work with various data types", {
  
  # Factor data
  factor_transactions <- list(
    as.factor(c("A", "B", "C")),
    as.factor(c("A", "D")),
    as.factor(c("B", "C"))
  )
  result_factor <- apriori_rules(factor_transactions, verbose = FALSE)
  expect_s3_class(result_factor, "association_rules")
  
  # Mixed character/factor data frame
  mixed_df <- data.frame(
    transaction = as.factor(c(1, 1, 2, 2, 3, 3)),
    item = c("A", "B", "A", "C", "B", "C"),
    stringsAsFactors = FALSE
  )
  result_mixed <- fp_growth_rules(mixed_df, verbose = FALSE)
  expect_s3_class(result_mixed, "association_rules")
})

cat("All association rules tests completed successfully!\n") 
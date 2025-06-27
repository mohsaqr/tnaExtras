test_that("compute_distance_matrix works with different methods", {
  # Test data
  data <- data.frame(
    T1 = c("A", "B", "A", "C"),
    T2 = c("B", "A", "B", "A"),
    T3 = c("C", "C", "A", "B")
  )
  
  # Test euclidean distance
  dist_euclidean <- compute_distance_matrix(data, "euclidean")
  expect_s3_class(dist_euclidean, "dist")
  expect_equal(length(dist_euclidean), 6)  # 4 choose 2
  
  # Test hamming distance
  dist_hamming <- compute_distance_matrix(data, "hamming")
  expect_s3_class(dist_hamming, "dist")
  expect_equal(length(dist_hamming), 6)
  
  # Test LCS distance
  dist_lcs <- compute_distance_matrix(data, "lcs")
  expect_s3_class(dist_lcs, "dist")
  expect_equal(length(dist_lcs), 6)
  
  # Test transition distance
  dist_transition <- compute_distance_matrix(data, "transition")
  expect_s3_class(dist_transition, "dist")
  expect_equal(length(dist_transition), 6)
  
  # Test stringdist methods
  dist_lv <- compute_distance_matrix(data, "lv")
  expect_s3_class(dist_lv, "dist")
  expect_equal(length(dist_lv), 6)
  
  # Test jaro-winkler with parameter
  dist_jw <- compute_distance_matrix(data, "jw", p = 0.2)
  expect_s3_class(dist_jw, "dist")
  expect_equal(length(dist_jw), 6)
})

test_that("cluster_sequences works with different methods", {
  # Test data
  data <- data.frame(
    T1 = c("A", "B", "A", "C", "A", "B"),
    T2 = c("B", "A", "B", "A", "C", "A"),
    T3 = c("C", "C", "A", "B", "B", "C")
  )
  
  # Test PAM clustering
  result_pam <- cluster_sequences(data, k = 2, distance_method = "euclidean", clustering_method = "pam")
  expect_type(result_pam, "list")
  expect_equal(result_pam$k, 2)
  expect_equal(length(result_pam$assignments), 6)
  expect_true(all(result_pam$assignments %in% 1:2))
  expect_s3_class(result_pam, "tna_cluster_result")
  
  # Test hierarchical clustering
  result_ward <- cluster_sequences(data, k = 2, distance_method = "hamming", clustering_method = "ward.D2")
  expect_type(result_ward, "list")
  expect_equal(result_ward$k, 2)
  expect_equal(length(result_ward$assignments), 6)
  expect_true(all(result_ward$assignments %in% 1:2))
  
  # Test with different distance methods
  result_lcs <- cluster_sequences(data, k = 2, distance_method = "lcs")
  expect_type(result_lcs, "list")
  expect_equal(result_lcs$k, 2)
})

test_that("cluster_mmm works correctly", {
  # Test data
  data <- data.frame(
    T1 = c("A", "B", "A", "C", "A", "B"),
    T2 = c("B", "A", "B", "A", "C", "A"),
    T3 = c("C", "C", "A", "B", "B", "C")
  )
  
  # Test basic MMM clustering
  result <- cluster_mmm(data, k = 2, n_starts = 2, max_iter = 50, verbose = FALSE, seed = 123)
  expect_s3_class(result, "mmm_result")
  expect_equal(result$k, 2)
  expect_equal(length(result$assignments), 6)
  expect_true(all(result$assignments %in% 1:2))
  expect_true(is.numeric(result$log_likelihood))
  expect_true(is.numeric(result$bic))
  expect_true(is.numeric(result$aic))
  expect_equal(length(result$models), 2)
  expect_equal(length(result$mixture_weights), 2)
  expect_true(abs(sum(result$mixture_weights) - 1) < 1e-10)
})

test_that("compare_clustering_methods works", {
  # Test data
  data <- data.frame(
    T1 = c("A", "B", "A", "C", "A", "B"),
    T2 = c("B", "A", "B", "A", "C", "A"),
    T3 = c("C", "C", "A", "B", "B", "C")
  )
  
  # Test comparison
  comparison <- compare_clustering_methods(data, k = 2, distance_method = "euclidean", 
                                          clustering_methods = c("pam", "ward.D2"))
  expect_s3_class(comparison, "data.frame")
  expect_equal(nrow(comparison), 2)
  expect_true("method" %in% names(comparison))
  expect_true("silhouette" %in% names(comparison))
  expect_true("cluster_balance" %in% names(comparison))
})

test_that("find_clusters_range works", {
  # Test data
  data <- data.frame(
    T1 = c("A", "B", "A", "C", "A", "B", "C", "A"),
    T2 = c("B", "A", "B", "A", "C", "A", "A", "C"),
    T3 = c("C", "C", "A", "B", "B", "C", "B", "B")
  )
  
  # Test range finding
  results <- find_clusters_range(data, distance_method = "euclidean", k_range = 2:3, clustering_method = "pam")
  expect_type(results, "list")
  expect_true("all_results" %in% names(results))
  expect_equal(length(results$all_results), 2)  # k=2 and k=3
  expect_true("k2" %in% names(results$all_results))
  expect_true("k3" %in% names(results$all_results))
})

test_that("analyze_sequences works", {
  # Test data
  data <- data.frame(
    T1 = c("A", "B", "A", "C"),
    T2 = c("B", "A", "B", "A"),
    T3 = c("C", "C", "A", "B")
  )
  
  # Test multiple methods
  results <- analyze_sequences(data, methods = c("euclidean", "hamming"))
  expect_type(results, "list")
  expect_equal(length(results), 2)
  expect_true("euclidean" %in% names(results))
  expect_true("hamming" %in% names(results))
  expect_s3_class(results$euclidean, "dist")
  expect_s3_class(results$hamming, "dist")
})

test_that("sequence_complete_analysis works", {
  # Test data
  data <- data.frame(
    T1 = c("A", "B", "A", "C", "A", "B"),
    T2 = c("B", "A", "B", "A", "C", "A"),
    T3 = c("C", "C", "A", "B", "B", "C")
  )
  
  # Test complete analysis
  analysis <- sequence_complete_analysis(data, k_range = 2:3, 
                                        distance_methods = c("euclidean", "hamming"),
                                        clustering_methods = c("pam"))
  expect_type(analysis, "list")
  expect_true("summary" %in% names(analysis))
  expect_s3_class(analysis$summary, "data.frame")
  expect_true(nrow(analysis$summary) >= 4)  # 2 k values × 2 distance methods × 1 clustering method
})

test_that("error handling works correctly", {
  # Test invalid inputs
  expect_error(compute_distance_matrix(c(1, 2, 3), "euclidean"))
  expect_error(compute_distance_matrix(data.frame(T1 = "A"), "euclidean"))  # Only 1 row
  expect_error(compute_distance_matrix(data.frame(T1 = c("A", "B")), "invalid_method"))
  
  # Test invalid k values
  data <- data.frame(T1 = c("A", "B"), T2 = c("B", "A"))
  expect_error(cluster_sequences(data, k = 1))
  expect_error(cluster_sequences(data, k = 3))  # k >= nrow
  expect_error(cluster_mmm(data, k = 1))
  expect_error(cluster_mmm(data, k = 3))
})

test_that("missing value handling works", {
  # Test data with missing values
  data <- data.frame(
    T1 = c("A", "B", NA, "C"),
    T2 = c("B", NA, "B", "A"),
    T3 = c(NA, "C", "A", "B")
  )
  
  # Test distance computation with NAs
  dist_euclidean <- compute_distance_matrix(data, "euclidean")
  expect_s3_class(dist_euclidean, "dist")
  expect_equal(length(dist_euclidean), 6)
  
  # Test clustering with NAs
  result <- cluster_sequences(data, k = 2, distance_method = "euclidean")
  expect_equal(length(result$assignments), 4)
  expect_true(all(result$assignments %in% 1:2))
})

test_that("parameter validation works", {
  data <- data.frame(
    T1 = c("A", "B", "A", "C"),
    T2 = c("B", "A", "B", "A"),
    T3 = c("C", "C", "A", "B")
  )
  
  # Test q parameter validation
  expect_error(compute_distance_matrix(data, "qgram", q = 0))
  expect_error(compute_distance_matrix(data, "qgram", q = -1))
  
  # Test p parameter validation  
  expect_error(compute_distance_matrix(data, "jw", p = -0.1))
  expect_error(compute_distance_matrix(data, "jw", p = 1.1))
  
  # Test cost parameter validation
  expect_error(compute_distance_matrix(data, "optimal_matching", substitution_cost = -1))
  expect_error(compute_distance_matrix(data, "optimal_matching", indel_cost = -1))
})

test_that("clustering with small datasets works", {
  # Minimal test data
  data <- data.frame(
    T1 = c("A", "B"),
    T2 = c("B", "A")
  )
  
  # This should work with k=2 (but will give warning about k being close to n)
  # We expect it to still work for the minimum case
  expect_error(cluster_sequences(data, k = 2))  # k must be < nrow, so k=2 for n=2 should fail
  
  # Test with 3 sequences
  data3 <- data.frame(
    T1 = c("A", "B", "C"),
    T2 = c("B", "A", "A")
  )
  
  result <- cluster_sequences(data3, k = 2)
  expect_equal(length(result$assignments), 3)
  expect_true(all(result$assignments %in% 1:2))
}) 
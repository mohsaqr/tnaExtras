test_that("compute_weighted_distance works (static and dynamic)", {
  data <- data.frame(
    T1 = c("A", "B", "C", "A"),
    T2 = c("B", "C", NA,  "B"),
    T3 = c("C", "A", "B", "C")
  )

  state_w <- matrix(c(
    2.0, 1.0, 1.5,
    1.0, 2.5, 1.0,
    1.5, 1.0, 2.0,
    1.0, 1.0, 1.0
  ), nrow = 4, byrow = TRUE)
  colnames(state_w) <- c("A","B","C")

  time_w <- c(0.5, 1.0, 1.5)

  d_static_eu <- compute_weighted_distance(data, method = "euclidean",
                                           time_weights = time_w, state_weights = state_w)
  expect_s3_class(d_static_eu, "dist")

  d_static_hm <- compute_weighted_distance(data, method = "hamming",
                                           time_weights = time_w, state_weights = state_w)
  expect_s3_class(d_static_hm, "dist")

  d_dynamic <- compute_weighted_distance(data, method = "euclidean",
                                         time_weights = time_w, state_weights = state_w,
                                         decay_rate = 0.5)
  expect_s3_class(d_dynamic, "dist")
})

test_that("cluster_weighted_sequences works (pam and hierarchical)", {
  data <- data.frame(
    T1 = c("A", "B", "C", "A", "A", "B"),
    T2 = c("B", "C", NA,  "B", "B", "A"),
    T3 = c("C", "A", "B", "C", "C", "C")
  )

  state_w <- matrix(c(
    2.0, 1.0, 1.5,
    1.0, 2.5, 1.0,
    1.5, 1.0, 2.0,
    1.0, 1.0, 1.0,
    2.5, 1.0, 1.0,
    1.0, 2.0, 1.5
  ), nrow = 6, byrow = TRUE)
  colnames(state_w) <- c("A","B","C")
  time_w <- c(0.5, 1.0, 1.5)

  res_pam <- cluster_weighted_sequences(
    data = data,
    k = 3,
    distance_method = "hamming",
    clustering_method = "pam",
    time_weights = time_w,
    state_weights = state_w,
    decay_rate = 0.75
  )
  expect_s3_class(res_pam, "tna_cluster_result")
  expect_equal(res_pam$k, 3)

  res_hier <- cluster_weighted_sequences(
    data = data,
    k = 3,
    distance_method = "euclidean",
    clustering_method = "ward.D2",
    time_weights = time_w,
    state_weights = state_w
  )
  expect_s3_class(res_hier, "tna_cluster_result")
  expect_equal(res_hier$k, 3)
})

test_that("weighted_clustering alias works", {
  data <- data.frame(
    T1 = c("A", "B", "C", "A", "A", "B"),
    T2 = c("B", "C", NA,  "B", "B", "A"),
    T3 = c("C", "A", "B", "C", "C", "C")
  )
  state_w <- matrix(rep(c(2,1,1.5), 6), nrow = 6, byrow = TRUE)
  colnames(state_w) <- c("A","B","C")
  time_w <- c(1,1,1)
  res <- weighted_clustering(data, k = 2, distance_method = "hamming",
                             time_weights = time_w, state_weights = state_w)
  expect_s3_class(res, "tna_cluster_result")
  expect_equal(res$k, 2)
})


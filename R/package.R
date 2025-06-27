#' @keywords internal
"_PACKAGE"

#' tnaCluster: Advanced Sequence Analysis and Clustering for Temporal Network Analysis
#'
#' The tnaCluster package provides advanced methods for sequence dissimilarity analysis 
#' and clustering, particularly designed for temporal network analysis and collaborative 
#' learning regulation data. The package implements multiple distance measures including 
#' an innovative Markov-based transition analysis with higher-order patterns, 
#' time-weighted transitions, and pattern complexity measures.
#'
#' @section Main Functions:
#'
#' The main functions for users are:
#' \itemize{
#'   \item \code{\link{compute_distance_matrix}} - Compute dissimilarity matrices using various methods
#'   \item \code{\link{cluster_sequences}} - Perform clustering with specified parameters
#'   \item \code{\link{cluster_mmm}} - Mixture Markov Model clustering with EM algorithm
#'   \item \code{\link{find_optimal_mmm}} - Enhanced model selection for MMM with comprehensive diagnostics
#'   \item \code{\link{find_clusters_range}} - Test multiple k values without optimization
#'   \item \code{\link{compare_clustering_methods}} - Compare different clustering algorithms
#'   \item \code{\link{analyze_sequences}} - Analyze multiple distance methods
#'   \item \code{\link{sequence_complete_analysis}} - Comprehensive analysis workflow
#' }
#'
#' @section Distance Methods:
#'
#' The package supports the following distance methods:
#' \itemize{
#'   \item \strong{Traditional Methods}
#'   \itemize{
#'     \item \code{euclidean} - Fast numeric encoding, good for exploration
#'     \item \code{hamming} - Position-wise differences, handles missing values
#'     \item \code{lcs} - Longest Common Subsequence similarity
#'     \item \code{start_position} - Based on first occurrence of states
#'     \item \code{transition} - Advanced Markov analysis with configurable weights
#'     \item \code{optimal_matching} - Optimal matching with dynamic programming
#'   }
#'   \item \strong{StringDist Methods}
#'   \itemize{
#'     \item \code{osa} - Optimal String Alignment
#'     \item \code{lv} - Levenshtein distance
#'     \item \code{dl} - Damerau-Levenshtein distance
#'     \item \code{qgram} - Q-gram distance (with q parameter)
#'     \item \code{jaro} - Jaro distance
#'     \item \code{jw} - Jaro-Winkler distance (with p parameter)
#'     \item \code{cosine} - Cosine distance on q-grams
#'     \item \code{jaccard} - Jaccard distance on q-grams
#'     \item \code{soundex} - Soundex phonetic distance
#'     \item \code{lcs_stringdist} - LCS via stringdist package
#'   }
#' }
#'
#' @section Clustering Methods:
#'
#' Supports multiple clustering algorithms:
#' \itemize{
#'   \item \strong{pam} - Partitioning Around Medoids (default, recommended)
#'   \item \strong{ward.D2} - Ward's method (recommended for hierarchical)
#'   \item \strong{complete} - Complete linkage, forms compact clusters
#'   \item \strong{average} - Average linkage (UPGMA), balanced approach
#'   \item \strong{single} - Single linkage, can create elongated clusters
#'   \item Additional hierarchical methods: ward.D, mcquitty, median, centroid
#' }
#'
#' @section Key Features:
#'
#' \itemize{
#'   \item \strong{Performance}: Vectorized operations and optimized algorithms
#'   \item \strong{Parallel Processing}: Robust parallel support with graceful fallbacks
#'   \item \strong{Flexibility}: 16 different distance methods with extensive parameterization
#'   \item \strong{Robustness}: Multiple random starts for MMM, comprehensive error handling
#'   \item \strong{Analysis}: Built-in comparison and evaluation functions
#'   \item \strong{Missing Data}: Robust handling of incomplete sequences
#'   \item \strong{Documentation}: Comprehensive examples and clear parameter descriptions
#' }
#'
#' @section Quick Start:
#'
#' \preformatted{
#' # Load example data
#' data <- data.frame(
#'   T1 = c("A", "B", "A", "C"),
#'   T2 = c("B", "A", "B", "A"),
#'   T3 = c("C", "C", "A", "B")
#' )
#' 
#' # Basic clustering
#' result <- cluster_sequences(data, k = 2)
#' 
#' # Compare methods
#' comparison <- compare_clustering_methods(data, k = 2)
#' 
#' # Enhanced MMM clustering with automatic k selection
#' optimal_mmm_result <- find_optimal_mmm(data, k_range = 2:4)
#' 
#' # Parallel processing (faster for multiple k values)
#' optimal_mmm_parallel <- find_optimal_mmm(data, k_range = 2:4, parallel = TRUE)
#' 
#' # Complete analysis with parallel processing
#' complete_analysis <- sequence_complete_analysis(data, parallel = TRUE)
#' print(optimal_mmm_result$best_result)
#' 
#' # MMM clustering with specified k
#' mmm_result <- cluster_mmm(data, k = 2, n_starts = 5)
#' }
#'
#' @author TNA Clustering Team
#' @docType package
#' @name tnaCluster-package
#' @aliases tnaCluster
#' @importFrom cluster pam silhouette
#' @importFrom stats cutree dist hclust setNames
#' @importFrom stringdist stringdist
NULL 
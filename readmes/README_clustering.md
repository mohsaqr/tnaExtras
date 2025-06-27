# tnaCluster: Advanced Sequence Analysis and Clustering

[![CRAN status](https://www.r-pkg.org/badges/version/tnaCluster)](https://CRAN.R-project.org/package=tnaCluster)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

## Overview

**tnaExtras** provides advanced methods for sequence dissimilarity analysis and clustering, particularly designed for temporal network analysis and collaborative learning regulation data. The package implements multiple distance measures including innovative Markov-based transition analysis with higher-order patterns, time-weighted transitions, and pattern complexity measures.

## Key Features

- **Multiple Distance Methods**: Choose from traditional methods (Euclidean, Hamming, LCS) to advanced approaches (Markov transition analysis, optimal matching)
- **StringDist Integration**: Complete support for all `stringdist` package methods
- **Mixture Markov Models**: Enhanced EM algorithm for probabilistic sequence clustering
- **Comprehensive Analysis**: Compare multiple distance and clustering methods
- **Optimized Performance**: Vectorized operations and efficient algorithms
- **Missing Value Support**: Robust handling of incomplete sequences

## Installation

```r
# Install from GitHub
devtools::install_github("mohsaqr/tnaExtras")
```

## Quick Start

```r
library(tnaExtras)

# Example sequence data
data <- data.frame(
  T1 = c("A", "B", "A", "C", "A", "B"),
  T2 = c("B", "A", "B", "A", "C", "A"),
  T3 = c("C", "C", "A", "B", "B", "C")
)

# Basic clustering with default settings (Euclidean distance + PAM)
result <- cluster_sequences(data, k = 2)
print(result)

# Mixture Markov Model clustering
mmm_result <- cluster_mmm(data, k = 2, n_starts = 5)
print(mmm_result)

# Compare multiple clustering methods
comparison <- compare_clustering_methods(data, k = 2)
print(comparison)
```

## Distance Methods

### Traditional Methods
- **`euclidean`**: Fast numeric encoding (recommended for exploration)
- **`hamming`**: Position-wise differences
- **`lcs`**: Longest Common Subsequence similarity
- **`start_position`**: First occurrence positions
- **`transition`**: Advanced Markov transition analysis
- **`optimal_matching`**: Optimal matching with dynamic programming

### StringDist Methods
- **`lv`**: Levenshtein distance
- **`jw`**: Jaro-Winkler distance (with prefix parameter)
- **`qgram`**: Q-gram distance (with q parameter)
- **`cosine`**: Cosine distance on q-grams
- **`jaccard`**: Jaccard distance on q-grams
- And more...

## Clustering Methods

- **`pam`**: Partitioning Around Medoids (default)
- **`ward.D2`**: Ward's method (recommended for hierarchical)
- **`complete`**: Complete linkage
- **`average`**: Average linkage (UPGMA)
- **`single`**: Single linkage

## Advanced Features

### Comprehensive Analysis

```r
# Test multiple methods and parameters
analysis <- cluster_complete_analysis(
  data, 
  k_range = 2:4,
  distance_methods = c("euclidean", "lcs", "transition"),
  clustering_methods = c("pam", "ward.D2")
)

# View summary of all combinations
print(analysis$summary)
```

### Distance Method Exploration

```r
# Compare different distance measures
distances <- analyze_distances(data, 
  methods = c("euclidean", "lv", "transition"))

# Use distances for further analysis
result_euclidean <- cluster_sequences(data, k = 2, "euclidean")
result_levenshtein <- cluster_sequences(data, k = 2, "lv")
```

### Advanced Transition Analysis

```r
# Fine-tune transition distance weights
result <- cluster_sequences(data, k = 2, 
  distance_method = "transition",
  first_order_weight = 0.4,
  second_order_weight = 0.3,
  time_weighted_weight = 0.2,
  complexity_weight = 0.1,
  persistence_weight = 0.1
)
```

## Examples

### Basic Workflow

```r
# 1. Load data
data <- your_sequence_data

# 2. Explore distance methods
distances <- analyze_distances(data, 
  methods = c("euclidean", "lcs", "transition"))

# 3. Compare clustering approaches
comparison <- compare_clustering_methods(data, k = 3)
print(comparison[1:3, ])  # Show top 3 methods

# 4. Apply best method
best_method <- comparison$method[1]
result <- cluster_sequences(data, k = 3, clustering_method = best_method)

# 5. Examine results
print(result$assignments)
print(result$silhouette)
```

### Mixture Markov Model Analysis

```r
# Fit MMM with multiple restarts for robustness
mmm_result <- cluster_mmm(data, 
  k = 3, 
  n_starts = 20,
  max_iter = 500,
  seed = 123
)

# Examine cluster quality
print(mmm_result$avg_posterior_per_cluster)
print(mmm_result$entropy)

# Access Markov models for each cluster
for (i in 1:3) {
  cat("Cluster", i, "transition matrix:\n")
  print(mmm_result$models[[i]]$transition)
  cat("\n")
}
```

## Performance Tips

1. **Start with Euclidean distance** for quick exploration
2. **Use multiple random starts** for MMM clustering on complex data
3. **Compare methods systematically** using `compare_clustering_methods()`
4. **Consider data size** when choosing between methods:
   - Small datasets (< 100 sequences): Any method
   - Medium datasets (100-1000): Euclidean, Hamming, LCS
   - Large datasets (> 1000): Euclidean, Hamming

## Citation

If you use tnaCluster in your research, please cite:

```
Author(s) (Year). tnaCluster: Advanced Sequence Analysis and Clustering 
for Temporal Network Analysis. R package version X.X.X.
```

## License

This package is licensed under GPL-3. See LICENSE file for details.

## Contributing

We welcome contributions! Please see our contributing guidelines and submit issues or pull requests on GitHub.

## Support

- Documentation: `?cluster_sequences`, `?cluster_mmm`, etc.
- Issues: [GitHub Issues](https://github.com/username/tnaCluster/issues)
- Examples: See package vignettes 
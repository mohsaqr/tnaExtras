# tnaExtras: Transition Network Analysis and Sequential Pattern Detection Extras

A comprehensive R package for analyzing temporal patterns in sequential data. This package provides tools for sequence comparison, pattern analysis, and computation of sequence-level indices with focus on complexity systems, Markov chain dynamics, and temporal network analysis.

## Features

- **Pattern Analysis**: Comprehensive analysis of pattern differences between groups using support, lift, confidence, and effect size measures
- **Sequence Comparison**: Advanced sequence comparison with statistical testing capabilities (chi-square, Fisher's exact tests)
- **Sequence Indices**: Computation of sequence-level statistics and complexity measures including entropy, diversity, and system dynamics
- **Statistical Testing**: Built-in statistical testing with multiple comparison corrections
- **Visualization**: Heatmap visualizations for pattern analysis results

## Installation

You can install the development version of tnaExtras from GitHub:

```r
# Install devtools if you haven't already
install.packages("devtools")

# Install tnaExtras from GitHub
devtools::install_github("mohsaqr/tnaExtras")
```

## Main Functions

### 1. `analyze_patterns()`

Comprehensive pattern analysis between two groups with multiple measures:

```r
library(tnaExtras)

# Create example data
data <- data.frame(
  T1 = c("plan", "discuss", "monitor", "plan", "discuss"),
  T2 = c("consensus", "emotion", "plan", "consensus", "emotion"),
  T3 = c("discuss", "plan", "consensus", "discuss", "plan"),
  Group = c("A", "B", "A", "B", "A")
)

# Run pattern analysis
result <- analyze_patterns(data, group_col = "Group")
print(result)
summary(result, measure = "support")
```

### 2. `compare_sequences()`

Advanced sequence comparison with statistical testing:

```r
# Basic discrimination analysis
result <- compare_sequences(data[1:3], data$Group, 
                          min_length = 2, max_length = 4)

# Statistical analysis with testing
result_stat <- compare_sequences(data[1:3], data$Group, 
                               statistical = TRUE,
                               correction = "bonferroni")
print(result_stat)
```

### 3. `compute_sequence_indices()`

Compute comprehensive sequence-level indices:

```r
# Define favorable states for your system
favorable_states <- c("consensus", "plan", "discuss")

# Compute indices
indices <- compute_sequence_indices(data, 
                                  group_col = "Group",
                                  favorable_states = favorable_states,
                                  return_group_summary = TRUE)

# Print summary
print_indices_summary(indices)
```

## Key Measures and Indices

### Pattern Analysis Measures

- **Support**: Frequency of patterns in each group
- **Lift**: How much more likely patterns are compared to chance
- **Confidence**: Conditional probability within groups
- **Effect Sizes**: Cohen's h, Cohen's d, Cramer's V, and other effect size measures

### Sequence Indices

- **Complexity Measures**: Entropy, diversity, transition complexity
- **System Dynamics**: Self-loop tendency, transition rates, feedback loops
- **Temporal Measures**: Initial state influence, attractor states, emergent patterns
- **Favorability Analysis**: Proportion of favorable states, stability measures

### Statistical Testing

- Chi-square tests for independence
- Fisher's exact tests for small samples
- Multiple comparison corrections (Bonferroni, Holm, BH, etc.)
- Automatic test selection based on expected frequencies

## Example Workflow

```r
library(tnaExtras)

# Load your sequential data
# Columns should represent time points, with an optional group column
data <- read.csv("your_sequence_data.csv")

# 1. Pattern Analysis
patterns <- analyze_patterns(data, 
                           group_col = "Group",
                           min_length = 2,
                           max_length = 5)

# 2. Statistical Sequence Comparison  
comparison <- compare_sequences(data, 
                             "Group",
                              statistical = TRUE)

# 3. Sequence Indices
indices <- compute_sequence_indices(data,
                                  group_col = "Group",
                                  return_group_summary = TRUE)

# View results
print(patterns)
print(comparison)
print_indices_summary(indices)
```

## Data Format

The package expects data in wide format where:
- Each row represents a sequence (individual, trial, etc.)
- Each column represents a time point
- An optional group column for between-group comparisons
- State values should be character strings

Example:
```
  T1        T2        T3        Group
  plan      discuss   consensus A
  monitor   emotion   plan      B
  discuss   plan      discuss   A
```

## Visualization

The package includes built-in visualization capabilities:
- Heatmaps for pattern discrimination
- Automatic color coding for over/under-representation
- Separate plots by subsequence length when requested

## Dependencies

- Base R (stats, graphics, grDevices, utils)
- Suggested: testthat, knitr, rmarkdown

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Citation

If you use this package in your research, please cite:

```
Author Name (2025). tnaExtras: Temporal Network Analysis and Sequential Pattern Detection Extras.
R package version 0.1.0. https://github.com/yourusername/tnaExtras
``` 

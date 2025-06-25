# tnaExtras: Transition Network Analysis and Sequential Pattern Detection Extras

A comprehensive R package for analyzing temporal patterns in sequential data. This package provides tools for sequence comparison, pattern analysis, and computation of sequence-level indices with focus on complexity systems, Markov chain dynamics, and temporal network analysis.

**NEW in v0.2.0: Multi-Group Support** - Now supports analysis across multiple groups (3+) in addition to the original two-group comparisons.

## Features

- **Pattern Analysis**: Comprehensive analysis of pattern differences between groups using support, lift, confidence, and effect size measures
  - **NEW**: Multi-group pattern analysis with `analyze_patterns_multi()`
- **Sequence Comparison**: Advanced sequence comparison with statistical testing capabilities (chi-square, Fisher's exact tests)
  - **NEW**: Multi-group sequence comparison with `compare_sequences_multi()`
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

### Multi-Group Analysis (NEW)

#### 1. `analyze_patterns_multi()`

Comprehensive pattern analysis across multiple groups:

```r
library(tnaExtras)

# Create example data with 3 groups
data <- data.frame(
  T1 = c("plan", "discuss", "monitor", "plan", "discuss", "evaluate"),
  T2 = c("consensus", "emotion", "plan", "consensus", "emotion", "plan"),
  T3 = c("discuss", "plan", "consensus", "discuss", "plan", "monitor"),
  Group = c("A", "B", "C", "A", "B", "C")
)

# Run multi-group pattern analysis
result <- analyze_patterns_multi(data, group_col = "Group")
print(result)
summary(result, measure = "support")
```

#### 2. `compare_sequences_multi()`

Advanced sequence comparison across multiple groups:

```r
# Multi-group discrimination analysis
result <- compare_sequences_multi(data, "Group", 
                                 min_length = 2, max_length = 4)
print(result)
summary(result)
```

### Traditional Two-Group Analysis

The original functions remain available for backward compatibility and still provide the most comprehensive analysis for two-group comparisons:

#### 1. `analyze_patterns()`

Comprehensive pattern analysis between two groups with multiple measures:

```r
# For two-group data
data_2groups <- data.frame(
  T1 = c("plan", "discuss", "monitor", "plan", "discuss"),
  T2 = c("consensus", "emotion", "plan", "consensus", "emotion"),
  T3 = c("discuss", "plan", "consensus", "discuss", "plan"),
  Group = c("A", "B", "A", "B", "A")
)

# Run pattern analysis (automatically detects 2 groups)
result <- analyze_patterns(data_2groups, group_col = "Group")
print(result)
summary(result, measure = "support")
```

#### 2. `compare_sequences()`

Advanced sequence comparison with statistical testing (two groups):

```r
# Basic discrimination analysis
result <- compare_sequences(data_2groups[1:3], data_2groups$Group, 
                          min_length = 2, max_length = 4)

# Statistical analysis with testing
result_stat <- compare_sequences(data_2groups[1:3], data_2groups$Group, 
                               statistical = TRUE,
                               correction = "bonferroni")
print(result_stat)
```

### 3. `compute_sequence_indices()`

Compute comprehensive sequence-level indices (supports any number of groups):

```r
# Define favorable states for your system
favorable_states <- c("consensus", "plan", "discuss")

# Compute indices (works with any number of groups)
indices <- compute_sequence_indices(data, 
                                  group_col = "Group",
                                  favorable_states = favorable_states,
                                  return_group_summary = TRUE)

# Print summary
print_indices_summary(indices)
```

## Multi-Group vs Two-Group Analysis

### When to Use Multi-Group Functions

- **3+ groups**: Use `analyze_patterns_multi()` and `compare_sequences_multi()`
- **Exploratory analysis**: Multi-group functions provide discrimination measures across all groups
- **Simple comparisons**: Multi-group functions are more straightforward for basic pattern detection

### When to Use Two-Group Functions

- **Exactly 2 groups**: Traditional functions provide more detailed statistical measures
- **Statistical testing**: Two-group functions include chi-square tests, Fisher's exact tests, effect sizes
- **Detailed analysis**: More comprehensive measures (lift, confidence, effect sizes) available
- **Research contexts**: When you need rigorous statistical comparisons

## Key Measures and Indices

### Multi-Group Pattern Analysis Measures

- **Support**: Frequency of patterns in each group
- **Support Range**: Difference between maximum and minimum support across groups
- **Support Variance**: Variability of pattern support across groups
- **Dominant Group**: Group with highest support for each pattern

### Two-Group Pattern Analysis Measures

- **Support**: Frequency of patterns in each group
- **Lift**: How much more likely patterns are compared to chance
- **Confidence**: Conditional probability within groups
- **Effect Sizes**: Cohen's h, Cohen's d, Cramer's V, and other effect size measures

### Sequence Indices (Universal)

- **Complexity Measures**: Entropy, diversity, transition complexity
- **System Dynamics**: Self-loop tendency, transition rates, feedback loops
- **Temporal Measures**: Initial state influence, attractor states, emergent patterns
- **Favorability Analysis**: Proportion of favorable states, stability measures

### Statistical Testing (Two-Group Only)

- Chi-square tests for independence
- Fisher's exact tests for small samples
- Multiple comparison corrections (Bonferroni, Holm, BH, etc.)
- Automatic test selection based on expected frequencies

## Example Workflows

### Multi-Group Workflow

```r
library(tnaExtras)

# Load your sequential data with multiple groups
data <- read.csv("your_sequence_data.csv")

# 1. Multi-Group Pattern Analysis
patterns <- analyze_patterns_multi(data, 
                                  group_col = "Group",
                                  min_length = 2,
                                  max_length = 5)

# 2. Multi-Group Sequence Comparison  
comparison <- compare_sequences_multi(data, 
                                     "Group",
                                     min_length = 2,
                                     max_length = 4)

# 3. Sequence Indices (works with any number of groups)
indices <- compute_sequence_indices(data,
                                  group_col = "Group",
                                  return_group_summary = TRUE)

# View results
print(patterns)
print(comparison)
print_indices_summary(indices)
```

### Two-Group Workflow

```r
library(tnaExtras)

# Load your sequential data with two groups
data <- read.csv("your_sequence_data.csv")

# 1. Comprehensive Pattern Analysis
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
  discuss   plan      discuss   C
  evaluate  plan      monitor   A
```

## Backward Compatibility

All existing functions remain unchanged and fully functional. The package automatically detects the number of groups:

- **2 groups**: Original functions work as before
- **3+ groups**: Original functions issue a warning and suggest multi-group alternatives
- **New functions**: Available for users who want explicit multi-group analysis

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
Mohammed Saqr (2025). tnaExtras: Transition Network Analysis and Sequential Pattern Detection Extras.
R package version 0.2.0. https://github.com/mohsaqr/tnaExtras
``` 

# tnaExtras: Transition Network Analysis and Sequential Pattern Detection Extras

A comprehensive R package for analyzing temporal patterns in sequential data. This package provides tools for sequence comparison, pattern analysis, and computation of sequence-level indices with focus on complexity systems, Markov chain dynamics, and temporal network analysis.

**NEW in v0.3.0: Association Rule Learning** - Complete implementation of Apriori and FP-Growth algorithms for discovering frequent patterns and association rules in transaction data.

**Enhanced in v0.2.0: Multi-Group Support** - Now supports analysis across multiple groups (3+) in addition to the original two-group comparisons.

**Enhanced in v0.3.0: Actual Group Names** - All analysis functions now display actual group names in results, tables, and graphs instead of generic "A", "B", "C" labels.

## Features

- **Association Rule Learning**: Complete implementation of Apriori and FP-Growth algorithms for mining frequent patterns and association rules
  - **NEW**: `apriori_rules()` and `fp_growth_rules()` functions
  - **NEW**: Comprehensive visualization suite with 5 different plot types
  - **NEW**: Advanced rule filtering, ranking, and analysis utilities
- **Pattern Analysis**: Comprehensive analysis of pattern differences between groups using support, lift, confidence, and effect size measures
  - **Enhanced**: Multi-group pattern analysis with `analyze_patterns_multi()`
- **Sequence Comparison**: Advanced sequence comparison with statistical testing capabilities (chi-square, Fisher's exact tests)
  - **Enhanced**: Multi-group sequence comparison with `compare_sequences_multi()`
- **Sequence Indices**: Computation of sequence-level statistics and complexity measures including entropy, diversity, and system dynamics
- **Statistical Testing**: Built-in statistical testing with multiple comparison corrections
- **Visualization**: Comprehensive visualization support including heatmaps, scatter plots, network diagrams, and quality metrics

## Installation

You can install the development version of tnaExtras from GitHub:

```r
# Install devtools if you haven't already
install.packages("devtools")

# Install tnaExtras from GitHub
devtools::install_github("mohsaqr/tnaExtras")
```

## Included Datasets

### Student Engagement Dataset (NEW in v0.3.0)

The package includes a comprehensive student engagement dataset (`engagement_data`) with:
- **1000 students** with sequential engagement patterns
- **25 time points** per student sequence  
- **3 engagement groups**: Low (260 students), Moderate (225 students), Engaged (515 students)
- **3 engagement states**: Active, Average, Disengaged
- Perfect for demonstrating all tnaExtras analysis capabilities

```r
# Load and explore the dataset
data(engagement_data)
summary(engagement_data)

# Quick analysis example
patterns <- analyze_patterns_multi(engagement_data, group_col = "Group")
rules <- apriori_rules(engagement_data[,2:26], min_support = 0.1)
```

## Main Functions

### Association Rule Learning (NEW in v0.3.0)

#### 1. `apriori_rules()`

Discovers association rules using the classic Apriori algorithm:

```r
library(tnaExtras)

# Learning activity sequences
transactions <- list(
  c("plan", "discuss", "execute", "reflect"),
  c("plan", "research", "analyze", "present"),
  c("discuss", "execute", "collaborate", "reflect"),
  c("plan", "discuss", "execute", "evaluate"),
  c("research", "analyze", "collaborate", "present")
)

# Run Apriori algorithm
rules <- apriori_rules(transactions, 
                      min_support = 0.2, 
                      min_confidence = 0.6,
                      min_lift = 1.0)

print(rules)
summary(rules)
```

#### 2. `fp_growth_rules()`

Efficient pattern mining using the FP-Growth algorithm:

```r
# Same learning activity data
rules_fp <- fp_growth_rules(transactions, 
                           min_support = 0.2, 
                           min_confidence = 0.6)

# Compare algorithms
compare_rule_algorithms(list(
  Apriori = rules,
  FP_Growth = rules_fp
))
```

#### 3. Association Rules Visualization

Multiple visualization options using base R graphics:

```r
# Scatter plot: Support vs Confidence, colored by Lift
plot_rules_scatter(rules)

# Network diagram showing item relationships
plot_rules_network(rules, layout = "force")

# Quality metrics distribution
plot_rules_quality(rules)

# Frequency bar chart
plot_itemset_frequency(rules)

# Matrix visualization
plot_rules_matrix(rules)

# Generic plot method
plot(rules, type = "scatter")
plot(rules, type = "network")
```

#### 4. Rule Analysis and Filtering

Advanced rule manipulation and analysis:

```r
# Filter rules by quality metrics
high_quality_rules <- filter_association_rules(rules, 
                                              min_confidence = 0.8,
                                              min_lift = 1.5)

# Rank rules by different metrics
top_support_rules <- rank_association_rules(rules, by = "support")
top_lift_rules <- rank_association_rules(rules, by = "lift")

# Extract rules containing specific learning actions
planning_rules <- extract_rules_by_item(rules, "plan")
discuss_antecedent_rules <- extract_rules_by_item(rules, "discuss", side = "antecedent")

# Find redundant rules
redundant <- find_redundant_rules(rules)
clean_rules <- rules$rules[!redundant, ]

# Calculate rule overlap/similarity
overlap_matrix <- calculate_rule_overlap(rules$rules, method = "jaccard")
```

#### 5. Data Format Support

Works with multiple input formats:

```r
# List format (each element is a transaction)
list_data <- list(c("A", "B"), c("B", "C"), c("A", "C"))

# Data frame format (transaction ID and item columns)
df_data <- data.frame(
  transaction = c(1, 1, 2, 2, 3, 3),
  item = c("A", "B", "B", "C", "A", "C")
)

# Binary matrix format
matrix_data <- matrix(c(1,1,0, 1,0,1, 0,1,1), nrow=3)
colnames(matrix_data) <- c("A", "B", "C")

# All formats work with both algorithms
rules_list <- apriori_rules(list_data)
rules_df <- fp_growth_rules(df_data)
rules_matrix <- apriori_rules(matrix_data)
```

#### 6. Export and Utility Functions

```r
# Export rules to different formats
export_association_rules(rules, "rules.csv", format = "csv")
export_association_rules(rules, "rules.json", format = "json")
export_association_rules(rules, "rules.txt", format = "txt")

# Calculate specific metrics
metrics <- calculate_rule_metrics(c("plan"), c("discuss"), transaction_matrix)
support <- calculate_itemset_support(c("plan", "discuss"), transaction_matrix)

# Convert rules back to transactions for further analysis
rule_transactions <- rules_to_transactions(rules$rules)
```

### Multi-Group Analysis (Enhanced in v0.2.0)

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
- Multiple comparison corrections (all R p.adjust methods: Holm, Hochberg, Hommel, Bonferroni, BH, BY, FDR, none)
- Automatic test selection based on expected frequencies

## Example Workflows

### Association Rule Learning Workflow

```r
library(tnaExtras)

# 1. Prepare your learning activity data
transactions <- list(
  c("plan", "discuss", "execute", "reflect"),
  c("plan", "research", "analyze", "present"),
  c("discuss", "execute", "collaborate", "reflect"),
  c("plan", "discuss", "execute", "evaluate"),
  c("research", "analyze", "collaborate", "present"),
  c("plan", "research", "execute"),
  c("discuss", "collaborate", "analyze", "present")
)

# 2. Mine association rules using Apriori
apriori_rules_result <- apriori_rules(transactions, 
                                     min_support = 0.2,
                                     min_confidence = 0.6,
                                     min_lift = 1.0)

# 3. Mine association rules using FP-Growth
fp_growth_rules_result <- fp_growth_rules(transactions, 
                                         min_support = 0.2,
                                         min_confidence = 0.6)

# 4. Compare algorithms
compare_rule_algorithms(list(
  Apriori = apriori_rules_result,
  FP_Growth = fp_growth_rules_result
))

# 5. Analyze and filter rules
high_quality_rules <- filter_association_rules(apriori_rules_result$rules,
                                              min_confidence = 0.8,
                                              min_lift = 1.5)

planning_rules <- extract_rules_by_item(apriori_rules_result$rules, "plan")

# 6. Visualize results
plot_rules_scatter(apriori_rules_result)
plot_rules_network(apriori_rules_result, top_n = 15)
plot_rules_quality(apriori_rules_result)

# 7. Export results
export_association_rules(apriori_rules_result, "learning_patterns.csv")
```

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

- **Base R only** (stats, graphics, grDevices, utils)
- **Zero external dependencies** for all core functionality including association rule learning
- Suggested: testthat, knitr, rmarkdown (for development and documentation only)

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Citation

If you use this package in your research, please cite:

```
Mohammed Saqr (2025). tnaExtras: Transition Network Analysis and Sequential Pattern Detection Extras.
R package version 0.3.0. https://github.com/mohsaqr/tnaExtras
``` 

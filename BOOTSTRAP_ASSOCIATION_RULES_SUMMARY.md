# Bootstrap Association Rules - Function Summary

## Overview

The `bootstrap_association_rules()` function performs **bootstrapped association rule mining** on transaction data to assess the **stability and reproducibility** of discovered rules. This approach helps identify rules that consistently appear across multiple random samples of the data, providing more robust and reliable patterns.

## Key Features

1. **Stability Assessment**: Identifies rules that appear consistently across bootstrap samples
2. **Confidence Intervals**: Provides statistical confidence intervals for rule metrics (support, confidence, lift)
3. **Algorithm Support**: Works with both Apriori and FP-Growth algorithms
4. **Parallel Processing**: Utilizes multi-core processing for faster computation
5. **Frequency Filtering**: Filters rules based on their occurrence frequency across samples

## Function Signature

```r
bootstrap_association_rules(
  transactions,
  algorithm = "apriori",
  n_reps = 100,
  min_frequency = 0.9,
  min_support = 0.1,
  min_confidence = 0.8,
  min_lift = 1.0,
  max_length = 10,
  parallel = TRUE,
  n_cores = NULL,
  verbose = TRUE
)
```

## Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `transactions` | list/matrix/data.frame | *required* | Transaction data in various formats |
| `algorithm` | character | `"apriori"` | Algorithm to use: `"apriori"` or `"fp_growth"` |
| `n_reps` | integer | `100` | Number of bootstrap replications |
| `min_frequency` | numeric | `0.9` | Minimum frequency of rule occurrence (0-1) to retain |
| `min_support` | numeric | `0.1` | Minimum support threshold for rules |
| `min_confidence` | numeric | `0.8` | Minimum confidence threshold for rules |
| `min_lift` | numeric | `1.0` | Minimum lift threshold for rules |
| `max_length` | integer | `10` | Maximum length of itemsets (Apriori only) |
| `parallel` | logical | `TRUE` | Whether to use parallel processing |
| `n_cores` | integer | `NULL` | Number of cores to use (defaults to `detectCores() - 1`) |
| `verbose` | logical | `TRUE` | Whether to show progress messages |

## Return Value

Returns an object of class `"bootstrapped_association_rules"` containing:

### Components

1. **`rules`**: Data frame with aggregated rule statistics
   - `antecedent`: Left-hand side of the rule
   - `consequent`: Right-hand side of the rule
   - `frequency`: Proportion of bootstrap samples containing the rule
   - `support_mean`: Mean support across samples
   - `support_lower`: Lower 95% confidence interval for support
   - `support_upper`: Upper 95% confidence interval for support
   - `confidence_mean`: Mean confidence across samples
   - `confidence_lower`: Lower 95% confidence interval for confidence
   - `confidence_upper`: Upper 95% confidence interval for confidence
   - `lift_mean`: Mean lift across samples
   - `lift_lower`: Lower 95% confidence interval for lift
   - `lift_upper`: Upper 95% confidence interval for lift

2. **`parameters`**: List of input parameters used

3. **`summary`**: Summary statistics
   - `n_transactions`: Number of transactions in original data
   - `n_successful_reps`: Number of successful bootstrap replications
   - `n_rules`: Number of stable rules found

## How It Works

1. **Bootstrap Sampling**: Creates `n_reps` bootstrap samples (sampling with replacement)
2. **Rule Mining**: Runs association rule mining on each sample using the specified algorithm
3. **Aggregation**: Combines rules from all samples and counts their frequency
4. **Filtering**: Retains only rules that appear in at least `min_frequency` proportion of samples
5. **Statistical Summary**: Calculates mean and confidence intervals for each rule's metrics

## Test Results (tna::group_regulation dataset)

### Dataset Information
- **Rows**: 2,000 sequences
- **Columns**: 26 time points
- **Items**: Group regulation activities (plan, discuss, consensus, monitor, etc.)

### Test Scenarios

#### Test 1: Apriori with Standard Parameters
- **Algorithm**: Apriori
- **Replications**: 50
- **Min Frequency**: 0.7 (70%)
- **Min Support**: 0.1
- **Min Confidence**: 0.6
- **Result**: **4,108 stable rules found**

#### Test 2: FP-Growth with Standard Parameters
- **Algorithm**: FP-Growth
- **Replications**: 50
- **Min Frequency**: 0.7 (70%)
- **Min Support**: 0.1
- **Min Confidence**: 0.6
- **Result**: **11,802 stable rules found**

#### Test 3: Strict Frequency Threshold
- **Algorithm**: Apriori
- **Min Frequency**: 0.9 (90%)
- **Other parameters**: Same as Test 1
- **Result**: **3,064 stable rules found**

#### Test 4: Lower Support Threshold
- **Algorithm**: Apriori
- **Min Support**: 0.05 (lower than standard)
- **Other parameters**: Same as Test 1
- **Result**: **4,349 stable rules found**

### Frequency Distribution Analysis

For Test 1 (Apriori, 70% frequency threshold):
- **70-80% frequency**: 639 rules
- **80-90% frequency**: 921 rules
- **90-100% frequency**: 2,548 rules (62% of all stable rules)

### Example Stable Rules (100% frequency)

Top rules that appeared in **all** bootstrap samples:

1. `adapt,cohesion => consensus,coregulate,synthesis`
   - Support: 0.449 [0.401 - 0.497]
   - Confidence: 1.000 [1.000 - 1.000]
   - Lift: 1.915 [1.708 - 2.122]

2. `adapt,cohesion => coregulate,discuss,synthesis`
   - Support: 0.449 [0.401 - 0.497]
   - Confidence: 1.000 [1.000 - 1.000]
   - Lift: 1.915 [1.708 - 2.122]

3. `adapt,cohesion => coregulate,emotion,synthesis`
   - Support: 0.449 [0.401 - 0.497]
   - Confidence: 1.000 [1.000 - 1.000]
   - Lift: 1.915 [1.708 - 2.122]

## Usage Examples

### Basic Usage

```r
library(tnaExtras)
library(tna)

# Load data
data(group_regulation, package = "tna")

# Run bootstrap with default parameters
results <- bootstrap_association_rules(
  transactions = group_regulation,
  algorithm = "apriori",
  n_reps = 100,
  min_frequency = 0.9,
  verbose = TRUE
)

# View results
print(results)
summary(results)
```

### Custom Parameters

```r
# More lenient frequency threshold, stricter confidence
results <- bootstrap_association_rules(
  transactions = group_regulation,
  algorithm = "fp_growth",
  n_reps = 50,
  min_frequency = 0.7,  # Rules must appear in 70% of samples
  min_support = 0.05,
  min_confidence = 0.8,
  min_lift = 1.5,
  parallel = TRUE,
  n_cores = 4
)
```

### Accessing Results

```r
# Get the stable rules
stable_rules <- results$rules

# Filter for very stable rules (>95% frequency)
very_stable <- stable_rules[stable_rules$frequency > 0.95, ]

# Get rules with high lift
high_lift <- stable_rules[stable_rules$lift_mean > 2.0, ]

# View confidence intervals
head(stable_rules[, c("antecedent", "consequent", 
                       "lift_mean", "lift_lower", "lift_upper")])
```

## Interpretation Guide

### Frequency
- **High frequency (>0.9)**: Rule is very stable and consistently appears
- **Medium frequency (0.7-0.9)**: Rule is moderately stable
- **Low frequency (<0.7)**: Rule may be sample-dependent

### Confidence Intervals
- **Narrow intervals**: Rule metrics are consistent across samples
- **Wide intervals**: Rule metrics vary significantly across samples
- Use intervals to assess reliability of rule strength

### Comparison with Single-Run Mining
- Bootstrap results are more robust than single-run mining
- Rules with high frequency are less likely to be spurious
- Confidence intervals provide uncertainty quantification

## Performance Considerations

1. **Computational Cost**: Bootstrap is computationally intensive
   - Use `parallel = TRUE` for faster execution
   - Reduce `n_reps` for quicker testing (but less reliable results)

2. **Memory Usage**: Large datasets with many rules can consume significant memory
   - Consider filtering with stricter thresholds
   - Process in batches if needed

3. **Recommended Settings**:
   - **Quick testing**: `n_reps = 50`, `parallel = TRUE`
   - **Production analysis**: `n_reps = 100-500`, `parallel = TRUE`
   - **High confidence**: `n_reps = 1000+`, `min_frequency = 0.95`

## Advantages

1. **Robustness**: Identifies stable patterns that aren't artifacts of sampling
2. **Uncertainty Quantification**: Provides confidence intervals for metrics
3. **Reproducibility**: Results are less sensitive to data variations
4. **Statistical Rigor**: Based on established bootstrap methodology

## Limitations

1. **Computational Cost**: Significantly slower than single-run mining
2. **Parameter Sensitivity**: Results depend on `min_frequency` threshold
3. **Large Rule Sets**: May produce many rules with lenient thresholds

## Best Practices

1. **Start Conservative**: Begin with high `min_frequency` (0.9) and adjust down
2. **Use Parallel Processing**: Always enable for datasets with >100 transactions
3. **Validate Results**: Check confidence intervals for rule reliability
4. **Compare Algorithms**: Try both Apriori and FP-Growth to see which works better
5. **Document Parameters**: Keep track of parameter settings for reproducibility

## Related Functions

- `apriori_rules()`: Single-run Apriori algorithm
- `fp_growth_rules()`: Single-run FP-Growth algorithm
- `compare_rule_algorithms()`: Compare results from different algorithms
- `filter_association_rules()`: Filter rules by various criteria

## References

- Efron, B., & Tibshirani, R. J. (1994). An Introduction to the Bootstrap. CRC press.
- Agrawal, R., & Srikant, R. (1994). Fast algorithms for mining association rules. VLDB.
- Han, J., Pei, J., & Yin, Y. (2000). Mining frequent patterns without candidate generation. SIGMOD.

---

**Package**: tnaExtras  
**Version**: 0.3.0  
**Last Updated**: 2025-12-06

# tnatest 0.1.0

## New Features

* Initial release of the tnatest package
* Added `analyze_patterns()` function for comprehensive pattern analysis between groups
* Added `compare_sequences()` function for advanced sequence comparison with statistical testing
* Added `compute_sequence_indices()` function for sequence-level statistics and complexity measures

## Functions Included

### Pattern Analysis
* `analyze_patterns()` - Main pattern analysis function with support, lift, confidence, and effect size measures
* Support for multiple pattern length ranges and frequency thresholds
* Print and summary methods for pattern analysis results

### Sequence Comparison  
* `compare_sequences()` - Statistical sequence comparison with chi-square and Fisher's exact tests
* Multiple comparison correction methods (Bonferroni, Holm, BH, etc.)
* Discrimination analysis for non-statistical comparisons
* Visualization support with heatmaps
* Print and summary methods for comparison results

### Sequence Indices
* `compute_sequence_indices()` - Comprehensive sequence-level statistics
* Complexity measures (entropy, diversity, transition complexity)
* System dynamics (self-loop tendency, transition rates, feedback loops)
* Temporal measures (initial state influence, attractor states, emergent patterns)
* Favorability analysis for user-defined favorable states
* Group-level summaries and individual sequence analysis
* `print_indices_summary()` function for formatted output

## Helper Functions

* Internal helper functions for data processing, n-gram extraction, and statistical computations
* Robust error handling and input validation
* Support for missing data and edge cases

## Testing

* Comprehensive test suite covering all main functions
* Edge case testing and performance validation
* Mathematical formula verification

## Documentation

* Complete roxygen2 documentation for all exported functions
* Usage examples and parameter descriptions
* Package-level documentation with workflow examples 
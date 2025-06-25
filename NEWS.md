# tnaExtras 0.2.0 (In Development)

## Major New Features

* **Multi-Group Comparisons:**
    * Added `compare_sequences_multiple()` function to perform sequence comparisons across three or more groups by conducting all unique pairwise comparisons.
    * Added `analyze_patterns_multiple()` function to perform pattern analysis (support, lift, confidence, effect sizes) across three or more groups by conducting all unique pairwise analyses.
* **S3 Methods for Multi-Group Results:**
    * Implemented `print()` and `summary()` methods for `compare_sequences_multiple_result` and `analyze_patterns_multiple_result` objects to provide structured and informative outputs.

## Refactoring & Enhancements

* **Modular Internal Logic:**
    * Refactored `compare_sequences()` and `analyze_patterns()` to serve as robust engines for pairwise comparisons. Their core computational logic has been extracted into internal helper functions.
    * This change makes the codebase more modular and maintainable.
* **Utility Functions:**
    * Introduced a new `R/utils.R` file containing helper functions for common tasks such as group input validation (`validate_group_input`), data splitting by group (`split_data_by_group`), generating group pairs (`generate_group_pairs`), sequence preprocessing (`preprocess_sequences_for_analysis`), and safer statistical testing (`safe_statistical_test`).
* **Improved Error Handling:** Enhanced error messages and validation across various functions to provide clearer feedback to users.
* **Dynamic Column Naming:** Pairwise comparison results now feature dynamic column names based on the actual group names being compared (e.g., `freq_GroupX`, `support_GroupY`), improving clarity in multi-group contexts.

---

# tnaExtras 0.1.0

## New Features

* Initial release of the tnaExtras package. Package name updated from `tnatest`.
* Added `analyze_patterns()` function for comprehensive pattern analysis between two groups.
* Added `compare_sequences()` function for advanced sequence comparison (two groups) with statistical testing.
* Added `compute_sequence_indices()` function for sequence-level statistics and complexity measures.

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
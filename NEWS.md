# tnaExtras 0.3.0

## Major New Features

### Association Rule Learning
* **NEW**: Complete association rule mining implementation with Apriori and FP-Growth algorithms
* **NEW**: `apriori_rules()` - Classic Apriori algorithm for frequent itemset mining and rule generation
* **NEW**: `fp_growth_rules()` - FP-Growth algorithm for efficient pattern mining
* **NEW**: `compare_rule_algorithms()` - Compare results from different algorithms
* **NEW**: Comprehensive visualization suite with 5 different plot types using base R graphics only

### Transaction Data Support
* **NEW**: `prepare_transactions()` - Convert various data formats (list, matrix, data.frame) to transaction format
* **NEW**: `create_transaction_matrix()` - Create binary transaction matrices for efficient processing
* **NEW**: Automatic format detection and validation

### Rule Analysis and Filtering
* **NEW**: `filter_association_rules()` - Advanced rule filtering by multiple criteria
* **NEW**: `rank_association_rules()` - Sort rules by different quality metrics
* **NEW**: `extract_rules_by_item()` - Extract rules containing specific items
* **NEW**: `find_redundant_rules()` - Identify and remove redundant rules
* **NEW**: `calculate_rule_overlap()` - Compute rule similarity matrices

### Visualization Functions
* **NEW**: `plot_rules_scatter()` - Scatter plots with support vs confidence, colored by lift
* **NEW**: `plot_rules_network()` - Network diagrams showing item relationships
* **NEW**: `plot_itemset_frequency()` - Bar charts of frequent itemsets
* **NEW**: `plot_rules_quality()` - Multi-panel quality metric distributions
* **NEW**: `plot_rules_matrix()` - Matrix visualization of antecedent-consequent relationships
* **NEW**: Generic `plot.association_rules()` method supporting multiple plot types

### Utility Functions
* **NEW**: `calculate_rule_metrics()` - Compute support, confidence, lift, and conviction
* **NEW**: `calculate_itemset_support()` - Calculate itemset frequency
* **NEW**: `rules_to_transactions()` - Convert rules back to transaction format
* **NEW**: `export_association_rules()` - Export rules to CSV, JSON, or text formats

### S3 Methods
* **NEW**: Complete S3 method suite for association_rules class
* **NEW**: `print.association_rules()` - Formatted rule output
* **NEW**: `summary.association_rules()` - Statistical summaries
* **NEW**: `head.association_rules()` - Display top rules

## Major Enhancements

### Actual Group Names in Analysis
* **ENHANCED**: All analysis functions now display actual group names throughout results
* **ENHANCED**: `analyze_patterns()` and `analyze_patterns_multi()` show actual group names in column headers (e.g., `support_Expert` instead of `support_A`)
* **ENHANCED**: `compare_sequences()` and `compare_sequences_multi()` display actual group names in tables and visualizations
* **ENHANCED**: `compute_sequence_indices()` preserves actual group names in summaries and outputs
* **ENHANCED**: Heatmaps and visualizations now use actual group names on axes and legends

## Demo and Documentation
* **NEW**: `demo_association_rules.R` - Comprehensive demo with 8 examples
* **NEW**: Market basket analysis examples
* **NEW**: Educational sequential pattern examples
* **NEW**: Parameter sensitivity analysis
* **NEW**: Complete documentation with examples for all functions

## Technical Implementation
* All algorithms implemented with **zero external dependencies** (base R only)
* Efficient memory usage and processing for large datasets
* Robust error handling and input validation
* Support for various data formats and edge cases

# tnaExtras 0.2.0

## Multi-Group Support

### New Multi-Group Functions
* **NEW**: `analyze_patterns_multi()` - Pattern analysis for 3+ groups with support range, variance, and dominant group metrics
* **NEW**: `compare_sequences_multi()` - Multi-group sequence comparison using coefficient of variation
* **NEW**: Enhanced `compute_sequence_indices()` - Now supports any number of groups

### Improved Original Functions
* Enhanced `analyze_patterns()` - Now detects 3+ groups and suggests multi-group alternatives
* Enhanced `compare_sequences()` - Warns when more than 2 groups detected
* Backward compatibility maintained for all existing functionality

### group_tna Object Support
* **NEW**: `is_group_tna()` - Check if object is a group_tna object
* **NEW**: `convert_group_tna()` - Convert group_tna objects to tnaExtras format
* **NEW**: Full integration with all main functions - automatic detection and conversion
* **NEW**: `create_mock_group_tna()` - Create test group_tna objects

### New S3 Methods
* `print.pattern_analysis_multi()` - Formatted output for multi-group pattern analysis
* `summary.pattern_analysis_multi()` - Statistical summaries for multi-group results
* `print.compare_sequences_multi()` - Formatted output for multi-group sequence comparison
* `summary.compare_sequences_multi()` - Statistical summaries for multi-group comparisons

### Enhanced Documentation
* Updated README with multi-group examples and workflows
* Clear guidance on when to use multi-group vs two-group functions
* Complete examples for both traditional and multi-group analysis

### Demo Files
* **NEW**: `demo_multi_group.R` - Comprehensive 4-group analysis example
* **NEW**: `demo_group_tna_support.R` - Complete group_tna functionality demonstration

### Testing
* **NEW**: 69 comprehensive tests for group_tna functionality
* **NEW**: Multi-group analysis tests with 3 and 4 group scenarios
* All existing tests maintained and passing

## Bug Fixes
* Fixed critical issue in `convert_group_tna()` function for proper matrix data handling
* Improved error messages for frequency threshold issues
* Enhanced parameter validation across all functions

# tnaExtras 0.1.0

## New Features

* Initial release of the tnaExtras package
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
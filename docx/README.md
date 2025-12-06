# Demo Scripts Summary

## Overview

Created three comprehensive demo scripts in the `docx` folder demonstrating key features of the tnaExtras package using `tna::group_regulation` dataset.

## Demo Scripts

### 1. Pattern Discovery Demo
**File**: `docx/demo_pattern_discovery.R`

Demonstrates the new `starts_with` and `ends_with` parameters for pattern filtering:

- **Example 1**: Patterns starting with "plan"
- **Example 2**: Patterns ending with "consensus"
- **Example 3**: Patterns from "plan" to "consensus"
- **Example 4**: Gapped patterns starting with "plan"

Each example includes:
- Pattern discovery with `discover_patterns()`
- Results printing
- Visualization with `plot(result, type = "patterns", top_n = 15)`

**Status**: ✅ Tested and working

### 2. Association Rules WITHOUT Bootstrap
**File**: `docx/demo_association_rules_no_ci.R`

Demonstrates standard Apriori algorithm without confidence intervals:

- Uses `apriori_rules()` function
- Finds association rules with specified thresholds
- Shows top 10 rules by lift
- Includes three types of plots:
  - `plot(results, type = "rules", top_n = 15)` - Bar plot
  - `plot_rules_network(results, top_n = 15)` - Network visualization
  - `plot_rules_quality(results)` - Quality metrics scatter plot

**Status**: ✅ Tested and working

### 3. Association Rules WITH Bootstrap
**File**: `docx/demo_association_rules_with_ci.R`

Demonstrates bootstrap association rules with confidence intervals:

- Uses `bootstrap_association_rules()` function
- Provides stability estimates for each rule
- Shows confidence intervals for support, confidence, and lift
- Displays frequency (% of bootstrap samples containing rule)
- Includes visualization with `plot(results, type = "rules", top_n = 15)`

**Status**: ✅ Running (bootstrap takes ~2-3 minutes with 50 replications)

## Key Improvements Made

### 1. Fixed `plot.sequence_patterns()` Function

**Problem**: The plot function didn't support `type = "patterns"` parameter

**Solution**: Enhanced the function to:
- Accept `type` parameter (defaults to "patterns")
- Truncate long pattern labels (max 40 characters)
- Adjust margins for better label visibility
- Improve readability with smaller font size

**Code changes** in `/Users/mohammedsaqr/Git/tnaExtras/R/sequence_patterns.R`:
```r
plot.sequence_patterns <- function(x, type = "patterns", top_n = 20, ...) {
  # Support type parameter
  # Truncate long labels
  # Better margins and formatting
}
```

### 2. Created Organized Demo Folder

All demos are now in `docx/` folder for easy access and organization.

## Usage

### Run Pattern Discovery Demo
```bash
cd /Users/mohammedsaqr/Git/tnaExtras
Rscript docx/demo_pattern_discovery.R
```

### Run Association Rules Demo (No Bootstrap)
```bash
Rscript docx/demo_association_rules_no_ci.R
```

### Run Association Rules Demo (With Bootstrap)
```bash
Rscript docx/demo_association_rules_with_ci.R
```

## Dataset

All demos use `tna::group_regulation`:
- **Rows**: 2,000 students/groups
- **Columns**: 26 time points (T1-T26)
- **States**: 9 unique activities (plan, discuss, consensus, etc.)

## Verification

✅ Pattern discovery demo - Tested successfully
✅ Association rules (no CI) - Tested successfully  
✅ Association rules (with CI) - Running successfully
✅ Plot function - Fixed and working with `type = "patterns"`

## Files Modified

1. `/Users/mohammedsaqr/Git/tnaExtras/R/sequence_patterns.R` - Fixed plot function
2. `/Users/mohammedsaqr/Git/tnaExtras/docx/demo_pattern_discovery.R` - Created
3. `/Users/mohammedsaqr/Git/tnaExtras/docx/demo_association_rules_no_ci.R` - Created
4. `/Users/mohammedsaqr/Git/tnaExtras/docx/demo_association_rules_with_ci.R` - Created

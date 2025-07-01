# Installation Instructions for tnatest Package

## Method 1: Install from Source (Current Method)

Since this package is not yet on CRAN or GitHub, follow these steps to install from source:

### Prerequisites

1. Make sure you have R installed (version 3.5.0 or higher recommended)
2. Install required packages:

```r
install.packages(c("devtools", "roxygen2", "testthat"))
```

### Installation Steps

1. **Option A: Using the build script**
   - Open R or RStudio
   - Set your working directory to the `tnaTest` folder
   - Run: `source("build_package.R")`
   - This will automatically build, check, and install the package

2. **Option B: Manual installation with devtools**
   ```r
   # Set working directory to the tnaTest folder
   setwd("path/to/tnaTest")
   
   # Install using devtools
   devtools::install()
   ```

3. **Option C: Build tarball and install**
   ```r
   # From the parent directory of tnaTest
   devtools::build("tnaTest")
   install.packages("tnatest_0.1.0.tar.gz", repos = NULL, type = "source")
   ```

### Verification

After installation, test that the package works:

```r
library(tnatest)

# Create test data
data <- data.frame(
  T1 = c("plan", "discuss", "monitor"),
  T2 = c("consensus", "emotion", "plan"),
  T3 = c("discuss", "plan", "consensus"),
  Group = c("A", "B", "A")
)

# Test main functions
result1 <- analyze_patterns(data, group_col = "Group")
result2 <- compare_sequences(data[1:3], data$Group)
result3 <- compute_sequence_indices(data, group_col = "Group")
```

## Method 2: Future GitHub Installation (Coming Soon)

Once the package is uploaded to GitHub, you will be able to install it with:

```r
# Install from GitHub (when available)
devtools::install_github("username/tnatest")
```

## Troubleshooting

### Common Issues

1. **"namespace 'xxx' is not available"**
   - Install missing dependencies: `install.packages(c("stats", "graphics", "grDevices", "utils"))`

2. **"could not find function 'xxx'"**
   - Make sure you loaded the package: `library(tnatest)`
   - Check that installation completed successfully

3. **Documentation not showing**
   - Run `devtools::document()` in the package directory
   - Reinstall the package

### Getting Help

- Check the package documentation: `help(package = "tnatest")`
- View function help: `?analyze_patterns`, `?compare_sequences`, `?compute_sequence_indices`
- Run the included tests: `devtools::test()`

## Package Structure

The package includes:
- **Main functions**: `analyze_patterns()`, `compare_sequences()`, `compute_sequence_indices()`
- **Helper functions**: Various internal functions for data processing
- **Methods**: Print and summary methods for result objects
- **Tests**: Comprehensive test suite in `tests/testthat/`
- **Documentation**: Complete function documentation with examples 
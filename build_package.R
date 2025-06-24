# Build and Install tnatest Package
# Run this script to build, check, and install the package

# Set working directory to package root
if (!basename(getwd()) == "tnaTest") {
  stop("Please run this script from the tnaTest directory")
}

# Load required packages
if (!require(devtools)) {
  install.packages("devtools")
  library(devtools)
}

if (!require(roxygen2)) {
  install.packages("roxygen2")
  library(roxygen2)
}

# Generate documentation
cat("Generating documentation...\n")
devtools::document()

# Check package
cat("Checking package...\n")
devtools::check()

# Build package
cat("Building package...\n")
devtools::build()

# Install package
cat("Installing package...\n")
devtools::install()

# Load and test
cat("Loading package for testing...\n")
library(tnatest)

# Run basic tests
cat("Running basic tests...\n")

# Test data
test_data <- data.frame(
  T1 = c("plan", "discuss", "monitor", "plan", "discuss"),
  T2 = c("consensus", "emotion", "plan", "consensus", "emotion"),
  T3 = c("discuss", "plan", "consensus", "discuss", "plan"),
  Group = c("A", "B", "A", "B", "A")
)

# Test analyze_patterns
cat("Testing analyze_patterns...\n")
tryCatch({
  result1 <- analyze_patterns(test_data, group_col = "Group")
  cat("✓ analyze_patterns works\n")
}, error = function(e) {
  cat("✗ analyze_patterns failed:", e$message, "\n")
})

# Test compare_sequences
cat("Testing compare_sequences...\n")
tryCatch({
  result2 <- compare_sequences(test_data[1:3], test_data$Group)
  cat("✓ compare_sequences works\n")
}, error = function(e) {
  cat("✗ compare_sequences failed:", e$message, "\n")
})

# Test compute_sequence_indices
cat("Testing compute_sequence_indices...\n")
tryCatch({
  result3 <- compute_sequence_indices(test_data, group_col = "Group")
  cat("✓ compute_sequence_indices works\n")
}, error = function(e) {
  cat("✗ compute_sequence_indices failed:", e$message, "\n")
})

cat("Package build and test complete!\n")
cat("Use 'library(tnatest)' to load the package.\n") 
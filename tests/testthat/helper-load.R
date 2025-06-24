# This helper file will attempt to load all functions from the R directory
# as if the package were loaded. This is a common approach for testing
# uninstalled packages.

# More robust way for devtools::test() or R CMD check:
# If functions are exported and package is "loaded" during tests,
# they should be available.
# However, for direct sourcing if not using testthat's full machinery:

# Get the list of .R files in the R directory
r_files <- list.files(path = "../../R", pattern = "\\.R$", full.names = TRUE)

# Source each file
for (r_file in r_files) {
  source(r_file)
}

# Alternatively, if you are using devtools for testing,
# devtools::load_all() might be called automatically or can be called here.
# For simplicity in this environment, explicit sourcing is used.
# Consider if devtools::load_all(quiet = TRUE) is better if devtools is available.
# try(devtools::load_all("../../.", quiet = TRUE), silent = TRUE)

# Load test data directly here if it's used across multiple test files for this context
# For now, test data is defined within each test file.

# Quick diagnostic to understand group_regulation structure
library(tna)
data(group_regulation, package = "tna")

cat("Dataset structure:\n")
str(group_regulation)

cat("\n\nFirst 3 rows:\n")
print(head(group_regulation, 3))

cat("\n\nColumn names:\n")
print(colnames(group_regulation))

cat("\n\nData class:\n")
print(class(group_regulation))

cat("\n\nIs it a data frame?\n")
print(is.data.frame(group_regulation))

cat("\n\nChecking what prepare_transactions does:\n")
# The issue is likely in how it detects the format
# For group_regulation, each ROW should be a transaction (student)
# Each COLUMN (T1, T2, ...) contains activities

# Let's manually create what we expect:
cat("\nManual transaction creation:\n")
manual_trans <- list()
for (i in 1:min(5, nrow(group_regulation))) {
    # Get all non-NA values from this row
    row_values <- as.character(group_regulation[i, ])
    row_values <- row_values[!is.na(row_values)]
    row_values <- unique(row_values)
    manual_trans[[i]] <- row_values
    cat(sprintf("Row %d: %s\n", i, paste(row_values, collapse = ", ")))
}

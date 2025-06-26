#!/usr/bin/env Rscript
# Prepare Engagement Dataset for tnaExtras Package
# ==================================================

# Read the raw CSV data
engagement_raw <- read.csv("data-raw/Grouped_engagement.csv", 
                          sep = ";", 
                          stringsAsFactors = FALSE,
                          na.strings = c("", "NA"))

# Clean and prepare the data
engagement_data <- engagement_raw

# Ensure consistent naming
colnames(engagement_data)[1] <- "Group"

# Convert Group to factor with meaningful levels
engagement_data$Group <- factor(engagement_data$Group, 
                               levels = c("Low", "Moderate", "Engaged"),
                               labels = c("Low", "Moderate", "Engaged"))

# Remove any completely empty rows
engagement_data <- engagement_data[!is.na(engagement_data$Group), ]

# Data summary
cat("Engagement Dataset Summary:\n")
cat("=========================\n")
cat("Total sequences:", nrow(engagement_data), "\n")
cat("Time points per sequence:", sum(grepl("^T", colnames(engagement_data))), "\n")
cat("Groups:\n")
print(table(engagement_data$Group))

cat("\nEngagement levels in data:\n")
engagement_levels <- unique(unlist(engagement_data[, grepl("^T", colnames(engagement_data))]))
engagement_levels <- engagement_levels[!is.na(engagement_levels)]
print(sort(engagement_levels))

cat("\nSample of data:\n")
print(head(engagement_data[, 1:6]))

# Save as package data
save(engagement_data, file = "data/engagement_data.rda", compress = "bzip2")

cat("\nDataset successfully prepared and saved to data/engagement_data.rda\n")
cat("This dataset contains student engagement sequences across 25 time points.\n")
cat("Perfect for tnaExtras sequential pattern analysis and association rule mining!\n") 
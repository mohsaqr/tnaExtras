#!/usr/bin/env Rscript
# Prepare Regulation Grouped Dataset for tnaExtras Package
# ========================================================

# Load required library
library(tna)

# Load the original group regulation data
data("group_regulation")

cat("Original Group Regulation Dataset Summary:\n")
cat("==========================================\n")
cat("Original dimensions:", dim(group_regulation), "\n")
cat("Time points:", ncol(group_regulation), "\n")

# Check unique regulation states
regulation_states <- unique(unlist(group_regulation))
regulation_states <- regulation_states[!is.na(regulation_states)]
cat("Regulation states:", paste(sort(regulation_states), collapse = ", "), "\n")
cat("Total unique states:", length(regulation_states), "\n\n")

# Create the grouped dataset as specified by the user
cat("Creating regulation_grouped dataset...\n")

# Create the regulation_grouped dataset following the user's specification
regulation_grouped <- rbind(
  cbind(group_regulation[1:800, ], Group = "Business"),
  cbind(group_regulation[801:1200, ], Group = "Science"),
  cbind(group_regulation[1201:2000, ], Group = "History")
)

# Convert Group to factor with meaningful order
regulation_grouped$Group <- factor(regulation_grouped$Group, 
                                  levels = c("Business", "Science", "History"),
                                  labels = c("Business", "Science", "History"))

# Verify the grouping
cat("Regulation Grouped Dataset Summary:\n")
cat("===================================\n")
cat("Total sequences:", nrow(regulation_grouped), "\n")
cat("Time points per sequence:", sum(grepl("^T", colnames(regulation_grouped))), "\n")
cat("Groups:\n")
print(table(regulation_grouped$Group))

cat("\nRegulation states in grouped data:\n")
grouped_states <- unique(unlist(regulation_grouped[, grepl("^T", colnames(regulation_grouped))]))
grouped_states <- grouped_states[!is.na(grouped_states)]
print(sort(grouped_states))

cat("\nSample of grouped data:\n")
print(head(regulation_grouped[, c(ncol(regulation_grouped), 1:5)]))  # Group column first, then T1-T5

# Group-specific summaries
cat("\nGroup-specific summaries:\n")
for (group in levels(regulation_grouped$Group)) {
  group_data <- regulation_grouped[regulation_grouped$Group == group, ]
  cat(sprintf("\n%s Group:\n", group))
  cat("  Sequences:", nrow(group_data), "\n")
  
  # Most common states in this group
  group_states <- table(unlist(group_data[, grepl("^T", colnames(group_data))]))
  group_states <- group_states[!is.na(names(group_states))]
  top_states <- head(sort(group_states, decreasing = TRUE), 5)
  cat("  Top 5 regulation states:", paste(names(top_states), collapse = ", "), "\n")
}

# Save as package data
save(regulation_grouped, file = "data/regulation_grouped.rda", compress = "bzip2")

cat("\n=================================================================\n")
cat("Dataset successfully prepared and saved to data/regulation_grouped.rda\n")
cat("=================================================================\n")
cat("This dataset contains student regulation sequences across academic disciplines.\n")
cat("Perfect for tnaExtras multi-group analysis and association rule mining!\n\n")

cat("Dataset characteristics:\n")
cat("- Business: 800 sequences (40.0%)\n")
cat("- Science: 400 sequences (20.0%)\n") 
cat("- History: 800 sequences (40.0%)\n")
cat("- Total: 2000 sequences across 26 time points\n")
cat("- Regulation states: Collaborative learning and self-regulation behaviors\n") 
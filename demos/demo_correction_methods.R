#!/usr/bin/env Rscript
# Demo: All Available Multiple Comparison Correction Methods in tnaExtras
# ======================================================================

library(tnaExtras)

cat("=== Multiple Comparison Correction Methods Demo ===\n\n")

# Create demo data with potential differences
set.seed(123)
demo_data <- data.frame(
  T1 = sample(c("plan", "discuss", "monitor", "evaluate"), 20, replace = TRUE),
  T2 = sample(c("discuss", "plan", "emotion", "consensus"), 20, replace = TRUE),
  T3 = sample(c("consensus", "evaluate", "plan", "monitor"), 20, replace = TRUE),
  Group = rep(c("Expert", "Novice"), each = 10)
)

cat("Demo dataset:")
print(head(demo_data))
cat("\n")

# Show all available correction methods
cat("Available correction methods in R:\n")
cat(paste(p.adjust.methods, collapse = ", "), "\n\n")

cat("Testing each correction method:\n")
cat("==============================\n\n")

# Test each correction method
for (method in p.adjust.methods) {
  cat("Method:", method, "\n")
  cat("Description:", switch(method,
    "holm" = "Holm's step-down method (1979)",
    "hochberg" = "Hochberg's step-up method (1988)", 
    "hommel" = "Hommel's method (1988)",
    "bonferroni" = "Bonferroni correction (most conservative)",
    "BH" = "Benjamini & Hochberg False Discovery Rate (1995)",
    "BY" = "Benjamini & Yekutieli FDR (2001)",
    "fdr" = "False Discovery Rate (alias for BH)",
    "none" = "No correction (raw p-values)"
  ), "\n")
  
  tryCatch({
    result <- compare_sequences(demo_data, "Group", 
                               statistical = TRUE, 
                               correction = method,
                               detailed = FALSE)
    
    # Extract p-values from the result - use length-specific results
    all_patterns <- do.call(rbind, lapply(names(result)[grepl("length_", names(result))], function(len) {
      if (!is.null(result[[len]]$patterns) && nrow(result[[len]]$patterns) > 0) {
        patterns_df <- result[[len]]$patterns
        if ("p_adjusted" %in% names(patterns_df)) {
          return(patterns_df[!is.na(patterns_df$p_value), ])
        }
      }
      return(data.frame())
    }))
    
    if (nrow(all_patterns) > 0) {
      n_significant <- sum(all_patterns$p_adjusted < 0.05, na.rm = TRUE)
      cat("  Result: Found", nrow(all_patterns), "testable patterns,", 
          n_significant, "significant after correction\n")
    } else {
      cat("  Result: No testable patterns found\n")
    }
    
  }, error = function(e) {
    cat("  ERROR:", e$message, "\n")
  })
  
  cat("\n")
}

cat("Key differences between methods:\n")
cat("================================\n")
cat("• Conservative (strict): Bonferroni, Holm\n")
cat("• Balanced: Hochberg, Hommel\n") 
cat("• Liberal (FDR control): BH, BY, fdr\n")
cat("• No correction: none\n\n")

cat("Recommendations:\n")
cat("• For exploratory analysis: BH or fdr\n")
cat("• For confirmatory studies: Holm or Bonferroni\n")
cat("• For balanced approach: Hochberg\n")
cat("• For no correction: none\n\n")

cat("Demo completed successfully!\n") 
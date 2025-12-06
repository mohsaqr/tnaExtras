# Quick bootstrap test with fixed data preparation
library(tnaExtras)
library(tna)

cat("Testing bootstrap with FIXED prepare_transactions\n")
cat("==================================================\n\n")

data(group_regulation, package = "tna")

# Quick test with fewer reps
set.seed(42)
results <- bootstrap_association_rules(
    transactions = group_regulation,
    algorithm = "apriori",
    n_reps = 20, # Fewer reps for quick test
    min_frequency = 0.7,
    min_support = 0.1,
    min_confidence = 0.6,
    min_lift = 1.0,
    max_length = 5,
    parallel = TRUE,
    verbose = TRUE
)

cat("\n\nResults:\n")
print(results)

if (!is.null(results) && nrow(results$rules) > 0) {
    cat("\n\nTop 10 rules with confidence intervals:\n")
    top10 <- head(results$rules, 10)
    for (i in 1:nrow(top10)) {
        cat(sprintf("\n%d. %s => %s\n", i, top10$antecedent[i], top10$consequent[i]))
        cat(sprintf(
            "   Frequency: %.2f (appeared in %.0f%% of samples)\n",
            top10$frequency[i], top10$frequency[i] * 100
        ))
        cat(sprintf(
            "   Support:    %.3f [%.3f - %.3f]\n",
            top10$support_mean[i], top10$support_lower[i], top10$support_upper[i]
        ))
        cat(sprintf(
            "   Confidence: %.3f [%.3f - %.3f]\n",
            top10$confidence_mean[i], top10$confidence_lower[i], top10$confidence_upper[i]
        ))
        cat(sprintf(
            "   Lift:       %.3f [%.3f - %.3f]\n",
            top10$lift_mean[i], top10$lift_lower[i], top10$lift_upper[i]
        ))
    }
}

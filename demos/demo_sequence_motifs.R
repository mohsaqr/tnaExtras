# ==============================================================================
# DEMO: Sequence Motif Analysis
# ==============================================================================
# This demo shows how to use find_patterns() for wildcard search
# and find_meta_paths() for type-level pattern discovery.
# ==============================================================================

# Load sample data
load("data/seqdata.rda")
seq_data <- seqdata[, -1]

cat("======================================================\n")
cat("SEQUENCE MOTIF ANALYSIS DEMO\n")
cat("======================================================\n\n")

# ------------------------------------------------------------------------------
# 1. WILDCARD PATTERN SEARCH
# ------------------------------------------------------------------------------
cat("=== 1. WILDCARD PATTERN SEARCH ===\n\n")

# Search for specific pattern with single wildcard (*)
result1 <- find_patterns(
  seq_data,
  pattern = c("plan", "*", "consensus"),  # Vector format
  verbose = FALSE
)
cat("Pattern: plan->*->consensus\n")
print(result1$patterns)
cat("\nInstances found:\n")
print(head(result1$instances, 10))

# Multi-wildcard pattern (**)
result2 <- find_patterns(
  seq_data,
  pattern = c("plan", "**", "plan"),  # plan returns to plan
  verbose = FALSE
)
cat("\n\nPattern: plan->**->plan (plan returns)\n")
print(result2$patterns)

# ------------------------------------------------------------------------------
# 2. AUTO-DISCOVER GAPPED PATTERNS
# ------------------------------------------------------------------------------
cat("\n=== 2. AUTO-DISCOVER GAPPED PATTERNS ===\n\n")

# Discover all gapped patterns automatically
gapped <- find_patterns(
  seq_data,
  pattern = NULL,  # Auto-discover
  min_gap = 1,
  max_gap = 2,
  min_support = 0.1,
  verbose = FALSE
)

cat("Gapped patterns found:", nrow(gapped$patterns), "\n\n")
print(head(gapped$patterns, 15))

# Filter by starting state
gapped_plan <- find_patterns(
  seq_data,
  pattern = NULL,
  min_gap = 1,
  max_gap = 2,
  start_state = "plan",
  verbose = FALSE
)
cat("\n\nGapped patterns starting with 'plan':", nrow(gapped_plan$patterns), "\n")
print(head(gapped_plan$patterns, 10))

# ------------------------------------------------------------------------------
# 3. META-PATH DISCOVERY
# ------------------------------------------------------------------------------
cat("\n=== 3. META-PATH DISCOVERY ===\n\n")

# Define node types
node_types <- list(
  cognitive = c("plan", "monitor", "adapt"),
  social = c("discuss", "consensus", "coregulate", "synthesis"),
  emotional = c("emotion", "cohesion")
)

# Auto-discover all meta-paths
meta <- find_meta_paths(
  seq_data,
  node_types = node_types,
  max_length = 3,
  min_support = 0.05
)

print(meta)
summary(meta)

# ------------------------------------------------------------------------------
# 4. SPECIFIC SCHEMA SEARCH
# ------------------------------------------------------------------------------
cat("\n=== 4. SPECIFIC SCHEMA SEARCH ===\n\n")

# Search for cognitive->social->cognitive pattern
meta_csc <- find_meta_paths(
  seq_data,
  node_types = node_types,
  schema = c("cognitive", "social", "cognitive"),  # Vector format
  verbose = FALSE
)

cat("Schema: cognitive->social->cognitive\n\n")
cat("Schema-level statistics:\n")
print(meta_csc$schemas)

cat("\nState-level patterns:\n")
print(head(meta_csc$patterns, 10))

# With wildcard
meta_return <- find_meta_paths(
  seq_data,
  node_types = node_types,
  schema = c("cognitive", "**", "cognitive"),  # Cognitive returns
  verbose = FALSE
)
cat("\n\nSchema: cognitive->**->cognitive (cognitive returns)\n")
print(meta_return$schemas)

# ------------------------------------------------------------------------------
# 5. FILTERING BY STATE AND SCHEMA
# ------------------------------------------------------------------------------
cat("\n=== 5. FILTERING BY STATE AND SCHEMA ===\n\n")

# Filter by starting state
meta_plan <- find_meta_paths(
  seq_data,
  node_types = node_types,
  max_length = 3,
  start_state = "plan",
  verbose = FALSE
)
cat("Patterns starting with 'plan':", nrow(meta_plan$patterns), "\n")
print(head(meta_plan$patterns[, c("pattern", "schema", "count", "support")], 10))

# Filter by starting schema type
meta_cog <- find_meta_paths(
  seq_data,
  node_types = node_types,
  max_length = 3,
  start_schema = "cognitive",
  verbose = FALSE
)
cat("\n\nPatterns with cognitive start:", nrow(meta_cog$patterns), "\n")
print(head(meta_cog$patterns[, c("pattern", "schema", "count", "support")], 10))

# Filter by ending state
meta_consensus <- find_meta_paths(
  seq_data,
  node_types = node_types,
  max_length = 3,
  end_state = "consensus",
  verbose = FALSE
)
cat("\n\nPatterns ending with 'consensus':", nrow(meta_consensus$patterns), "\n")
print(head(meta_consensus$patterns[, c("pattern", "schema", "count", "support")], 10))

# Combine filters
meta_combined <- find_meta_paths(
  seq_data,
  node_types = node_types,
  max_length = 3,
  start_state = "plan",
  end_schema = "social",
  verbose = FALSE
)
cat("\n\nPatterns: start='plan', end_schema='social':", nrow(meta_combined$patterns), "\n")
print(head(meta_combined$patterns[, c("pattern", "schema", "count", "support")], 10))

# ------------------------------------------------------------------------------
# 6. TYPE TRANSITIONS
# ------------------------------------------------------------------------------
cat("\n=== 6. TYPE TRANSITIONS ===\n\n")

cat("Type-to-type transition probabilities:\n")
print(meta$type_transitions)

cat("\n======================================================\n")
cat("DEMO COMPLETE\n")
cat("======================================================\n")

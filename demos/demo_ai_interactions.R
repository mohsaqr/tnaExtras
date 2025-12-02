# ==============================================================================
# AI-Student Interaction Sequence Analysis
# ==============================================================================
# Analysis of interaction patterns between AI roles and student intents
# using the unified tnaExtras pattern analysis functions.
# ==============================================================================

options(scipen = 999)  # Disable scientific notation

# ------------------------------------------------------------------------------
# DATA SETUP
# ------------------------------------------------------------------------------

# Load sequence data
seq_data <- preparedRC$sequence_data

# Define interaction types for meta-path analysis
node_types <- list(
  student_intents = c("Instruct", "Seek.Help", "Inquire", "Debug", "Socialize"),
  ai_roles = c("Tutor", "Analyzer", "Coder", "Supporter", "Sycophant")
)

# ==============================================================================
# 1. N-GRAM PATTERN ANALYSIS
# ==============================================================================
# Explore basic sequential patterns (3-5 states long)

patterns <- explore_patterns(
  seq_data,
  type = "ngrams",
  min_length = 3,
  max_length = 5,
  min_support = 0.05,
  correction = "fdr"
)

print(patterns)
summary(patterns)

# Extract only significant patterns
sig_ngrams <- significant_patterns(patterns)
print(sig_ngrams)

# ==============================================================================
# 2. ABSTRACT STRUCTURAL PATTERNS
# ==============================================================================
# Detect returns (A->*->A), repetitions (A->A), oscillations (A->B->A->B)

abstract <- explore_patterns(
  seq_data,
  type = "abstract",
  correction = "fdr"
)

print(abstract)

# ==============================================================================
# 3. GAP-CONSTRAINED PATTERNS
# ==============================================================================
# Find patterns with gaps (e.g., A->*->B, A->*->*->B)

gapped <- find_patterns(
  seq_data,
  pattern = NULL,  # Auto-discover
  max_gap = 2,
  min_support = 0.05,
  correction = "fdr"
)

print(gapped)

# Filter gapped patterns: student_intents to student_intents transitions
sig_student_gaps <- gapped$patterns |>
  dplyr::filter(significant == TRUE) |>
  dplyr::filter(grepl(paste(node_types$student_intents, collapse = "|"), pattern))

print(sig_student_gaps)

# ==============================================================================
# 4. META-PATH ANALYSIS (Type-Level Patterns)
# ==============================================================================
# Discover patterns across node types (student_intents <-> ai_roles)

# 4a. Auto-discover all meta-paths
meta <- find_meta_paths(
  seq_data,
  node_types = node_types,
  min_length = 2,
  max_length = 3,
  min_support = 0.05,
  correction = "fdr"
)

print(meta)
summary(meta)

# View type-to-type transition probabilities
print(meta$type_transitions)

# 4b. Longer meta-paths (3-4 types)
meta_long <- find_meta_paths(
  seq_data,
  node_types = node_types,
  min_length = 3,
  max_length = 4,
  min_support = 0.05,
  correction = "fdr"
)

print(meta_long)
summary(meta_long)

# ==============================================================================
# 5. SPECIFIC SCHEMA SEARCH
# ==============================================================================
# Search for student->AI->student interaction pattern

meta_sis <- find_meta_paths(
  seq_data,
  node_types = node_types,
  schema = c("student_intents", "ai_roles", "student_intents"),
  correction = "fdr"
)

print(meta_sis)
summary(meta_sis)

# View actual state-level instances
print(meta_sis$patterns)

# ==============================================================================
# 6. FILTERED ANALYSIS
# ==============================================================================

# Patterns starting with specific student intent
patterns_from_instruct <- find_meta_paths(
  seq_data,
  node_types = node_types,
  max_length = 3,
  start_state = "Instruct",
  correction = "fdr"
)

print(patterns_from_instruct)

# Patterns ending with specific AI role
patterns_to_tutor <- find_meta_paths(
  seq_data,
  node_types = node_types,
  max_length = 3,
  end_state = "Tutor",
  correction = "fdr"
)

print(patterns_to_tutor)

# ==============================================================================
# 7. VISUALIZATION
# ==============================================================================

# Plot top patterns
plot(patterns, top_n = 20)

# ==============================================================================
# 8. EXTRACT KEY FINDINGS
# ==============================================================================

# Significant n-gram patterns
sig_patterns <- significant_patterns(patterns)

# Significant meta-path schemas
sig_schemas <- meta$schemas[meta$schemas$significant, ]

# Significant state-level patterns from meta-paths
sig_meta_patterns <- meta$patterns[meta$patterns$significant, ]

# Summary of findings
cat("\n=== KEY FINDINGS ===\n")
cat("Significant n-gram patterns:", nrow(sig_patterns), "\n")
cat("Significant meta-path schemas:", nrow(sig_schemas), "\n")
cat("Significant state patterns:", nrow(sig_meta_patterns), "\n")


#' Student Self-Regulation Sequences by Academic Discipline Dataset
#'
#' A dataset containing sequential self-regulation patterns of 2000 students across 26 time points,
#' grouped by academic discipline. Students are from three different academic programs: Business,
#' Science, and History. This dataset is derived from the \code{group_regulation} dataset in the
#' \code{tna} package and is ideal for demonstrating multi-group analysis, sequential pattern mining,
#' and association rule learning in educational contexts.
#'
#' @format A data frame with 2000 rows and 27 columns:
#' \describe{
#'   \item{T1-T26}{Character variables representing self-regulation state at each time point:
#'     \itemize{
#'       \item "adapt" - Adapting strategies or approaches
#'       \item "cohesion" - Building group cohesion and collaboration
#'       \item "consensus" - Seeking consensus and agreement
#'       \item "coregulate" - Co-regulating with peers
#'       \item "discuss" - Engaging in discussion and dialogue
#'       \item "emotion" - Managing emotions and affect
#'       \item "monitor" - Monitoring progress and performance
#'       \item "plan" - Planning strategies and approaches
#'       \item "synthesis" - Synthesizing and integrating information
#'     }
#'   }
#'   \item{Group}{Factor with 3 levels: "Business", "Science", "History" - represents the student's academic discipline}
#' }
#'
#' @details
#' This dataset represents student self-regulation patterns collected across three academic disciplines
#' over 26 consecutive time points. Each row represents one student's regulation trajectory, and each 
#' time point (T1-T26) captures their self-regulation behavior at that moment.
#'
#' **Group Distribution:**
#' \itemize{
#'   \item Business: 800 students (40.0%)
#'   \item Science: 400 students (20.0%) 
#'   \item History: 800 students (40.0%)
#' }
#'
#' **Self-Regulation States:**
#' \itemize{
#'   \item \strong{adapt}: Adapting learning strategies based on feedback
#'   \item \strong{cohesion}: Building group cohesion and collaborative relationships
#'   \item \strong{consensus}: Seeking consensus and shared understanding
#'   \item \strong{coregulate}: Co-regulating learning with peers
#'   \item \strong{discuss}: Engaging in academic discussion and dialogue
#'   \item \strong{emotion}: Managing emotions and affective states
#'   \item \strong{monitor}: Monitoring learning progress and performance
#'   \item \strong{plan}: Planning learning strategies and approaches
#'   \item \strong{synthesis}: Synthesizing and integrating information
#' }
#'
#' **Most Common Regulation Patterns by Discipline:**
#' \itemize{
#'   \item \strong{Business}: consensus → plan → discuss → emotion → cohesion
#'   \item \strong{Science}: consensus → plan → discuss → emotion → coregulate
#'   \item \strong{History}: plan → consensus → discuss → emotion → coregulate
#' }
#'
#' **Potential Analyses:**
#' \itemize{
#'   \item Multi-group pattern analysis with \code{analyze_patterns_multi()}
#'   \item Academic discipline comparison with \code{compare_sequences_multi()}
#'   \item Self-regulation rule mining with \code{apriori_rules()} and \code{fp_growth_rules()}
#'   \item Sequence complexity analysis with \code{compute_sequence_indices()}
#'   \item Temporal regulation analysis and transition network modeling
#' }
#'
#' @examples
#' # Load the dataset
#' data(regulation_grouped)
#' 
#' # Basic summary
#' summary(regulation_grouped)
#' table(regulation_grouped$Group)
#' 
#' # Check regulation state distribution
#' regulation_summary <- table(unlist(regulation_grouped[, 1:26]))
#' print(regulation_summary)
#' 
#' # Example: Multi-group regulation analysis
#' \dontrun{
#' library(tnaExtras)
#' 
#' # Analyze regulation patterns across academic disciplines
#' patterns <- analyze_patterns_multi(regulation_grouped, 
#'                                   group_col = "Group",
#'                                   min_length = 2, 
#'                                   max_length = 4)
#' print(patterns)
#' 
#' # Compare regulation sequences between disciplines
#' comparison <- compare_sequences_multi(regulation_grouped, "Group")
#' print(comparison)
#' 
#' # Association rule mining for regulation transitions
#' # Convert to transaction format for rule mining
#' regulation_transactions <- lapply(1:nrow(regulation_grouped), function(i) {
#'   states <- as.character(regulation_grouped[i, 1:26])
#'   states[!is.na(states)]  # Remove NA values
#' })
#' 
#' # Mine frequent regulation patterns
#' rules <- apriori_rules(regulation_transactions, 
#'                       min_support = 0.05, 
#'                       min_confidence = 0.6)
#' print(rules)
#' 
#' # Visualize regulation pattern networks
#' plot_rules_network(rules, top_n = 20)
#' 
#' # Focus on specific discipline comparisons
#' business_science <- regulation_grouped[regulation_grouped$Group %in% c("Business", "Science"), ]
#' business_science$Group <- droplevels(business_science$Group)
#' 
#' # Detailed two-group analysis with statistical testing
#' detailed_analysis <- analyze_patterns(business_science, group_col = "Group", min_length = 2)
#' statistical_comparison <- compare_sequences(business_science[, 1:26], 
#'                                           business_science$Group, 
#'                                           statistical = TRUE, 
#'                                           correction = "BH")
#' 
#' # Compute regulation complexity indices
#' indices <- compute_sequence_indices(regulation_grouped,
#'                                   group_col = "Group", 
#'                                   favorable_states = c("plan", "monitor", "coregulate"),
#'                                   return_group_summary = TRUE)
#' print_indices_summary(indices)
#' 
#' # Discipline-specific regulation analysis
#' business_students <- regulation_grouped[regulation_grouped$Group == "Business", ]
#' business_patterns <- analyze_patterns_multi(business_students, 
#'                                           group_col = "Group", 
#'                                           min_length = 3)
#' 
#' # Extract specific regulation transitions
#' planning_rules <- extract_rules_by_item(rules, "plan")
#' consensus_rules <- extract_rules_by_item(rules, "consensus", side = "consequent")
#' }
#'
#' # Quick discipline comparison
#' by(regulation_grouped[, 1:5], regulation_grouped$Group, summary)
#'
#' @source Derived from the group_regulation dataset in the tna package.
#'   Original data represents student self-regulation behaviors in collaborative
#'   learning environments across different academic disciplines.
#'
#' @references
#' Saqr, M., & López-Pernas, S. (2021). The longitudinal association between
#' engagement and achievement varies by time, students' profiles, and achievement
#' state: A full program study. Computers & Education, 199, 104787.
#'
#' @seealso \code{\link{engagement_data}} for student engagement sequences,
#'   \code{\link[tna]{group_regulation}} for the original dataset
#'
#' @keywords datasets education regulation sequential temporal self-regulation
"regulation_grouped" 
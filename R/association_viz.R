# ==============================================================================
# ASSOCIATION RULES VISUALIZATION
# ==============================================================================

#' Plot Association Rules Scatter Plot
#'
#' Creates a scatter plot of association rules with support vs confidence,
#' colored by lift values. Uses only base R graphics.
#'
#' @param rules Association rules object or data frame
#' @param x_metric X-axis metric ("support", "confidence", "lift", "conviction")
#' @param y_metric Y-axis metric ("support", "confidence", "lift", "conviction") 
#' @param color_metric Metric to use for coloring points ("lift", "conviction", "count")
#' @param top_n Number of top rules to display (NULL for all)
#' @param main Plot title
#' @param ... Additional arguments passed to plot()
#' @export
#' @examples
#' # Create sample rules
#' transactions <- list(
#'   c("bread", "milk", "eggs"),
#'   c("bread", "butter"),
#'   c("milk", "eggs", "cheese"),
#'   c("bread", "milk", "butter")
#' )
#' rules <- apriori_rules(transactions, min_support = 0.25)
#' plot_rules_scatter(rules)
plot_rules_scatter <- function(rules, x_metric = "support", y_metric = "confidence", 
                              color_metric = "lift", top_n = NULL, 
                              main = "Association Rules", ...) {
  
  # Handle association_rules object
  if (inherits(rules, "association_rules")) {
    rules_df <- rules$rules
  } else if (is.data.frame(rules)) {
    rules_df <- rules
  } else {
    stop("rules must be an association_rules object or data frame")
  }
  
  if (nrow(rules_df) == 0) {
    plot(1, 1, type = "n", main = "No rules to plot")
    text(1, 1, "No association rules found", cex = 1.2)
    return(invisible(NULL))
  }
  
  # Validate metrics
  available_metrics <- names(rules_df)
  if (!x_metric %in% available_metrics) {
    stop("x_metric '", x_metric, "' not found. Available: ", paste(available_metrics, collapse = ", "))
  }
  if (!y_metric %in% available_metrics) {
    stop("y_metric '", y_metric, "' not found. Available: ", paste(available_metrics, collapse = ", "))
  }
  if (!color_metric %in% available_metrics) {
    stop("color_metric '", color_metric, "' not found. Available: ", paste(available_metrics, collapse = ", "))
  }
  
  # Select top rules if specified
  if (!is.null(top_n) && top_n < nrow(rules_df)) {
    rules_df <- head(rules_df[order(rules_df[[color_metric]], decreasing = TRUE), ], top_n)
  }
  
  # Extract values
  x_vals <- rules_df[[x_metric]]
  y_vals <- rules_df[[y_metric]]
  color_vals <- rules_df[[color_metric]]
  
  # Create color palette
  n_colors <- 100
  color_palette <- colorRampPalette(c("blue", "green", "yellow", "red"))(n_colors)
  
  # Map color values to palette
  color_breaks <- seq(min(color_vals), max(color_vals), length.out = n_colors + 1)
  color_indices <- cut(color_vals, breaks = color_breaks, labels = FALSE, include.lowest = TRUE)
  point_colors <- color_palette[color_indices]
  
  # Create scatter plot
  plot(x_vals, y_vals, 
       col = point_colors,
       pch = 19,
       cex = 1.2,
       xlab = paste(toupper(substring(x_metric, 1, 1)), substring(x_metric, 2), sep = ""),
       ylab = paste(toupper(substring(y_metric, 1, 1)), substring(y_metric, 2), sep = ""),
       main = main,
       ...)
  
  # Add grid
  grid(col = "lightgray", lty = "dotted")
  
  # Add color legend
  legend_x <- par("usr")[2] - 0.25 * (par("usr")[2] - par("usr")[1])
  legend_y_top <- par("usr")[4] - 0.1 * (par("usr")[4] - par("usr")[3])
  legend_y_bottom <- par("usr")[4] - 0.4 * (par("usr")[4] - par("usr")[3])
  
  # Legend colors
  legend_colors <- color_palette[seq(1, n_colors, length.out = 5)]
  legend_values <- round(seq(min(color_vals), max(color_vals), length.out = 5), 3)
  
  legend("topright", 
         legend = legend_values,
         fill = legend_colors,
         title = paste(toupper(substring(color_metric, 1, 1)), substring(color_metric, 2), sep = ""),
         cex = 0.8)
  
  # Add information
  mtext(paste("Rules:", nrow(rules_df)), side = 3, adj = 0, cex = 0.8)
  
  return(invisible(rules_df))
}

#' Plot Rules Network Diagram
#'
#' Creates a simple network visualization showing relationships between items.
#' Uses base R graphics to show nodes (items) and edges (rules).
#'
#' @param rules Association rules object or data frame
#' @param top_n Number of top rules to visualize
#' @param layout Layout method ("circle", "random", "force")
#' @param node_size Size scaling factor for nodes
#' @param edge_width Width scaling factor for edges
#' @param main Plot title
#' @export
plot_rules_network <- function(rules, top_n = 20, layout = "circle", 
                              node_size = 1, edge_width = 1, 
                              main = "Association Rules Network") {
  
  # Handle association_rules object
  if (inherits(rules, "association_rules")) {
    rules_df <- rules$rules
  } else if (is.data.frame(rules)) {
    rules_df <- rules
  } else {
    stop("rules must be an association_rules object or data frame")
  }
  
  if (nrow(rules_df) == 0) {
    plot(1, 1, type = "n", main = "No rules to plot")
    text(1, 1, "No association rules found", cex = 1.2)
    return(invisible(NULL))
  }
  
  # Select top rules
  if (top_n < nrow(rules_df)) {
    rules_df <- head(rules_df[order(rules_df$lift, decreasing = TRUE), ], top_n)
  }
  
  # Extract all unique items
  all_antecedents <- unlist(strsplit(rules_df$antecedent, ","))
  all_consequents <- unlist(strsplit(rules_df$consequent, ","))
  all_items <- unique(c(all_antecedents, all_consequents))
  n_items <- length(all_items)
  
  if (n_items == 0) {
    plot(1, 1, type = "n", main = "No items to plot")
    return(invisible(NULL))
  }
  
  # Calculate item frequencies
  item_freq <- table(c(all_antecedents, all_consequents))
  
  # Generate node positions
  if (layout == "circle") {
    angles <- seq(0, 2 * pi, length.out = n_items + 1)[1:n_items]
    x_pos <- cos(angles)
    y_pos <- sin(angles)
  } else if (layout == "random") {
    set.seed(123)  # For reproducibility
    x_pos <- runif(n_items, -1, 1)
    y_pos <- runif(n_items, -1, 1)
  } else if (layout == "force") {
    # Simple force-directed layout
    x_pos <- runif(n_items, -1, 1)
    y_pos <- runif(n_items, -1, 1)
    
    # Simple force simulation (basic implementation)
    for (iter in 1:50) {
      for (i in 1:n_items) {
        fx <- fy <- 0
        
        # Repulsive forces
        for (j in 1:n_items) {
          if (i != j) {
            dx <- x_pos[i] - x_pos[j]
            dy <- y_pos[i] - y_pos[j]
            dist <- sqrt(dx^2 + dy^2) + 0.01
            force <- 0.1 / dist^2
            fx <- fx + force * dx / dist
            fy <- fy + force * dy / dist
          }
        }
        
        # Update positions
        x_pos[i] <- x_pos[i] + fx * 0.1
        y_pos[i] <- y_pos[i] + fy * 0.1
        
        # Keep in bounds
        x_pos[i] <- max(-1, min(1, x_pos[i]))
        y_pos[i] <- max(-1, min(1, y_pos[i]))
      }
    }
  } else {
    stop("layout must be 'circle', 'random', or 'force'")
  }
  
  # Create plot
  plot(x_pos, y_pos, type = "n", 
       xlim = c(-1.3, 1.3), ylim = c(-1.3, 1.3),
       xlab = "", ylab = "", main = main,
       axes = FALSE)
  
  # Draw edges (rules)
  for (i in 1:nrow(rules_df)) {
    antecedent_items <- strsplit(rules_df$antecedent[i], ",")[[1]]
    consequent_items <- strsplit(rules_df$consequent[i], ",")[[1]]
    
    # Draw edges from each antecedent to each consequent
    for (ant in antecedent_items) {
      for (con in consequent_items) {
        ant_idx <- which(all_items == ant)
        con_idx <- which(all_items == con)
        
        if (length(ant_idx) > 0 && length(con_idx) > 0) {
          # Edge color based on lift
          lift_val <- rules_df$lift[i]
          edge_color <- if (lift_val > 2) "red" else if (lift_val > 1.5) "orange" else "gray"
          
          # Edge width based on confidence
          conf_val <- rules_df$confidence[i]
          edge_lwd <- max(1, conf_val * 3 * edge_width)
          
          # Draw arrow
          arrows(x_pos[ant_idx], y_pos[ant_idx], 
                x_pos[con_idx], y_pos[con_idx],
                col = edge_color, lwd = edge_lwd, length = 0.1)
        }
      }
    }
  }
  
  # Draw nodes
  node_colors <- rainbow(n_items, alpha = 0.7)
  node_sizes <- sqrt(as.numeric(item_freq[all_items])) * node_size + 1
  
  points(x_pos, y_pos, pch = 19, col = node_colors, cex = node_sizes)
  points(x_pos, y_pos, pch = 1, col = "black", cex = node_sizes, lwd = 1.5)
  
  # Add labels
  text(x_pos, y_pos, all_items, cex = 0.7, col = "black", font = 2)
  
  # Add legend
  legend("bottomright", 
         legend = c("High Lift (>2)", "Medium Lift (1.5-2)", "Low Lift (<1.5)"),
         col = c("red", "orange", "gray"), 
         lty = 1, lwd = 2, cex = 0.8,
         title = "Rule Strength")
  
  return(invisible(list(items = all_items, x = x_pos, y = y_pos)))
}

#' Plot Itemset Frequency Bar Chart
#'
#' Creates a bar chart showing the frequency of frequent itemsets.
#'
#' @param rules Association rules object or data frame
#' @param top_n Number of top itemsets to show
#' @param type Type of itemsets to show ("antecedent", "consequent", "both")
#' @param main Plot title
#' @export
plot_itemset_frequency <- function(rules, top_n = 20, type = "both", 
                                  main = "Itemset Frequencies") {
  
  # Handle association_rules object
  if (inherits(rules, "association_rules")) {
    rules_df <- rules$rules
  } else if (is.data.frame(rules)) {
    rules_df <- rules
  } else {
    stop("rules must be an association_rules object or data frame")
  }
  
  if (nrow(rules_df) == 0) {
    plot(1, 1, type = "n", main = "No rules to plot")
    text(1, 1, "No association rules found", cex = 1.2)
    return(invisible(NULL))
  }
  
  # Extract itemsets
  itemsets <- character(0)
  
  if (type %in% c("antecedent", "both")) {
    itemsets <- c(itemsets, rules_df$antecedent)
  }
  
  if (type %in% c("consequent", "both")) {
    itemsets <- c(itemsets, rules_df$consequent)
  }
  
  # Count frequencies
  itemset_freq <- table(itemsets)
  itemset_freq <- sort(itemset_freq, decreasing = TRUE)
  
  # Select top itemsets
  if (top_n < length(itemset_freq)) {
    itemset_freq <- itemset_freq[1:top_n]
  }
  
  # Create bar plot
  par(mar = c(8, 4, 4, 2))  # Increase bottom margin for labels
  
  bp <- barplot(itemset_freq, 
                main = main,
                ylab = "Frequency",
                col = rainbow(length(itemset_freq), alpha = 0.7),
                border = "black",
                las = 2,  # Rotate labels
                cex.names = 0.8)
  
  # Add value labels on bars
  text(bp, itemset_freq + max(itemset_freq) * 0.01, 
       labels = itemset_freq, cex = 0.8, pos = 3)
  
  # Add grid
  grid(NA, NULL, col = "lightgray", lty = "dotted")
  
  # Redraw bars over grid
  barplot(itemset_freq, 
          col = rainbow(length(itemset_freq), alpha = 0.7),
          border = "black",
          las = 2,
          cex.names = 0.8,
          add = TRUE)
  
  par(mar = c(5, 4, 4, 2))  # Reset margins
  
  return(invisible(itemset_freq))
}

#' Plot Rules Quality Metrics
#'
#' Creates a multi-panel plot showing distributions of rule quality metrics.
#'
#' @param rules Association rules object or data frame
#' @param metrics Vector of metrics to plot
#' @param main Main title for the plot
#' @export
plot_rules_quality <- function(rules, metrics = c("support", "confidence", "lift"), 
                              main = "Rule Quality Metrics") {
  
  # Handle association_rules object
  if (inherits(rules, "association_rules")) {
    rules_df <- rules$rules
  } else if (is.data.frame(rules)) {
    rules_df <- rules
  } else {
    stop("rules must be an association_rules object or data frame")
  }
  
  if (nrow(rules_df) == 0) {
    plot(1, 1, type = "n", main = "No rules to plot")
    text(1, 1, "No association rules found", cex = 1.2)
    return(invisible(NULL))
  }
  
  # Validate metrics
  available_metrics <- names(rules_df)
  metrics <- metrics[metrics %in% available_metrics]
  
  if (length(metrics) == 0) {
    stop("No valid metrics found. Available: ", paste(available_metrics, collapse = ", "))
  }
  
  # Set up multi-panel plot
  n_metrics <- length(metrics)
  n_cols <- ceiling(sqrt(n_metrics))
  n_rows <- ceiling(n_metrics / n_cols)
  
  par(mfrow = c(n_rows, n_cols))
  
  # Create histogram for each metric
  for (metric in metrics) {
    values <- rules_df[[metric]]
    
    hist(values, 
         main = paste(toupper(substring(metric, 1, 1)), substring(metric, 2), sep = ""),
         xlab = metric,
         ylab = "Frequency",
         col = "lightblue",
         border = "black",
         breaks = min(20, max(5, length(unique(values)))),
         las = 1)
    
    # Add vertical line for mean
    abline(v = mean(values), col = "red", lwd = 2, lty = 2)
    
    # Add statistics
    text(par("usr")[1] + 0.7 * (par("usr")[2] - par("usr")[1]),
         par("usr")[4] - 0.1 * (par("usr")[4] - par("usr")[3]),
         paste("Mean:", round(mean(values), 3)),
         cex = 0.8)
    
    text(par("usr")[1] + 0.7 * (par("usr")[2] - par("usr")[1]),
         par("usr")[4] - 0.2 * (par("usr")[4] - par("usr")[3]),
         paste("Median:", round(median(values), 3)),
         cex = 0.8)
  }
  
  # Reset layout
  par(mfrow = c(1, 1))
  
  # Add overall title
  mtext(main, outer = TRUE, cex = 1.2, font = 2, line = -2)
  
  return(invisible(rules_df))
}

#' Plot Rules Matrix
#'
#' Creates a matrix visualization of rules showing antecedent-consequent relationships.
#'
#' @param rules Association rules object or data frame
#' @param top_n Number of top rules to show
#' @param color_by Metric to use for coloring ("support", "confidence", "lift")
#' @param main Plot title
#' @export
plot_rules_matrix <- function(rules, top_n = 50, color_by = "lift", 
                             main = "Rules Matrix") {
  
  # Handle association_rules object
  if (inherits(rules, "association_rules")) {
    rules_df <- rules$rules
  } else if (is.data.frame(rules)) {
    rules_df <- rules
  } else {
    stop("rules must be an association_rules object or data frame")
  }
  
  if (nrow(rules_df) == 0) {
    plot(1, 1, type = "n", main = "No rules to plot")
    text(1, 1, "No association rules found", cex = 1.2)
    return(invisible(NULL))
  }
  
  # Select top rules
  if (top_n < nrow(rules_df)) {
    rules_df <- head(rules_df[order(rules_df[[color_by]], decreasing = TRUE), ], top_n)
  }
  
  # Get unique antecedents and consequents
  unique_antecedents <- unique(rules_df$antecedent)
  unique_consequents <- unique(rules_df$consequent)
  
  # Create matrix
  rule_matrix <- matrix(0, nrow = length(unique_antecedents), 
                       ncol = length(unique_consequents))
  rownames(rule_matrix) <- unique_antecedents
  colnames(rule_matrix) <- unique_consequents
  
  # Fill matrix with color values
  for (i in 1:nrow(rules_df)) {
    ant_idx <- which(unique_antecedents == rules_df$antecedent[i])
    con_idx <- which(unique_consequents == rules_df$consequent[i])
    rule_matrix[ant_idx, con_idx] <- rules_df[[color_by]][i]
  }
  
  # Create color palette
  max_val <- max(rule_matrix)
  min_val <- min(rule_matrix[rule_matrix > 0])
  n_colors <- 100
  colors <- colorRampPalette(c("white", "yellow", "orange", "red"))(n_colors)
  
  # Create plot
  par(mar = c(8, 8, 4, 2))
  
  image(1:ncol(rule_matrix), 1:nrow(rule_matrix), 
        t(rule_matrix[nrow(rule_matrix):1, ]),
        col = colors,
        xlab = "Consequent", ylab = "Antecedent",
        main = main,
        axes = FALSE)
  
  # Add axes
  axis(1, at = 1:ncol(rule_matrix), labels = colnames(rule_matrix), 
       las = 2, cex.axis = 0.8)
  axis(2, at = 1:nrow(rule_matrix), labels = rev(rownames(rule_matrix)), 
       las = 2, cex.axis = 0.8)
  
  # Add color legend
  legend_vals <- seq(min_val, max_val, length.out = 5)
  legend("topright", 
         legend = round(legend_vals, 3),
         fill = colors[round(seq(1, n_colors, length.out = 5))],
         title = paste(toupper(substring(color_by, 1, 1)), substring(color_by, 2), sep = ""))
  
  par(mar = c(5, 4, 4, 2))  # Reset margins
  
  return(invisible(rule_matrix))
}

#' Plot Method for Association Rules
#'
#' Generic plot method for association rules objects.
#'
#' @param x Association rules object
#' @param type Type of plot ("scatter", "network", "frequency", "quality", "matrix")
#' @param ... Additional arguments passed to specific plot functions
#' @export
plot.association_rules <- function(x, type = "scatter", ...) {
  
  if (!type %in% c("scatter", "network", "frequency", "quality", "matrix")) {
    stop("type must be one of: 'scatter', 'network', 'frequency', 'quality', 'matrix'")
  }
  
  if (type == "scatter") {
    plot_rules_scatter(x, ...)
  } else if (type == "network") {
    plot_rules_network(x, ...)
  } else if (type == "frequency") {
    plot_itemset_frequency(x, ...)
  } else if (type == "quality") {
    plot_rules_quality(x, ...)
  } else if (type == "matrix") {
    plot_rules_matrix(x, ...)
  }
} 
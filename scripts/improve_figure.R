# Function for remaking figures to look nice

# Plot remove title
improve_graph <- function(plot, x_label, y_label, n_ticks = 3) {
  # Get current data range
  plot_data <- ggplot_build(plot)$data[[1]]
  x_range <- range(plot_data$x, na.rm = TRUE)
  y_range <- range(plot_data$y, na.rm = TRUE)
  
  # Calculate breaks based on the data range
  x_breaks <- pretty(x_range, n = n_ticks)
  y_breaks <- pretty(y_range, n = n_ticks)
  
  # First remove any existing stat_cor elements
  plot$layers <- plot$layers[
    !sapply(plot$layers, function(x) 
      inherits(x$stat, "StatCor") || 
        (inherits(x$geom, "GeomText") && grepl("R\\s*=", x$data$label[1]))
    )
  ]
  
  # Apply improvements
  improved_plot <- plot +
    labs(
      title = "",  # Remove title
      x = ifelse(is.null(x_label), plot$labels$x, x_label),  # Use custom x label if provided
      y = ifelse(is.null(y_label), plot$labels$y, y_label)   # Use custom y label if provided
    ) + 
    theme(
      text = element_text(family = "Times New Roman", size = 10),
      axis.text = element_text(family = "Times New Roman", size = 10),
      axis.title = element_text(family = "Times New Roman", size = 10),
      plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"),  # Tighter margins
      plot.title = element_text(family = "Times New Roman", size = 8, hjust = 0.5, vjust = -1),  # Centered title
    ) +
    scale_x_continuous(breaks = x_breaks) +
    scale_y_continuous(breaks = y_breaks)
  
  # Extract correlation value from the data
  if (requireNamespace("ggpubr", quietly = TRUE)) {
    cor_test <- with(
      ggplot2::layer_data(plot),
      cor.test(x, y, method = "pearson")
    )
    r_value <- round(cor_test$estimate, 2)
    p_value <- cor_test$p.value
    
    # Format p-value in scientific notation if very small
    if (p_value < 0.001) {
      p_text <- format(p_value, scientific = TRUE, digits = 2)
    } else {
      p_text <- format(p_value, digits = 3)
    }
    
    # Add correlation as a title outside the grid
    improved_plot <- improved_plot + 
      ggtitle(paste0("R = ", r_value, ", p = ", p_text))
  }
  
  # Override the point size for all geom_point layers
  for (i in seq_along(improved_plot$layers)) {
    if ("GeomPoint" %in% class(improved_plot$layers[[i]]$geom)) {
      improved_plot$layers[[i]]$aes_params$size <- 0.1  # Make points very small
    }
    # Make line thinner for any line layers
    if ("GeomLine" %in% class(improved_plot$layers[[i]]$geom) || 
        "GeomSmooth" %in% class(improved_plot$layers[[i]]$geom) ||
        "GeomAbline" %in% class(improved_plot$layers[[i]]$geom)) {
      improved_plot$layers[[i]]$aes_params$size <- 0.25  # Make lines thinner
    }
  }
  
  return(improved_plot)
}
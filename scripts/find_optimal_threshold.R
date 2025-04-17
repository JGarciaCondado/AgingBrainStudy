# Function for finding delta optimum seperation of longitudinal cogntivie trajectories

find_optimal_threshold <- function(data, 
                                   threshold_range,
                                   outcome,
                                   measure,
                                   fixed_formula,
                                   random_formula) {
  
  # Function to fit model for a given threshold
  fit_threshold_model <- function(cut, data) {
    # Create binary amyloid status
    data$delta_status <- ifelse(data[[measure]] >= cut, 1, 0)
    
    # Remove rows with missing outcome
    data <- data[!is.na(data[[outcome]]), ]
    
    # Fit linear mixed effects model
    model <- lme(
      fixed = as.formula(paste(outcome, "~", fixed_formula)),
      random = as.formula(random_formula),
      data = data,
      na.action = 'na.omit',
      method = c('ML'),
      control = lmeControl(opt = "optim", maxIter = 500, msMaxIter = 1000)
    )
    
    return(AIC(model))
  }
  
  
  # Calculate AIC for each threshold
  results <- data.frame(
    threshold = threshold_range,
    aic = sapply(threshold_range, fit_threshold_model, data = data)
  )
  
  # Find optimal threshold
  optimal_threshold <- results$threshold[which.min(results$aic)]
  
  # Create plot
  plot <- ggplot(results, aes(x = threshold, y = aic)) +
    geom_point(color = '#377eb8', size = 2) + 
    geom_line(color = '#377eb8', size = 1.5) +
    geom_vline(xintercept = optimal_threshold, 
               color = "black", 
               linetype = 'longdash') +
    geom_text(aes(x = optimal_threshold + 0.06, 
                  y = max(aic) - (max(aic) - min(aic))/10),
              label = sprintf('Delta = %.2f', optimal_threshold)) +
    labs(x = 'Delta threshold (years)',
         y = "AIC") +
    theme_bw() + 
    theme(
      axis.title = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 12, colour = "#333333"),
      panel.grid = element_blank()
    )
  
  # Return results
  return(list(
    optimal_threshold = optimal_threshold,
    aic_values = results,
    plot = plot
  ))
}
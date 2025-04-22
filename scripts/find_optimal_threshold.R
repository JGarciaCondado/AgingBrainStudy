# Function for finding delta optimum separation of longitudinal cognitive trajectories
# with thresholds expressed in standard deviation units

find_optimal_threshold <- function(data, 
                                   threshold_range,
                                   outcome,
                                   measure,
                                   fixed_formula,
                                   random_formula) {
  
  # Calculate standard deviation of measure for conversion
  measure_sd <- sd(data[[measure]], na.rm = TRUE)
  measure_mean <- mean(data[[measure]], na.rm = TRUE)
  
  # Function to fit model for a given threshold
  fit_threshold_model <- function(cut, data) {
    # Create binary status
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
  
  # Add column for SD units
  results$threshold_sd <- (results$threshold - measure_mean) / measure_sd
  
  # Find optimal threshold
  optimal_threshold <- results$threshold[which.min(results$aic)]
  optimal_threshold_sd <- (optimal_threshold - measure_mean) / measure_sd
  
  # Create plot in standard deviation units
  plot <- ggplot(results, aes(x = threshold_sd, y = aic)) +
    geom_point(color = '#377eb8', size = 2) + 
    geom_line(color = '#377eb8', size = 1.5) +
    geom_vline(xintercept = optimal_threshold_sd, 
               color = "black", 
               linetype = 'longdash') +
    geom_text(aes(x = optimal_threshold_sd + 0.2, 
                  y = max(aic) - (max(aic) - min(aic))/10),
              label = sprintf('Delta = %.2f SD', optimal_threshold_sd)) +
    labs(x = 'Delta threshold (standard deviations)',
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
    optimal_threshold_sd = optimal_threshold_sd,
    aic_values = results,
    plot = plot,
    measure_mean = measure_mean,
    measure_sd = measure_sd
  ))
}
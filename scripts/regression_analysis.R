# File with all required scripts
library(ppcor)

regression_analysis <- function(data, x, y, covars) {
  
  # Filter nan values
  data <- data[complete.cases(data[, c(x, y, covars)]), ]
  
  # Scatter plot with correlation
  scatter_plot <- ggplot(data, aes(x = !!sym(y), y = !!sym(x))) +
    geom_point() +
    geom_smooth(method = 'lm', formula = y ~ x) +
    stat_cor(method = 'pearson') +
    labs(title = paste(x, "vs", y),
         x = y,
         y = x)
  
  # Create a formula for the regression model
  covar_formula <- paste(covars, collapse = " + ")
  full_formula <- as.formula(paste(x, "~", y, "+", covar_formula))
  
  # Linear regression model
  ols_results <- lm(as.formula(paste(x, "~", y)), data = data)
  ols_full_results <- lm(full_formula, data = data)
  
  # Model summary
  model_summary <- summary(ols_results)
  model_full_summary <- summary(ols_full_results)
  
  # Plot model
  pred_plot <- plot_model(ols_results, type = c('pred'), terms = c(y), show.values = TRUE) +
    geom_point(data = data, aes(x = !!sym(y), y = !!sym(x)), 
               inherit.aes = F) + 
    labs(title = paste(x, "vs", y),
         x = y,
         y = x)
  
  # Partial correlation
  covariates <- data[, covars, drop = FALSE]
  covariate_matrix <- model.matrix(~ . - 1, data = covariates)
  partial_corr <- pcor.test(data[x], data[y], covariate_matrix, method = "pearson")
  
  # Return a list with all results
  return(list(
    scatter_plot = scatter_plot,
    simple_model = ols_results,
    regression_model = ols_full_results,
    model_summary = model_summary,
    model_full_summary = model_full_summary,
    predicted_values_plot = pred_plot,
    partial_correlation = partial_corr
  ))
}
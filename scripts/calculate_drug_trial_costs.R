# Function to calculate drug trial costs
calculate_drug_trial_costs <- function(
    years,                     # Total years of the trial
    participants_required,     # Number of participants required in the trial
    screen_pass_rate,       # Percentage of participants that screen fail (0-100)
    screen_cost,               # Cost per screening (default value)
    mri_cost,                  # Cost per MRI screening (default value)
    treatment_cost_per_year,   # Treatment cost per participant per year (default value)
    has_mri_screening = FALSE, # Whether there are MRI screenings (TRUE/FALSE)
    mri_pass_rate = 0      # Percentage of participants that don't pass MRI screening (0-100)
) {
  
  # Calculate number of participants to screen
  if (has_mri_screening) {
    # If using MRI screening, account for both regular screening and MRI screening failures
    participants_to_screen = participants_required / (screen_pass_rate * mri_pass_rate)
  } else {
    # Without MRI, only account for regular screening failures
    participants_to_screen = participants_required / screen_pass_rate
  }
  
  # Calculate screening costs
  if (has_mri_screening) {
    # For participants who pass initial screening, they also get an MRI
    participants_passing_initial = participants_to_screen * screen_pass_rate
    screening_cost = (participants_to_screen * screen_cost) + (participants_passing_initial * mri_cost)
  } else {
    screening_cost = participants_to_screen * screen_cost
  }
  
  # Calculate treatment costs
  # Assuming half of the enrolled participants receive the treatment (other half is control)
  treatment_cost = treatment_cost_per_year * years * participants_required / 2
  
  # Calculate total costs
  total_cost = screening_cost + treatment_cost
  
  # Create and return results list
  results <- list(
    screening_cost = screening_cost,
    participants_to_screen = participants_to_screen,
    treatment_cost = treatment_cost,
    total_cost = total_cost
  )
  
  return(results)
}
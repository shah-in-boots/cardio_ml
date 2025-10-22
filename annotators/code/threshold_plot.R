# raw probabilities -------------------------------------------------------
# Experimental
make_predictions_integer <- function(predictions, threshold = 0.5) {
  #' @param predictions A single ECG of probabilities

  # Dimensions
  n_timesteps <- dim(predictions)[1]
  n_classes   <- dim(predictions)[2]  # should be 4
  n_leads     <- dim(predictions)[3]
  
  # Initialize with 1s
  predictions_integer <- matrix(0, nrow = n_timesteps, ncol = n_leads)
  
  # Extract the slice for this sample: timesteps × classes × leads
  
  # Loop over leads
  for (lead in 1:n_leads) {
    # Extract probabilities for this lead: timesteps × classes
    lead_probs <- predictions[,,lead]
    
    # Focus only on classes 2:4 (P, QRS, T)
    wave_probs <- lead_probs[, 2:4, drop = FALSE]
    
    # For each timestep, check which classes exceed threshold
    above_thresh <- wave_probs > threshold
    
    # Count how many classes exceed threshold per timestep
    n_above <- rowSums(above_thresh)
    
    # Case 1: exactly one class above threshold → assign that class
    one_class_idx <- which(n_above == 1)
    if (length(one_class_idx) > 0) {
      # Find which class was above threshold
      class_ids <- apply(above_thresh[one_class_idx, , drop = FALSE], 1, which.max)
      predictions_integer[one_class_idx, lead] <- class_ids
    }
    
    # Case 2: multiple classes above threshold → assign 4
    multi_class_idx <- which(n_above > 1)
    if (length(multi_class_idx) > 0) {
      predictions_integer[multi_class_idx, lead] <- 4
    }
  }
  
  return(predictions_integer)
}

make_predictions_integer_abs <- function(predictions) {
  # Dimensions
  n_timesteps <- dim(predictions)[1]
  n_classes   <- dim(predictions)[2]  # should be 4
  n_leads     <- dim(predictions)[3]
  
  # Initialize with 0s (or 1s if you prefer "no wave" as default)
  predictions_integer <- matrix(0, nrow = n_timesteps, ncol = n_leads)
  
  # Loop over leads
  for (lead in 1:n_leads) {
    # Extract probabilities for this lead: timesteps × classes
    lead_probs <- predictions[,,lead]
    
    # Find the class with the maximum probability at each timestep
    # which.max() returns the index (1–4) of the max class
    max_classes <- apply(lead_probs, 1, which.max)
    
    # Assign directly
    predictions_integer[, lead] <- max_classes - 1
  }
  
  return(predictions_integer)
}



# predict ---------------------------------------------------------------------

raw_preds <- predict_ecgs_raw(input = ecgs)

# plot --------------------------------------------------------------------
idx <- 2
lead <- 6
predictions_integer <- make_predictions_integer(raw_preds[idx,,,],
                                                threshold = 0.5)

old_predictions_integer <- make_predictions_integer_abs(raw_preds[idx,,,])


  
lead_name_list <- c('I','II','III','AVR','AVL','AVF','V1','V2','V3','V4','V5','V6')
lead_name <- lead_name_list[lead]

sig <- ecg_filter(ecgs[[idx]]$signal[[lead+1]])
plot_func2(sig,predictions_integer[,lead])
# plot_func2(sig,old_predictions_integer[,lead])

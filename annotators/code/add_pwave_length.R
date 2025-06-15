library(keras)
load('../models/model_log.RData')
load('../uih_ecgs.RData')
source('annotator_prep_functions.R')

normalize <- TRUE


for (row in 350:499) {
  model_name <- model_log$name[row]
  # Extract the part after "cnn_" and before "2025"
  lead <- sub(".*bilstm_([^_]+)_2025.*", "\\1", model_name)
  model <- load_model_tf(paste0('../models/',model_name,'.h5'))

  # Create matrix
  uih_samples <- do.call(rbind, lapply(1:length(uih_ecgs), function(idx)
    uih_ecgs[[idx]]$signal[[lead]]))
  filtered <- array(0, dim(uih_samples))
  for (i in 1:nrow(uih_samples)) {
    filtered[i, ] <- ecg_filter(uih_samples[i, ])
  }
  for (i in 1:nrow(uih_samples)) {
    filtered[i, ] <- (filtered[i, ] - min(filtered[i, ])) / (max(filtered[i, ]) - min(filtered[i, ])) * 100
  }
  
  # Predict
  uih_predictions <- model %>% predict(filtered)
  uih_predictions_integer <- array(0, c(nrow(uih_predictions), ncol(uih_predictions)))
  for (i in 1:nrow(uih_predictions)) {
    uih_predictions_integer[i, ] <- max.col(uih_predictions[i, , ])
  }
  #convert from dimension value 1,2,3,4 to 0,1,2,3
  uih_predictions_integer <- uih_predictions_integer - 1
  
  for (i in 1:nrow(uih_predictions_integer)) {
    sample <- uih_predictions_integer[i,]
    rle_obj <- rle(sample)  # Run-length encoding
    p_wave_lengths <- rle_obj$lengths[rle_obj$values == 1]
  }
  p_wave_length <- median(p_wave_lengths,na.rm = TRUE)
  
  for (i in 1:nrow(uih_predictions_integer)) {
    sample <- uih_predictions_integer[i,]
    rle_obj <- rle(sample)  # Run-length encoding
    t_wave_lengths <- rle_obj$lengths[rle_obj$values == 3]
  }
  t_wave_length <- median(t_wave_lengths,na.rm = TRUE)
  
  model_log$p_wave_length[row] <- p_wave_length
  model_log$t_wave_length[row] <- t_wave_length
  
  print(paste(row))
}

save(model_log,file='../models/model_log.RData')
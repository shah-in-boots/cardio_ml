# Prep --------------------------------------------------------------------
source('annotator_prep_functions.R')
annotator_style <- 2
lead <- 3
rounds <- 4
max_noise = 0.05 # 0.05, 0.03
dilate_range <- c(0.03, 0.05)
filter <- FALSE
epochs_to_save <- c(15,20,25) # list which epochs you want to save
bilstm_layers <- 200 # original: 200
normalize <- TRUE
number_of_derivs = 2

out <- prep_ludb(
  lead = lead,
  annotator_style = 2,
  rounds = rounds,
  max_noise = max_noise,
  dilate_range = dilate_range,
  filter = filter,
  normalize = normalize,
  number_of_derivs = number_of_derivs,
  mask_value = -1
)

training_signal <- out$training_signal
training_annotations <- out$training_annotations
testing_signal <- out$testing_signal
testing_annotations <- out$testing_annotations



# Build Model -------------------------------------------------------------------
library(keras)
library(tensorflow)

# Model parameters
input_length <- 5000                # number of time steps in the ECG
num_channels <- number_of_derivs+1  # single lead input
units_lstm <- bilstm_layers         # LSTM units
mask_value <- out$mask_value        # if you wish to mask zero values in the input
activation <- "softmax"             # final activation for multi-class predictions
model_type <- 'bilstm_cnn'
kernel <- c(5,3)                    # size of kernel for each convolutional layer
filters <- c(32,64)                 # number of filters for each convolutional layer


# Build and name model
epochs <- max(epochs_to_save)
num_classes <- length(unique(as.vector(training_annotations))) # classes: no wave, P, QRS, T
leads <- c('I','II','III','AVR','AVL','AVF','V1','V2','V3','V4','V5','V6')
date_time <- format(Sys.time(), "%Y%m%d_%H%M%S")
model_name_template <- paste0(model_type, '_', leads[lead], '_', date_time, '_epoch_')
model_name_path <- paste0("../models/", model_name_template)

# Save model at epochs of choice:
epoch_save_callback <- callback_lambda(
  on_epoch_end = function(epoch, logs) {
    if (epoch %in% (epochs_to_save - 1)) {
      # Save at epochs 20, 30, 40 (0-based index)
      model_name <- paste0(model_name_path, epoch + 1, ".h5")
      save_model_tf(model, model_name)
      print(paste("Saved model at epoch", epoch + 1))
    }
  }
)

# Define the model input
inputs <- layer_input(shape = c(input_length, num_channels), dtype = 'float32') |>
  layer_masking(mask_value = mask_value) # Masking input values equal to 0

# -- Convolutional Feature Extraction Block --
# This block applies 1D convolutions to extract local features from the ECG slice.
cnn_block <- inputs %>%
  # First convolution: using a larger kernel can help capture the broader wave morphology.
  layer_conv_1d(
    filters = filters[1],
    kernel_size = kernel[1],
    activation = 'relu',
    padding = "same"
  ) %>%
  # Max pooling to downsample and focus on the most prominent features.
  layer_max_pooling_1d(pool_size = 2) %>%
  # Second convolution: deeper features at a possibly finer resolution.
  layer_conv_1d(
    filters = filters[2],
    kernel_size = kernel[2],
    activation = 'relu',
    padding = "same"
  ) %>%
  layer_max_pooling_1d(pool_size = 2)

# -- Recurrent Block for Temporal Modeling --
# The output of the convolutional block is still sequential.
# The bidirectional LSTM helps capture signals from both past and future contexts.
lstm_block <- cnn_block %>%
  bidirectional(layer_lstm(
    units = units_lstm,
    return_sequences = TRUE,  # important to output a sequence for per-time-step predictions
    activation = 'tanh'
  ))

# -- Upsampling Layer --
upsampled_block <- lstm_block %>%
  layer_upsampling_1d(size = 4)  # Upsample by a factor of 2 to restore time steps

# -- Output Layer --
# A dense layer produces classification probabilities at every time step.
outputs <- upsampled_block %>%
  layer_dense(units = num_classes, activation = activation, name = "predictions")

# Define, compile, and summarize the model
model <- keras_model(inputs = inputs, outputs = outputs, name = "CNN_LSTM_ECG_Model")

model %>% compile(
  optimizer = optimizer_adam(),
  loss = "sparse_categorical_crossentropy",
  metrics = "accuracy"
)

# Train model -------------------------------------------------------------
# Print to command line:
cat("lead:",lead,'/',leads[lead],"\n",
    "model_type:",model_type,"\n",
    "model_name:",model_name_path,"\n",
    "bilstm_layers:",bilstm_layers,"\n",
    "epochs:",epochs,"\n",
    "activation:",activation,"\n",
    "rounds:",rounds,"\n",
    "max_noise:",max_noise,"\n",
    "dilate_range:",dilate_range,"\n",
    "annotator_style:",annotator_style,"\n",
    "filter:",filter,"\n",
    "normalize:",out$normalize,"\n"
)

# Train
start <- Sys.time()
history <- model |> fit(
  training_signal,
  training_annotations,
  epochs = epochs,
  validation_data = list(testing_signal, testing_annotations),
  verbose = 2,
  callbacks = list(epoch_save_callback)
)

# output_name <- paste0("../models/",model_name)
# save_model_tf(model, output_name)

end <- Sys.time()
time_spent <- end - start
cat("time_spent:", time_spent, "\n")
# bilstm: bilstm_layers, activation
# cnn_bilstm_attn: dropout, bilstm_layers
# UNET: dropout, filters
# **specify training/testing samples from out$

# Add to

# Test model on testing LUDB samples--------------------------------------------------------------
new_row <- c()
load('../uih_ecgs.RData')

for (checkpoint_number in 1:length(epochs_to_save)) {
  # Load model
  model_name <- paste0(model_name_template, epochs_to_save[checkpoint_number])
  model_path <- paste0('../models/', model_name)
  
  model <- load_model_tf(paste0(model_path, '.h5'))
  
  predictions <- model %>% predict(testing_signal)
  predictions_integer <- array(0, c(nrow(predictions), ncol(predictions)))
  for (i in 1:nrow(predictions)) {
    predictions_integer[i, ] <- max.col(predictions[i, , ])
  }
  #convert from dimension value 1,2,3,4 to 0,1,2,3
  predictions_integer <- predictions_integer - 1
  
  
  # LUDB test sample Analysis
  confusion <- confusion_analysis(predictions = predictions_integer, actual = testing_annotations)
  wave_counter <- check_ann_prog_ludb_helper(predictions_integer = predictions_integer,
                                             testing_annotations = testing_annotations)
  
  
  # UIH test samples Analysis (nonlabelled) ---------------------------------
  # note: UIH lead order is AVF AVL AVR. Therefore, use names to match (ie leads[lead] ), **NOT** number (ie lead)
  uih_samples <- do.call(rbind, lapply(1:length(uih_ecgs), function(idx)
    uih_ecgs[[idx]]$signal[[leads[lead]]]))
  
  filtered <- array(0, dim(uih_samples))
  for (i in 1:nrow(uih_samples)) {
    filtered[i, ] <- ecg_filter(uih_samples[i, ])
  }
  
  if (normalize) {
    for (i in 1:nrow(uih_samples)) {
      filtered[i, ] <- (filtered[i, ] - min(filtered[i, ])) / (max(filtered[i, ]) - min(filtered[i, ])) * 100
    }
  }
  
  uih_input <- array(NA,c(dim(filtered),number_of_derivs+1))
  if (number_of_derivs > 0) {
    for (i in 1:nrow(filtered)) {
      
      signal <- add_derivs(signal = filtered[i,],
                           number_of_derivs = number_of_derivs,
                           mask_value = mask_value)
      uih_input[i,,] <- signal
    }
  } else {
    uih_input[,,1] <- filtered
  }

  # Predict
  uih_predictions <- model %>% predict(uih_input)
  uih_predictions_integer <- array(0, c(nrow(uih_predictions), ncol(uih_predictions)))
  for (i in 1:nrow(uih_predictions)) {
    uih_predictions_integer[i, ] <- max.col(uih_predictions[i, , ])
  }
  #convert from dimension value 1,2,3,4 to 0,1,2,3
  uih_predictions_integer <- uih_predictions_integer - 1
  
  
  
  # Check P>QRS>T progression.
  #   prog_count_revised: fills in 2 p-waves with a gap of <20 indicies. Same for QRS, T
  prog_count <- array(NA, nrow(uih_predictions_integer))
  prog_count_revised <- array(NA, nrow(uih_predictions_integer))
  
  uih_predictions_integer_revised <- array(0, dim(uih_predictions_integer))
  for (i in 1:nrow(uih_predictions_integer)) {
    uih_predictions_integer_revised[i, ] <- fill_wave_gaps(uih_predictions_integer[i, ], 20)
  }
  
  for (i in 1:nrow(uih_predictions_integer)) {
    prog_count[i] <- check_ann_prog(annotation = uih_predictions_integer[i, ])
    prog_count_revised[i] <- check_ann_prog(uih_predictions_integer_revised[i, ])
  }
  
  uih_prog <- sum(prog_count == 0) / length(prog_count)
  uih_prog_revised <- sum(prog_count_revised == 0) / length(prog_count_revised)
  
  
  # Confidence: average max probability value across all time points
  confidence <- round(mean(apply(uih_predictions, c(1, 2), max)), 4)
  
  
  # Check wave progression using EGM::find_RPeaks()
  library(dplyr)
  uih_Rpeaks <- do.call(rbind, lapply(1:nrow(uih_input[,,1]), function(idx) {
    round(check_ann_prog_RPeaks(uih_input[idx, ,1], uih_predictions_integer[idx, ]),
          2)
  }))
  
  Rpeaks_prog <- data.frame(
    QRS = sum(uih_Rpeaks$missed_QRS + uih_Rpeaks$duplicate_QRS) / nrow(uih_Rpeaks),
    P = sum(uih_Rpeaks$missed_P + uih_Rpeaks$duplicate_P) / nrow(uih_Rpeaks),
    T = sum(uih_Rpeaks$missed_T + uih_Rpeaks$duplicate_T) / nrow(uih_Rpeaks)
  )
  
  Rpeaks_prog$total <- sum(Rpeaks_prog)
  print(Rpeaks_prog)
  str(Rpeaks_prog)
  
  Rpeaks_prog_P <- Rpeaks_prog$P
  Rpeaks_prog_QRS <- Rpeaks_prog$QRS
  Rpeaks_prog_T <- Rpeaks_prog$T
  Rpeaks_prog_total <- Rpeaks_prog$total
  
  # Find median P, T wave length:
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
  
  print(data.frame(
    name = model_name,
    type = model_type,
    ann_style = annotator_style,
    lead = lead,
    bilstm_layers = bilstm_layers,
    dropout = NA,
    filters = NA,
    epochs = epochs_to_save[checkpoint_number],
    time = round(time_spent * (epochs_to_save[checkpoint_number] / epochs), 2),
    training_samples = I(list(out$training_samples)),
    confusion = I(list(confusion)),
    dilate_range = I(list(dilate_range)),
    max_noise = max_noise,
    rounds = rounds,
    filter = filter,
    wave_counter = I(list(wave_counter)),
    normalize = normalize,
    uih_prog = uih_prog,
    uih_prog_revised = uih_prog_revised,
    confidence = confidence,
    Rpeaks_prog_P = Rpeaks_prog_P,  # Ensure proper reference
    Rpeaks_prog_QRS = Rpeaks_prog_QRS,
    Rpeaks_prog_T = Rpeaks_prog_T,
    Rpeaks_prog_total = Rpeaks_prog_total,
    p_wave_length = p_wave_length,
    t_wave_length = t_wave_length
  ))
  
  # Add new row for model_log
  new_row <- rbind(
    new_row,
    data.frame(
      name = model_name,
      type = model_type,
      ann_style = annotator_style,
      lead = lead,
      bilstm_layers = bilstm_layers,
      dropout = NA,
      filters = I(list(filters)),
      epochs = epochs_to_save[checkpoint_number],
      time = round(time_spent * (epochs_to_save[checkpoint_number] / epochs), 2),
      training_samples = I(list(out$training_samples)),
      confusion = I(list(confusion)),
      dilate_range = I(list(dilate_range)),
      max_noise = max_noise,
      rounds = rounds,
      filter = filter,
      wave_counter = I(list(wave_counter)),
      normalize = normalize,
      uih_prog = uih_prog,
      uih_prog_revised = uih_prog_revised,
      confidence = confidence,
      Rpeaks_prog_P = Rpeaks_prog_P,  # Ensure proper reference
      Rpeaks_prog_QRS = Rpeaks_prog_QRS,
      Rpeaks_prog_T = Rpeaks_prog_T,
      Rpeaks_prog_total = Rpeaks_prog_total,
      p_wave_length = p_wave_length,
      t_wave_length = t_wave_length,
      kernel = I(list(kernel)),
      derivs = number_of_derivs
    )
  )
  
}

print(new_row)
cat("New row names:","\n")
print(names(new_row))
print(names(new_row), quote = TRUE) 

cat("Class of normalize and uih_prog:","\n")
print(class(normalize))
print(class(uih_prog))

cat("Finished training and analyzing. Adding data to model log...", "\n")

# Write to log ------------------------------------------------------------
# Add lock to avoid overwriting
library(filelock)
logFile <- '../models/model_log.RData'
lockFile <- paste0(logFile, ".lock")  # Create a lock file
# Acquire the lock before loading or saving
lock <- lock(lockFile)
load(logFile)
print(colnames(model_log))

model_log <- rbind(model_log, new_row)
save(model_log, file = logFile)

# Release the lock
unlock(lock)

cat("Data added successfully. Job complete.", "\n") 
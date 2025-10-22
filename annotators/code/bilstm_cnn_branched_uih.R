# Prep --------------------------------------------------------------------
source('annotator_prep_functions.R')
library(abind)

# Input data parameters:
annotator_style <- 2
lead <- 7
rounds <- c(4,8)                    # number of times the dataset is duplicated. 1st index: LUDB set. 2nd index: UIH set.
max_noise = 0.01 # 0.05, 0.03       # maximum noise that can be introduced to the signal. RECOMMEND: 0.01
dilate_range <- c(0.03, 0.05)       # factor to compress/dilate
filter <- FALSE                     # whether to filter the data on input
epochs_to_save <- c(15,20,25,30,35) # list which epochs you want to save
bilstm_layers <- 200                # Layers in bilstm block. original: 200
normalize <- TRUE                   # whether to scale signal 0 to 100. Recommend: TRUE
number_of_derivs = 2                # number of derivatives to include on input. 0 would be just the signal vector.

# Model parameters
input_length <- 5000                # number of time steps in the ECG
activation <- "softmax"             # final activation for multi-class predictions
kernel <- c(7, 31, 81)              # size of kernel for each convolutional layer. Base: c(5,3)
filters <- c(32, 32, 32)            # number of filters for each convolutional layer. For this model type, its the filters for the convolutional layers
proj_filters <- 128                 # number of filters for the projection layers, which is used after the parallel convolutional layers are combined 
model_type <- 'bilstm_cnn_branched_expanded'         # expanded: training data includes LUDB and supplemental UIH ECGs
mask_value <- -1

# Naming
leads <- c('I','II','III','AVR','AVL','AVF','V1','V2','V3','V4','V5','V6')
date_time <- format(Sys.time(), "%Y%m%d_%H%M%S")
model_name_template <- paste0(model_type, '_', leads[lead], '_', date_time, '_epoch_')
model_name_path <- paste0("../models/", model_name_template)

# Define the UIH ECGs which have been annotated
load('../uih_training_leads.RData')
uih_ecg_range <- uih_training_leads[[leads[lead]]] # selects ECGs which were pre-annotated


# Build LUDB input data
out <- prep_ludb(
  lead = lead,
  annotator_style = 2,
  rounds = rounds[1],
  max_noise = max_noise,
  dilate_range = dilate_range,
  filter = filter,
  normalize = normalize,
  number_of_derivs = number_of_derivs,
  mask_value = mask_value,
  data='ludb'
)
num_classes <- length(unique(as.vector(out$training_annotations))) # classes: no wave, P, QRS, T

# Add UIH input data
load('../deid_uih_ecgs.RData')
input_uih_ecgs <- deid_uih_ecgs[uih_ecg_range]
out_uih <- prep_ludb(
  lead = lead,
  annotator_style = 2,
  rounds = rounds[2],
  max_noise = max_noise,
  dilate_range = dilate_range,
  filter = TRUE,
  normalize = normalize,
  number_of_derivs = number_of_derivs,
  mask_value = mask_value,
  data = input_uih_ecgs
)

# Combine datasets:
training_signal <- abind(out$training_signal,out_uih$training_signal,along = 1)
training_annotations <- abind(out$training_annotations,out_uih$training_annotations,along = 1)
testing_signal <- abind(out$testing_signal,out_uih$testing_signal,along = 1)
testing_annotations <- abind(out$testing_annotations,out_uih$testing_annotations,along = 1)

training_samples <- data.frame(samples = c(out$training_samples,out_uih$training_samples),
                               dataset = c(rep("ludb",length(out$training_samples)),rep("uih",length(out_uih$training_samples))))
testing_samples <- data.frame(samples = c(out$testing_samples,out_uih$testing_samples),
                              dataset = c(rep("ludb",length(out$testing_samples)),rep("uih",length(out_uih$testing_samples))))
# Build model
model <- build_bilstm_cnn_branched(
  bilstm_layers = bilstm_layers,
  number_of_derivs = number_of_derivs,
  mask_value = mask_value,
  kernel = kernel,
  branch_filters = filters,
  proj_filters = proj_filters,
  activation = activation,
  num_classes = num_classes
)


# Define callbacks
epochs <- max(epochs_to_save)
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
  epochs =   epochs,
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
# load('../uih_ecgs.RData')

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
  ludb_training_indices <- which(testing_samples$dataset == 'ludb')
  uih_training_indices <- which(testing_samples$dataset == 'uih')
  
  confusion <- confusion_analysis(predictions = predictions_integer[ludb_training_indices,], 
                                  actual = testing_annotations[ludb_training_indices,])
  wave_counter <- check_ann_prog_ludb_helper(predictions_integer = predictions_integer[ludb_training_indices,],
                                             testing_annotations = testing_annotations[ludb_training_indices,])
  
  confusion_uih <- confusion_analysis(predictions = predictions_integer[uih_training_indices,], 
                                      actual = testing_annotations[uih_training_indices,])
  wave_counter_uih <- check_ann_prog_ludb_helper(predictions_integer = predictions_integer[uih_training_indices,],
                                                 testing_annotations = testing_annotations[uih_training_indices,])
  
  # UIH test samples Analysis (nonlabelled) ---------------------------------
  # note: UIH lead order is AVF AVL AVR. Therefore, use names to match (ie leads[lead] ), **NOT** number (ie lead)
  uih_samples <- do.call(rbind, lapply(1:length(deid_uih_ecgs), function(idx)
    deid_uih_ecgs[[idx]]$signal[[leads[lead]]]))
  
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
    training_samples = I(list(training_samples)),
    confusion = I(list(confusion)),
    dilate_range = I(list(dilate_range)),
    max_noise = max_noise,
    rounds = I(list(rounds)),
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
      training_samples = I(list(training_samples)),
      confusion = I(list(confusion)),
      confusion_uih = I(list(confusion_uih)),
      dilate_range = I(list(dilate_range)),
      max_noise = max_noise,
      rounds = I(list(rounds)),
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
# Build training and testing  sets---------------------------------------
#' @descrption Input: custom list of LUDB annotations and signal. Splits into training/testing. 
#' Takes random continuous segment between 2 and 8 seconds, zeros the rest. Randomly shifts the segment across the vector
#' Then normalizes the signal for each sample.
#' 
#' @param dilate_range For training set: 2 element vector. Defines the range the samples will be randomly. Recommend: c(0.03,0.05) (ie 3 to 5%)
#' 
#' @param max_noise For training set: maximum standard deviation of the added Gaussian noise. ECGs are scaled from 0 to 1. Recommend: between 0.01 to 0.07
#' 
#' @param rounds number of times the training dataset will be duplicated. Note: a value of 2 will double the length of the training set. A value of 1 will not change the length. Recommend: 3 or 4
#' 
#' @param filter if the signal should be run thru a bandpass filter to remove baseline wander, high frequency noise. Recommend: FALSE
#' 
#' @param number_of_derivs if/how many derivatives should be included. If set to 1: input would be the original signal AND the first derviative (matrix form)
prep_ludb <- function(lead,
                      annotator_style = 2,
                      split = 0.7,
                      dilate_range = c(0.03,0.05),
                      max_noise = 0.05,
                      rounds = 1,
                      filter,
                      number_of_derivs,
                      mask_value = -1,
                      normalize = TRUE) {
  library(stats)
  library(signal)
  
  # Load LUDB set
  load('../ludb_set.RData')
  
  # Assign the annotation style for ML input:
  use_func <- paste0('ann_wfdb2continuous', annotator_style)
  
  #         ann_wfdb2continuous1: 1 0 0 0 1 0 0 0 2 0 0 2 ...
  #         ann_wfdb2continuous2: 1 1 1 1 1 0 0 0 2 2 2 2 ...
  #         ann_wfdb2continuous3: 1 2 2 2 3 0 0 0 4 5 5 6 ...
  
  # Split into train/test
  split <- 0.7 # % of samples going to training set
  
  leads <- c("i","ii","iii","avr","avl","avf","v1","v2","v3","v4","v5","v6")
  
  length <- 5000
  sampling_frequency <- 500
  annotation_function <- get(use_func) # Get the function object
  
  # Split sets:
  samples <- length(ludb_set)
  sample_size <- ceiling(split * samples)
  training_samples <- sample(samples, sample_size)
  testing_samples <- 1:samples
  testing_samples <- testing_samples[-training_samples]
  
  training_set <- ludb_set[c(training_samples)]
  testing_set <- ludb_set[c(testing_samples)]
  
  training_signal      <- array(NA, c(length(training_set)*rounds, length, number_of_derivs+1))
  training_annotations <- array(NA, c(length(training_set)*rounds, length))
  testing_signal       <- array(NA, c(length(testing_set), length, number_of_derivs+1))
  testing_annotations  <- array(NA, c(length(testing_set), length))
  
  # Prep for calculating indices for time range between 2 and 8 seconds in training set
  start_index <- 2 * sampling_frequency + 1
  end_index <- 8 * sampling_frequency
  interval_length <- 4 * sampling_frequency
  
  
  # Training Set
  for (round in 1:rounds) {
    for (sample in 1:length(training_set)) {
      # **Take a random 4 second interval of signal**:
      
      # Randomly select the starting point for the 4 second interval
      random_start <- sample(start_index:(end_index - interval_length), 1)
      # Define the 4-second interval
      random_interval <- random_start:(random_start + interval_length - 1)
      
      signal <- training_set[[sample]]$signal[[leads[lead]]]
      ann <- annotation_function(training_set[[sample]]$annotation[[leads[lead]]])
      
      if (filter) {
        signal <- ecg_filter(signal = signal)
      }
      
      # Set indices outside of interval to 0
      # signal[-random_interval] <- 0
      # ann[-random_interval] <- 0
      
      # Random interval signal/ann vectors
      signal_iso <- signal[random_interval]
      ann_iso <- ann[random_interval]
      
      
      # Compress signal:
      # Generate a random compression/dilation factor
      
      if (round != 1) {
        compression_factor <- 1 + sample(c(-1, 1), size = 1) * runif(1, min = dilate_range[1], max = dilate_range[2]) # Random factor between 1 +/- (3 to 5)%
        new_length <- round(length(signal_iso) * compression_factor) # Define new length for the signal based on compression/dilation factor
        # Resample the ECG signal to the new length
        compressed_signal <- resample(signal_iso, p = new_length, q = length(signal_iso))
        # Cannot use same interpolation to do ann vector, which contains *integers* (not decimals). ie, must ensure termini of QRS, T don't become 1s or 2s if rounding
        new_indices <- seq(from = 1,
                           to = length(ann_iso),
                           length.out = new_length) # Create new indices from 1 to length(ann_iso) with the new length
        compressed_ann <- ann_iso[round(new_indices)] # For each new index, use the nearest original annotation
      } else {
        # if on the 1st round, do not dilate/compress, or add noise
        compressed_signal <- signal_iso
        compressed_ann <- ann_iso
        new_length <- length(compressed_signal)
      }
      
      
      
      
      # Generate Gaussian random noise
      # Compress signal:
      # Generate a random compression/dilation factor
      if (round != 1) {
        # if on the 1st round, do not dilate/compress, or add noise
        sigma <- runif(1, min = 0, max = max_noise) # vary the noise SD from 0 to max_noise parameter
        range <- range(compressed_signal)
        noise <- rnorm(length(compressed_signal),
                       mean = 0,
                       sd = sigma * abs(range[1] - range[2])) # create random noise
        compressed_signal <- compressed_signal + noise # Add the noise to the original ECG signal
      } 
      
      # Normalize signal
      if (normalize) {
        compressed_signal <- (compressed_signal - min(compressed_signal)) / (max(compressed_signal) - min(compressed_signal)) * 100
      }
      
      # Now, randomly place the random signal/ann interval:
      # Determine a random starting point for the shifted interval
      vector_length <- length(signal)
      shifted_start <- sample(1:(vector_length - new_length), 1)
      
      # Generate the new shifted indices
      shifted_indices <- shifted_start:(shifted_start + new_length - 1)
      
      # Create new vectors initialized with zeros
      shifted_signal <- array(mask_value, c(vector_length,number_of_derivs+1))
      shifted_ann <- rep(0, vector_length)
      
      # Add the derivatives as needed
      if (number_of_derivs > 0) {
        compressed_signal <- add_derivs(signal = compressed_signal,
                                        number_of_derivs = number_of_derivs,
                                        mask_value = mask_value)
      }
      
      # Shift the non-zero values to the new interval
      shifted_signal[shifted_indices,] <- compressed_signal
      shifted_ann[shifted_indices] <- compressed_ann
 
      
      # Add sample to matrix
      index <- sample + (round-1)*length(training_set)
      training_signal[index, , ] <- shifted_signal # signal
      training_annotations[index, ] <- shifted_ann # ann
    }
  }
  
  # Testing Set (no need to select a random interval)
  for (sample in 1:length(testing_set)) {
    signal <- testing_set[[sample]]$signal[[leads[lead]]]
    ann <- annotation_function(testing_set[[sample]]$annotation[[leads[lead]]])
    
    # Filter signal:
    signal <- ecg_filter(signal = signal)
    
    # Normalize signal
    if (normalize) {
      signal <- (signal - min(signal)) / (max(signal) - min(signal)) * 100
    }
    
    if (number_of_derivs > 0) {
      signal <- add_derivs(signal = signal,
                                      number_of_derivs = number_of_derivs,
                                      mask_value = mask_value)
    }
    
    # Set values of first and last 1000 indices to 0 (could also do 1st and last annotated indices)
    signal[c(1:1000, 4001:5000),] <- mask_value
    ann[c(1:1000, 4001:5000)] <- 0
    
    # Add sample to matrix
    testing_signal[sample, , ] <- signal
    testing_annotations[sample, ] <- ann
  }
  
  # Prepare output
  output <- list(
    training_signal = training_signal,
    training_annotations = training_annotations,
    testing_signal = testing_signal,
    testing_annotations = testing_annotations,
    training_samples = training_samples,
    testing_samples = testing_samples,
    mask_value = mask_value,
    normalize = normalize
  )
  
  return(output)
  
}

# ann_wfdb2continuous1 ----------------------------------------------------
#' @description format [1 0 0 0 1 0 0 0 0 2 0 0 0 2]
ann_wfdb2continuous1 <- function(object) {
  if (any(class(object) == 'egm')) {
    length <- nrow(object$signal)
    wfdb_ann <- object$annotation
  } else if (any(class(object) == 'annotation_table')) {
    wfdb_ann <- object
  }
  
  output <- rep(0,length)
  
  
  # P wave --------------------------
  # Find P indices, then verify onset and offset
  p_ind <- which(wfdb_ann$type == 'p')
  p_ind <- p_ind[wfdb_ann$type[p_ind - 1] == '(' & wfdb_ann$type[p_ind + 1] == ')']
  
  # For every p_ind, create array of integers like: p_onset : p_offset
  p_continuous <- unlist(lapply(p_ind, function(ind) {
    c(wfdb_ann$sample[(ind - 1)],wfdb_ann$sample[(ind + 1)])
  }))
  
  output[p_continuous] <- 1
  
  
  # QRS -------------------------- 
  # Find QRS indices, then verify onset and offset
  qrs_ind <- which(wfdb_ann$type == 'N')
  qrs_ind <- qrs_ind[wfdb_ann$type[qrs_ind - 1] == '(' & wfdb_ann$type[qrs_ind + 1] == ')']
  
  # For every qrs_ind, create array of integers like: qrs_onset : qrs_offset
  qrs_continuous <- unlist(lapply(qrs_ind, function(ind) {
    c(wfdb_ann$sample[(ind - 1)],wfdb_ann$sample[(ind + 1)])
  }))
  
  output[qrs_continuous] <- 2
  
  
  # T wave --------------------------  
  # Find T indices, then verify onset and offset
  t_ind <- which(wfdb_ann$type == 't')
  t_ind <- t_ind[wfdb_ann$type[t_ind - 1] == '(' & wfdb_ann$type[t_ind + 1] == ')']
  
  # For every t_ind, create array of integers like: t_onset : t_offset
  t_continuous <- unlist(lapply(t_ind, function(ind) {
    c(wfdb_ann$sample[(ind - 1)],wfdb_ann$sample[(ind + 1)])
  }))
  
  output[t_continuous] <- 3
  
  
  return(output)
  
}

# ann_wfdb2continuous2 ----------------------------------------------------
#' @description format [1 1 1 1 0 0 0 0 2 2 2]
#' @param object Either 'egm' or 'annotation_table'. If annotation table, should specify 'length'
#' 
#' @param length Length of corresponding ECG signal. Not necessary if 'object' is of class 'egm'. Assumed to be 5000.
ann_wfdb2continuous2 <- function(object, length=5000) {
  
  if (any(class(object) == 'egm')) {
    length <- nrow(object$signal)
    wfdb_ann <- object$annotation
  # } else if (any(class(object) == 'annotation_table')) {
  #   wfdb_ann <- object
  } else {
    wfdb_ann <- object
  }
  
  if (nrow(object) == 0) {
    output <- rep(NA,length)
    return(output)
    warning('Annotation file is empty')
    break
  }
  
  output <- rep(0,length)
  
  
  # Find P indices, then verify onset and offset
  p_ind <- which(wfdb_ann$type == 'p')
  p_ind <- p_ind[wfdb_ann$type[p_ind - 1] == '(' & wfdb_ann$type[p_ind + 1] == ')']
  
  # For every p_ind, create array of integers like: p_onset : p_offset
  p_continuous <- unlist(lapply(p_ind, function(ind) {
    wfdb_ann$sample[(ind - 1)]:wfdb_ann$sample[(ind + 1)]
  }))
  
  output[p_continuous] <- 1
  
  
  # Find QRS indices, then verify onset and offset
  qrs_ind <- which(wfdb_ann$type == 'N')
  qrs_ind <- qrs_ind[wfdb_ann$type[qrs_ind - 1] == '(' & wfdb_ann$type[qrs_ind + 1] == ')']
  
  # For every qrs_ind, create array of integers like: qrs_onset : qrs_offset
  qrs_continuous <- unlist(lapply(qrs_ind, function(ind) {
    wfdb_ann$sample[(ind - 1)]:wfdb_ann$sample[(ind + 1)]
  }))
  
  output[qrs_continuous] <- 2
  
  
  # Find T indices, then verify onset and offset
  t_ind <- which(wfdb_ann$type == 't')
  t_ind <- t_ind[wfdb_ann$type[t_ind - 1] == '(' & wfdb_ann$type[t_ind + 1] == ')']
  
  # For every t_ind, create array of integers like: t_onset : t_offset
  t_continuous <- unlist(lapply(t_ind, function(ind) {
    wfdb_ann$sample[(ind - 1)]:wfdb_ann$sample[(ind + 1)]
  }))
  
  output[t_continuous] <- 3
  
  
  return(output)
}


# ann_wfdb2continuous3 ----------------------------------------------------
#' @description format [1 2 2 3 0 0 0 0 4 5 5 6]
ann_wfdb2continuous3 <- function(object) {
  if (any(class(object) == 'egm')) {
    length <- nrow(object$signal)
    wfdb_ann <- object$annotation
  } else if (any(class(object) == 'annotation_table')) {
    wfdb_ann <- object
  }
  
  output <- rep(0,length)
  
  
  # Find P indices, then verify onset and offset
  p_ind <- which(wfdb_ann$type == 'p')
  p_ind <- p_ind[wfdb_ann$type[p_ind - 1] == '(' & wfdb_ann$type[p_ind + 1] == ')']
  
  # For every p_ind, create array of integers like: p_onset : p_offset
  p_continuous <- unlist(lapply(p_ind, function(ind) {
    wfdb_ann$sample[(ind - 1)]:wfdb_ann$sample[(ind + 1)]
  }))
  # Note P onset and offset values
  p_onset <- wfdb_ann$sample[p_ind-1]
  p_offset <- wfdb_ann$sample[p_ind+1]
  
  output[p_continuous] <- 2
  output[p_onset] <- 1
  output[p_offset] <- 3
  
  
  # Find QRS indices, then verify onset and offset
  qrs_ind <- which(wfdb_ann$type == 'N')
  qrs_ind <- qrs_ind[wfdb_ann$type[qrs_ind - 1] == '(' & wfdb_ann$type[qrs_ind + 1] == ')']
  
  # For every qrs_ind, create array of integers like: qrs_onset : qrs_offset
  qrs_continuous <- unlist(lapply(qrs_ind, function(ind) {
    wfdb_ann$sample[(ind - 1)]:wfdb_ann$sample[(ind + 1)]
  }))
  # Note QRS onset and offset values
  qrs_onset <- wfdb_ann$sample[qrs_ind-1]
  qrs_offset <- wfdb_ann$sample[qrs_ind+1]
  
  output[qrs_continuous] <- 5
  output[qrs_onset] <- 4
  output[qrs_offset] <- 6
  
  
  # Find T indices, then verify onset and offset
  t_ind <- which(wfdb_ann$type == 't')
  t_ind <- t_ind[wfdb_ann$type[t_ind - 1] == '(' & wfdb_ann$type[t_ind + 1] == ')']
  
  # For every t_ind, create array of integers like: t_onset : t_offset
  t_continuous <- unlist(lapply(t_ind, function(ind) {
    wfdb_ann$sample[(ind - 1)]:wfdb_ann$sample[(ind + 1)]
  }))
  # Note T onset and offset values
  t_onset <- wfdb_ann$sample[t_ind-1]
  t_offset <- wfdb_ann$sample[t_ind+1]
  
  output[t_continuous] <- 8
  output[t_onset] <- 7
  output[t_offset] <- 9
  
  return(output)
  
}


# Filter ------------------------------------------------------------------
ecg_filter <- function(signal, frequency = 500, low = 0.5, high = 40) {
  library(signal)
  if (is.vector(signal)) {
    signal <- matrix(signal, nrow = 1)
  }
  
  # Apply the filter row-wise
  t(apply(signal, 1, function(row) {
    butterFilter <- butter(4, c(low / (frequency / 2), high / (frequency / 2)), type = "pass")
    
    # Eliminate phase distortion
    filtfilt(butterFilter, row)
    
  }))
}


# Difference/Derivative functions -----------------------------------------
add_derivs <- function(signal, number_of_derivs = 2,mask_value=0) {
  if (number_of_derivs == 0) {
    return(signal)
    break
  }
  # If input is of size 1 x 5000, change to vector:
  if (length(dim(signal)) > 1) {
    signal <- c(signal)
    
  }
  
  # Do first derivative manually
  diff_signal <- cbind(signal, c(mask_value, diff(signal)))
  
  if (number_of_derivs > 1) {
    for (i in 2:number_of_derivs) {
      last_deriv <- c(mask_value, diff(diff_signal[, ncol(diff_signal)]))
      diff_signal <- cbind(diff_signal, last_deriv)
    }
  }
  
  return(diff_signal)
}


# Confusion matrix analysis ----------------------------------------------------------------
confusion_analysis <- function(predictions = predictions_integer,
                               actual = testing_annotations) {
  # Quick analysis of predicted values
  
  library(caret)
  levels <- sort(unique(as.vector(actual))) # number of value types must be the same when factoring
  
  confusion_matrix <- confusionMatrix(
    data      = factor(predictions, levels = levels),
    reference = factor(actual, levels = levels)
  )
  
  confusion_table <- t(t(confusion_matrix$table) / colSums(confusion_matrix$table)) * 100
  
  print(round(confusion_matrix$byClass,3))
  print(round(confusion_table, 1))
  # return(confusion_matrix)
  
  output <- list(round(confusion_matrix$byClass,2), round(confusion_table,1))
  names(output) <- c('byClass','table')
  
  # output <- list(confusion_matrix,confusion_table)
  names(output) <- c('matrix','table')
  return(output)
  
}


# LUDB wave analysis --------------------------------------------------------

# Count number of P, QRS, T waves in each sample, cross compare
# Output: which samples / number of samples which do not have the right number of waves

# For each wave, compare onset / offset times
check_ann_prog_ludb <- function(predicted, ludb) {
  # Check annotation progression against LUDB (ie against a ground truth)
  
  # Input: predicted and LUDB annotations. Either matrix form or LUDB
  # Psuedocode:
  #   move into wfdb format as needed
  #   **use LUDB R-R intervals as ground truth reference markers**
  #   define a tolerance window- different for QRS vs p/t...
  #     Could do a tolerance test for both onset offset, 
  #     or could fine the midpoint between them, then compare length?
  
  #   verify predicted QRS are within tolerance 
  #     also verify number of QRS are the same 
  #     if extra QRS outside the tolerance, note. 
  #     if multiple QRS within the tolerance, note
  
  #   Check p and t waves
  #     for each R-R: 
  #       verify correct number and p/t-waves. If extra/lacking, note
  #       verify they're within tolerance. If outside tolerance, note
  #       should do both tasks for indicies prior to 1st R, after last R too
  
  # If annotation is in matrix format (not WFDB), convert to wfdb
  if (class(predicted) == 'numeric') {
    predicted <- ann_continuous2wfdb(predicted)
  }
  
  if (class(ludb) == 'numeric') {
    ludb <- ann_continuous2wfdb(ludb)
  }
  
  # To compare ground truths, remove predicted annotations prior to first LUDB annotation, after last LUDB annotation
  
  
  frequency <- 500
  # QRS tolerance: avg duration = 120 ms. --> tolerance = 60 ms
  qrs_tolerance <- 60 / 1000 * frequency
  p_tolerance <- 60 / 1000 * frequency
  t_tolerance <- 60 / 1000 * frequency
  
  
  tolerance <- c(qrs_tolerance,p_tolerance,t_tolerance)
  waves <- c('N','p','t')
  # Initialize points counter (for each wave type)
  points <- array(0,3)
  duplicate_points <- array(0,3)
  
  # For each wave type (QRS/P/T)
  for (wave in 1:length(waves)) {
    # Define a time stamp for a wave via the mean of the onset and offset
    pred_wave <- unlist(lapply(which(predicted$type == waves[wave]), function(idx) {
      round(mean(predicted$sample[idx + 1], predicted$sample[idx + 1]))
    }))
    
    ludb_wave <- unlist(lapply(which(ludb$type == waves[wave]), function(idx) {
      round(mean(ludb$sample[idx + 1], ludb$sample[idx + 1]))
    }))
    
    # Initialize a vector to track paired ground truth indices
    paired_ludb <- rep(FALSE, length(ludb_wave)) # FALSE means unpaired initially

    
    # Loop through each predicted QRS index
    for (i in seq_along(pred_wave)) {
      pred_time <- pred_wave[i]
      
      # Find all ground truth indices within the tolerance range
      valid_pairs <- which(abs(ludb_wave - pred_time) <= qrs_tolerance)
      
      if (length(valid_pairs) == 0) {
        # No valid pairing found, mark as unpaired
        points[wave] <- points[wave] + 1
        
        # If wave is duplicated, add to separate count:
        # Find index in original wfdb dataframe:
        row <- which(predicted$type == waves[wave])[i]
        # Safeguard for undefined indices
        valid_row_minus_3 <- (row - 3 >= 1) 
        valid_row_plus_3 <- (row + 3 <= nrow(predicted))
        # Modified condition with checks for valid indices
        if ((valid_row_minus_3 && predicted$type[row - 3] == waves[wave]) ||
            (valid_row_plus_3 && predicted$type[row + 3] == waves[wave])) {
          duplicate_points[wave] <- duplicate_points[wave] + 1
        }
        
        
      } else {
        # Check if any of the valid pairs are already paired
        paired_indices <- which(paired_ludb[valid_pairs])
        unpaired_indices <- setdiff(valid_pairs, paired_indices)
        
        if (length(unpaired_indices) > 0) {
          # Pair to the first unpaired ground truth event
          paired_ludb[unpaired_indices[1]] <- TRUE
        } else {
          # If all valid pairs are already paired, add +1 point
          points[wave] <- points[wave] + 1
          
          # If wave is duplicated, add to separate count:
          # Find index in original wfdb dataframe:
          row <- which(predicted$type == waves[wave])[i]
          # Safeguard for undefined indices
          valid_row_minus_2 <- (row - 2 >= 1) 
          valid_row_plus_2 <- (row + 2 <= nrow(predicted))
          # Modified condition with checks for valid indices
          if ((valid_row_minus_2 && predicted$type[row - 2] == waves[wave]) ||
              (valid_row_plus_2 && predicted$type[row + 2] == waves[wave])) {
            duplicate_points[wave] <- duplicate_points[wave] + 1
          }
          
        }
      }
    }
    
    # Unneeded- paired_ludb
    # # Count additional pairings to the same ground truth
    # points[wave] <- points[wave] + sum(paired_ludb > 1)
    
  }
  
    
  points_df <- data.frame(
    qrs_unpair = points[1] - duplicate_points[1], 
    p_unpair = points[2] - duplicate_points[2], 
    t_unpair = points[3] - duplicate_points[3],
    qrs_dupl = duplicate_points[1], 
    p_dupl = duplicate_points[2], 
    t_dupl = duplicate_points[3]
  )

  return(points_df)
  
}

check_ann_prog_ludb_helper <- function(testing_annotations = testing_annotations, 
                                       predictions_integer = predictions_integer) {
  table <- c()
  for (sample in 1:nrow(testing_annotations)) {
    table <- rbind(table,check_ann_prog_ludb(ludb = testing_annotations[sample,], 
                                               predicted = predictions_integer[sample,])) 
  }
  
  samples_w_missed_wave <- data.frame(
    p = round((sum(table$p_unpair) + sum(table$p_dupl)) / nrow(table), 3),
    qrs = round((sum(table$qrs_unpair) + sum(table$qrs_dupl)) / nrow(table), 3),
    t = round((sum(table$t_unpair) + sum(table$t_dupl)) / nrow(table), 3)
  )
  
  missed_waves_per_sample <- data.frame(
    p = round(length(which(table$p_unpair != 0 | table$p_dupl != 0)) / nrow(table), 3),
    qrs = round(length(which(table$qrs_unpair != 0 |table$qrs_dupl != 0)) / nrow(table), 3),
    t = round(length(which(table$t_unpair != 0 | table$t_dupl != 0)) / nrow(table), 3)
  )
  
  
  output <- list(samples_w_missed_wave, missed_waves_per_sample)
  names(output) <- c('samples_w_missed_wave','missed_waves_per_sample')
  print(output)
  return(output)
}

ann_continuous2wfdb <- function(annotations, Fs = 500) {
  # Not up to date, largely defunct for now. Used to convert index to index
  # annotation (value for each time point) to wave onset/peak/offset format. 
  # WFDB format
  
  if (sum(is.na(annotations)) == length(annotations)) {
    annotation_table <- data.frame(time = c(), 
                                   sample = c(), 
                                   type = c(), 
                                   subtype = c(), 
                                   channel = c(), 
                                   number = c()
                                   )
    return(annotation_table)
    break
  }
  
  if (any(annotations == 1)) {
    p_waves <- which(annotations == 1)
    p_change <- ((p_waves[-1] - p_waves[1:(length(p_waves) - 1)]))
    new_pwave <- c(1, ((which(p_change != 1) + 1)))
    new_pwave <- c(new_pwave, length(p_waves) + 1) # account for last value
    
    p_on <- c()
    p_off <- c()
    p_max <- c()
    
    for (i in 1:(length(new_pwave) - 1)) {
      p_on <- c(p_on, p_waves[new_pwave[i]])
      p_off <- c(p_off, p_waves[new_pwave[i + 1] - 1])
      p_max <- c(p_max, round(mean(c(p_on[i], p_off[i])))) # temporary, will need to find peak point
    }
  }
    
  if (any(annotations == 2)) {
    qrs_waves <- which(annotations == 2)
    qrs_change <- ((qrs_waves[-1] - qrs_waves[1:(length(qrs_waves) - 1)]))
    new_qrswave <- c(1, ((which(qrs_change != 1) + 1)))
    new_qrswave <- c(new_qrswave, length(qrs_waves) + 1) # account for last value
    
    qrs_on <- c()
    qrs_off <- c()
    qrs_max <- c()
    
    for (i in 1:(length(new_qrswave) - 1)) {
      qrs_on <- c(qrs_on, qrs_waves[new_qrswave[i]])
      qrs_off <- c(qrs_off, qrs_waves[new_qrswave[i + 1] - 1])
      qrs_max <- c(qrs_max, round(mean(c(
        qrs_on[i], qrs_off[i]
      )))) # temporary, will need to find peak point
    }
  }
  
  if (any(annotations == 3)) {
    t_waves <- which(annotations == 3)
    t_change <- ((t_waves[-1] - t_waves[1:(length(t_waves) - 1)]))
    new_twave <- c(1, ((which(t_change != 1) + 1)))
    new_twave <- c(new_twave, length(t_waves) + 1) # account for last value
    
    t_on <- c()
    t_off <- c()
    t_max <- c()
    
    for (i in 1:(length(new_twave) - 1)) {
      t_on <- c(t_on, t_waves[new_twave[i]])
      t_off <- c(t_off, t_waves[new_twave[i + 1] - 1])
      t_max <- c(t_max, round(mean(c(t_on[i], t_off[i])))) # temporary, will need to find peak point
    }
  }
  
  t <- (1 : length(annotations) ) / Fs
  
  
  # Build Table
  type <- c()

  sample <- c()
  if (exists('p_on')) {
    type <- c(type,array(c("(","p",")"),length(p_on)*3))
    for (i in 1:length(p_on)) {
      sample <- c(sample, p_on[i], p_max[i], p_off[i])
    }
  }
  if (exists('qrs_on')) {
    type <- c(type,array(c("(","N",")"),length(qrs_on)*3))
    for (i in 1:length(qrs_on)) {
      sample <- c(sample, qrs_on[i], qrs_max[i], qrs_off[i])
    }
  }
  
  if (exists('t_on')) {
    type <- c(type,array(c("(","t",")"),length(t_on)*3))
    for (i in 1:length(t_on)) {
      sample <- c(sample, t_on[i], t_max[i], t_off[i])
    }
  }
  
  time <- t[sample]
  subtype <- array(0,length(sample))
  channel <- array(0,length(sample))
  number <- array(0,length(sample))
  
  annotation_table <- data.frame(
    time = time,
    sample = sample,
    type = type,
    subtype = subtype,
    channel = channel,
    number = number)
  
  class(annotation_table) <- c('data.frame','annotation_table')
  
  if (nrow(annotation_table > 0)) {
    annotation_table <- annotation_table[order(annotation_table$sample), ]
    row.names(annotation_table) <- NULL
  }
  
  return(annotation_table)
}



# UIH relative analysis ---------------------------------------------------
# ie- functions when we don't have an answer key

# Check if waves go: P > QRS > T > P...
check_ann_prog <- function(annotation) {
  ann <- ann_continuous2wfdb(annotation)
  
  # Remove first and last 0.5s
  ann <- ann[(ann$sample <= 4750) & (ann$sample >= 250),]
  
  progression <- ann$type[ann$type %in% c('p','N','t')]
  
  # Create logical checks
  check_p <- which(progression == "p")  # Indices of "p"
  check_N <- which(progression == "N")  # Indices of "N"
  check_t <- which(progression == "t")  # Indices of "t"
  
  # Ensure "p" is followed by "N"
  valid_pN <- all(progression[check_p[check_p < length(progression)] + 1] == "N")
  
  # Ensure "N" is followed by "t"
  valid_Nt <- all(progression[check_N[check_N < length(progression)] + 1] == "t")
  
  # Ensure "t" is followed by "p"
  valid_tp <- all(progression[check_t[check_t < length(progression)] + 1] == "p")
  
  # Final check
  if (valid_pN & valid_Nt & valid_tp & length(progression > 1)) {
    # Check if statements, and check the annotation has more than 1 wave
    
    # print("The progression follows the expected pattern.")
    output <- 1
  } else {
    # print("The progression does NOT follow the expected pattern.")
    output <- 0
  }
  return(output)
}

# If there are two waves of the same type (ie pwave) which are within 10 indicies,
# and those indicies are labelled as 0 (non-waves), convert the 0s to pwave values
fill_wave_gaps <- function(annotation = predictions_integer, max_gap = 10) {
  # Find nonzero indices
  nonzero_indices <- which(annotation != 0)
  
  # Iterate through the nonzero indices
  for (i in seq_along(nonzero_indices)) {
    if (i == length(nonzero_indices)) break  # Stop if at last index
    
    current_val <- annotation[nonzero_indices[i]]
    next_idx <- nonzero_indices[i + 1]
    
    # Check if the next wave is of the same value and gap is ≤ max_gap
    if (annotation[next_idx] == current_val && (next_idx - nonzero_indices[i] <= max_gap)) {
      annotation[(nonzero_indices[i] + 1):(next_idx - 1)] <- current_val
    }
  }
  
  return(annotation)
}

check_ann_prog_RPeaks <- function(signal,
                                  predictions_integer, 
                                  pr_tolerance = 200,
                                  rt_tolerance = 170,
                                  terminal_boundary = 250
) {
  # pr_tolerance <- 0.2 * 500 * 2 # 0.2 sec, fudge factor of 100%
  # rt_tolerance <- 0.170 * 500 * 2 # ST: 0.120 sec, RS: ~0.05 sec, fudge factor of 100%
  
  Rpeaks_lower_bound <- terminal_boundary
  Rpeaks_upper_bound <- 5000 - terminal_boundary
  
  # Convert predictions_integer sample to wfdb format
  ann <- ann_continuous2wfdb(predictions_integer)
  sig <- signal
  
  ## Identify and match QRS: ##
  
  Rpeaks <- EGM::detect_QRS(ecg_filter(sig), frequency = 500)
  # Remove Rpeaks if they are too close to the beginning/end
  Rpeaks <- Rpeaks[Rpeaks > Rpeaks_lower_bound &
                     Rpeaks < Rpeaks_upper_bound]
  
  # Find wave indicies/onsets/offsets
  QRS_wave_indices <- which(ann$type == "N")
  QRS_wave_onsets <- ann$sample[QRS_wave_indices - 1]
  QRS_wave_offsets <- ann$sample[QRS_wave_indices + 1]
  p_indices <- which(ann$type == "p")
  p_wave_onsets <- ann$sample[p_indices - 1]
  t_indices <- which(ann$type == "t")
  t_wave_onsets <- ann$sample[t_indices - 1]
  
  # If QRS onset is <250, offset > 4750, remove
  if (any(QRS_wave_onsets < Rpeaks_lower_bound)) {
    logic <- QRS_wave_onsets < Rpeaks_lower_bound
    QRS_wave_onsets <- QRS_wave_onsets[!logic]
    QRS_wave_offsets <- QRS_wave_offsets[!logic]
  }
  if (any(QRS_wave_offsets > Rpeaks_upper_bound)) {
    logic <- QRS_wave_offsets > Rpeaks_upper_bound
    QRS_wave_onsets <- QRS_wave_onsets[!logic]
    QRS_wave_offsets <- QRS_wave_offsets[!logic]
  }
  
  matched_Rpeaks <- list()
  unmatched_Rpeaks <- Rpeaks  # Start with all Rpeaks, remove matched ones later
  unmatched_QRS_waves <- data.frame(onset = integer(), offset = integer())
  
  # Find the matched/unmatched QRS/Rpeaks
  for (i in seq_along(QRS_wave_onsets)) {
    # Find Rpeaks within onset/offset range
    contained <- Rpeaks[Rpeaks >= QRS_wave_onsets[i] &
                          Rpeaks <= QRS_wave_offsets[i]]
    
    if (length(contained) > 0) {
      matched_Rpeaks[[i]] <- contained
      unmatched_Rpeaks <- setdiff(unmatched_Rpeaks, contained)  # Remove matched Rpeaks
    } else {
      unmatched_QRS_waves <- rbind(
        unmatched_QRS_waves,
        data.frame(onset = QRS_wave_onsets[i], offset = QRS_wave_offsets[i])
      )
    }
  }
  
  # matched_waves <- data.frame(Rpeaks = unlist(matched_Rpeaks))
  matched_waves <- data.frame(Rpeaks = unlist(Rpeaks))
  matched_waves$p <- NA
  matched_waves$t <- NA
  
  # Find missed pwaves, twaves
  
  for (i in seq_along(matched_waves$Rpeaks)) {
    peak <- matched_waves$Rpeaks[i]
    
    # Check for P-wave onset within pr_tolerance
    p_match <- p_wave_onsets[p_wave_onsets >= (peak - pr_tolerance) &
                               p_wave_onsets < peak]
    
    # If there are multiple matches, pick the closest one
    if (length(p_match >= 1)) {
      matched_waves$p[i] <- max(p_match)
    }
    
    # Check for T-wave onset within rt_tolerance
    t_match <- t_wave_onsets[t_wave_onsets <= (peak + rt_tolerance) &
                               t_wave_onsets > peak]
    
    # If there are multiple matches, pick the closest one
    if (length(t_match >= 1)) {
      matched_waves$t[i] <- min(t_match)
    }
  }
  
  # Find duplicated pwaves, twaves
  duplicate_pwaves <- c()
  duplicate_twaves <- c()
  
  for (i in seq_along(matched_waves$Rpeaks)) {
    if (i > 1) {
      # Skip first peak
      prev_peak <- matched_waves$Rpeaks[i - 1]
      curr_peak <- matched_waves$Rpeaks[i]
      pwaves <- p_wave_onsets[p_wave_onsets > prev_peak &
                                p_wave_onsets < curr_peak]
      
      if (length(pwaves) > 1) {
        unmatched_p <- setdiff(pwaves, matched_waves$p[i])
        duplicate_pwaves <- c(duplicate_pwaves, unmatched_p)
      }
    }
    
    if (i < nrow(matched_waves)) {
      # Skip last peak
      curr_peak <- matched_waves$Rpeaks[i]
      next_peak <- matched_waves$Rpeaks[i + 1]
      twaves <- t_wave_onsets[t_wave_onsets > curr_peak &
                                t_wave_onsets < next_peak]
      
      if (length(twaves) > 1) {
        unmatched_t <- setdiff(twaves, matched_waves$t[i])
        duplicate_twaves <- c(duplicate_twaves, unmatched_t)
      }
    }
  }
  
  output <- data.frame(
    missed_QRS = length(unmatched_Rpeaks) / nrow(matched_waves),
    missed_P = sum(is.na(matched_waves$p)) / nrow(matched_waves),
    missed_T = sum(is.na(matched_waves$t)) / nrow(matched_waves),
    duplicate_QRS = nrow(unmatched_QRS_waves) / nrow(matched_waves),
    duplicate_P = length(duplicate_pwaves) / nrow(matched_waves),
    duplicate_T = length(duplicate_twaves) / nrow(matched_waves)
  )
  
  return(output)
  # print(unmatched_Rpeaks)
  # print(unmatched_QRS_waves)
  # print(duplicate_pwaves)
  # print(duplicate_twaves)
  # matched_waves
  
  # Pseudocode:  
  # Need to adjust for waves which occur too close to beginning/end: 
  #   eliminate Rpeaks if <250, >4750
  #   there should be one labelled p-wave onset prior to  QRS (earliest: index ~150)
  #     ignore waves prior to this p-wave
  #   there should be one labelled t-wave onset after the QRS
  #     ignore waves after this t-wave
  
  # For "normal" cases not at the beginning/end:
  #   Check if there is p-wave before, t-wave after
  #   Note if there are duplicates:
  #     more than 1 p-wave before the prior QRS
  #     more than 1 t-wave after  the QRS
  #   Note if there are unmatched:
  #     no p-wave prior to QRS
  #     no t-wave after QRS
  
  # Toggle based on Rpeak- so it still looks for a p-wave even if QRS is missed, for example
  
  # Outputs:
  #   % missed P, QRS, T
  #   % duplicated P, QRS, T
  
}
# Custom plotting function ------------------------------------------------
plot_func <- function(y, 
                      color = 0,
                      linewidth=0.5, 
                      pointsize = 1.5, 
                      ylim = NULL, 
                      plotly = 'yes', 
                      x) {
  # Custom plotting function, Can toggle between ggplot (plotly = 'no') vs
  # plotly 
  library(ggplot2)
  library(plotly)
  color <- c(color)
  y <- c(y)
  
  color[color == 1] <- 'p'
  color[color == 2] <- 'N'
  color[color == 3] <- 't'
  
  x <- 1:length(y)
  
  frame <- data.frame(Time = x, Signal = y)
  plot <-
    ggplot(frame, aes(Time, Signal, color = color)) +
    geom_path(linewidth = linewidth, aes(group = 1)) + geom_point(size = pointsize) +
    scale_x_continuous(breaks = seq(0, 10, 1)) + 
    theme(legend.position = "none") + 
    theme(legend.title = element_blank()) +
    theme_bw() + 
    coord_cartesian(ylim = ylim)
  
  if (plotly == 'yes') {
    plot <- ggplotly(plot)
  }
  
  # original linewidth = 0.25
  return(plot)
}


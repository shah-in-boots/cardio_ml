# Predict -----------------------------------------------------------------
library(keras)
library(tools)
library(EGM)
source('annotator_prep_functions.R')
load("/mmfs1/projects/cardio_darbar_chi/common/cohorts/wes_ml/annotator_models/uih_ecgs.RData")
load('../models/model_log.RData')

row <- 15
# for (row in 11:nrow(model_log)) {
lead <- 1

# Note: uih_ecgs dataset was generated via generate_uih_test_ecgs.R. Sinus ECGs only

signal <- do.call(rbind,lapply(1:length(uih_ecgs), function(idx) uih_ecgs[[idx]]$signal[[lead+1]]))

model <- load_model_tf(paste0('../models/',model_log$name[row]))

filtered <- array(0,dim(signal))
for (i in 1:nrow(signal)) {
  filtered[i,] <- ecg_filter(signal[i,])
}

# Predict
predictions <- model %>% predict(filtered)
predictions_integer <- array(0,c(nrow(predictions),ncol(predictions))) 
for (i in 1: nrow(predictions)) {
  predictions_integer[i,] <- max.col(predictions[i,,])
}
#convert from dimension value 1,2,3,4 to 0,1,2,3
predictions_integer <- predictions_integer - 1


# Stats -------------------------------------------------------------------
prog_count <- array(NA,nrow(predictions_integer))
prog_count_revised <- array(NA,nrow(predictions_integer))

predictions_integer_revised <- array(0,dim(predictions_integer))
for (i in 1:nrow(predictions_integer)) {
  predictions_integer_revised[i,] <- fill_wave_gaps(predictions_integer[i,],20)
}

for (i in 1:nrow(predictions_integer)) {
  prog_count[i] <- check_ann_prog(annotation = predictions_integer[i,])
  prog_count_revised[i] <- check_ann_prog(predictions_integer_revised[i,])
}

fraction <- sum(prog_count == 0) / length(prog_count)
revised_fraction <- sum(prog_count_revised == 0) / length(prog_count_revised)

# Find the average probability of each time point
confidence <- round(mean(apply(predictions, c(1, 2), max)),4)

# model_log$confidence[row] <- confidence
# model_log$uih_prog[row] <- fraction
# model_log$uih_prog_revised[row] <- revised_fraction
# print(paste(row))
# }


# Compare to detect_QRS() -------------------------------------------------
library(EGM)
idx <- 1
# For a given sample

for (idx in 1:nrow(signal)) {
  
  Rpeaks <- EGM::detect_QRS(ecg_filter(signal[idx, ]), frequency = 500)
  # Convert predictions_integer sample to wfdb format
  # Remove annotations and Rpeaks which occur <250 and >4750
  
  # Match Rpeaks which fall within the ML QRS interval. Check for duplicates / unmatched
  
  # Of the Rpeaks, check for adjacent P and T waves
  
  
  # For pwaves:
  #   Is there a pwave prior to the R peak (within tolerance)
  #     What is the PR interval of that wave
}

# idx <- idx+1
annotations <- array(0,5000)
annotations[EGM::detect_QRS(ecg_filter(signal[idx,]), frequency = 500)] <- 1 
plot_func(ecg_filter(signal[idx,]),annotations)


# plot --------------------------------------------------------------------
sample <- 1
raw <- plot_func(signal[sample,],predictions_integer[sample,])
filt <- plot_func(filtered[sample,],predictions_integer[sample,])

ann <- ann_continuous2wfdb(predictions_integer[sample,])
ann[ann$type %in% c('p','N','t'),]

filt
# subplot(raw,filt,nrows = 2)

# Ideas for determining accuracy: compare to find_Rpeaks, and check p > QRS > t > p 
#   Rpeaks anchors the QRS (global), and p qrs t checks local progression

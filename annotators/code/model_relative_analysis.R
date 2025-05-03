# Predict -----------------------------------------------------------------
library(keras)
library(tools)
library(EGM)
source('annotator_prep_functions.R')
load("/mmfs1/projects/cardio_darbar_chi/common/cohorts/wes_ml/annotator_models/uih_ecgs.RData")
load('../models/model_log.RData')

row <- 12
# for (row in 11:nrow(model_log)) {
lead <- 1

# Note: uih_ecgs dataset was generated via generate_uih_test_ecgs.R. Sinus ECGs only

signal <- do.call(rbind,lapply(1:length(uih_ecgs), function(idx) uih_ecgs[[idx]]$signal[[lead+1]]))

# for LOCAL testing only. Suppress/delete when on cluster:
dir <- '../../AF_DM/wfdb/matched_controls/'
testfiles <- file_path_sans_ext(list.files(path = dir, pattern = '.hea')[1:100])
signal <- do.call(rbind, lapply(testfiles, function(file) {
  read_signal(record = file, record_dir = dir)$I
}))
# ---


model <- load_model_tf(paste0('../models/',model_log$name[row]))

filtered <- array(0,dim(signal))
for (i in 1:nrow(signal)) {
  filtered[i,] <- ecg_filter(signal[i,])
}

# Predict
uih_predictions <- model %>% predict(filtered)
uih_predictions_integer <- array(0,c(nrow(uih_predictions),ncol(uih_predictions))) 
for (i in 1: nrow(uih_predictions)) {
  uih_predictions_integer[i,] <- max.col(uih_predictions[i,,])
}
#convert from dimension value 1,2,3,4 to 0,1,2,3
uih_predictions_integer <- uih_predictions_integer - 1


# Stats -------------------------------------------------------------------
prog_count <- array(NA,nrow(uih_predictions_integer))
prog_count_revised <- array(NA,nrow(uih_predictions_integer))

uih_predictions_integer_revised <- array(0,dim(uih_predictions_integer))
for (i in 1:nrow(uih_predictions_integer)) {
  uih_predictions_integer_revised[i,] <- fill_wave_gaps(uih_predictions_integer[i,],20)
}

for (i in 1:nrow(uih_predictions_integer)) {
  prog_count[i] <- check_ann_prog(annotation = uih_predictions_integer[i,])
  prog_count_revised[i] <- check_ann_prog(uih_predictions_integer_revised[i,])
}

uih_prog <- sum(prog_count == 0) / length(prog_count)
uih_prog_revised <- sum(prog_count_revised == 0) / length(prog_count_revised)

confidence <- round(mean(apply(uih_predictions, c(1, 2), max)),4)

# model_log$confidence[row] <- confidence
# model_log$uih_prog[row] <- fraction
# model_log$uih_prog_revised[row] <- revised_fraction
# print(paste(row))
# }


# Compare to detect_QRS() -------------------------------------------------
library(dplyr)
library(EGM)
# For a given sample


uih_progression <- do.call(rbind,lapply(1:nrow(signal), function(idx) {
  round(check_ann_prog_RPeaks(signal[idx,], uih_predictions_integer[idx,]),2)
}))

uih_Rpeaks <- do.call(rbind,lapply(1:nrow(signal), function(idx) {
  round(check_ann_prog_RPeaks(signal[idx,], uih_predictions_integer[idx,]),2)
}))

Rpeaks_prog <- data.frame(
  QRS = sum(uih_Rpeaks$missed_QRS + uih_Rpeaks$duplicate_QRS)/nrow(uih_Rpeaks),
  P = sum(uih_Rpeaks$missed_P + uih_Rpeaks$duplicate_P)/nrow(uih_Rpeaks),
  T = sum(uih_Rpeaks$missed_T + uih_Rpeaks$duplicate_T)/nrow(uih_Rpeaks)
)

model_log$Rpeaks_prog_P[row] = Rpeaks_prog$P
model_log$Rpeaks_prog_QRS[row] = Rpeaks_prog$QRS
model_log$Rpeaks_prog_T[row] = Rpeaks_prog$T
model_log$Rpeaks_prot_total[row] = sum(Rpeaks_prog)

# idx <- 92
# ann_plot <- uih_predictions_integer_revised[idx,]
# Rpeaks <- detect_QRS(signal[idx,],500)
# ann_plot[Rpeaks] <- 4
# plot_func(signal[idx,],ann_plot)
# ann_continuous2wfdb(uih_predictions_integer_revised[idx,])


# UIH  ----------------------------------------------------------------


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

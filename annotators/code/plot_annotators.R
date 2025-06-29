library(keras)
source('annotator_prep_functions.R')
load('../models/model_log.RData')
load('../uih_ecgs.RData')


lead <- 12
rows <- c(409,663) # 564

normalize <- TRUE
leads <- c('I','II','III','AVR','AVL','AVF','V1','V2','V3','V4','V5','V6')

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

# Predict
for (row in 1:length(rows)) {
  model <- load_model_tf(paste0('../models/',model_log$name[rows[row]],'.h5'))
  uih_predictions <- model %>% predict(filtered)
  uih_predictions_integer <- array(0, c(nrow(uih_predictions), ncol(uih_predictions)))
  for (i in 1:nrow(uih_predictions)) {
    uih_predictions_integer[i, ] <- max.col(uih_predictions[i, , ])
  }
  #convert from dimension value 1,2,3,4 to 0,1,2,3
  uih_predictions_integer <- uih_predictions_integer - 1
  if (row == 1) {
    predictions1 <- uih_predictions_integer
  }
  
  if (row == 2) {
    predictions2 <- uih_predictions_integer
  }
}
idx=0


# start <- 1
# subplot(lapply(seq(start, start + 5), function(idx) {
#   plot_func(uih_samples[idx, ], uih_predictions_integer[idx, ])
# }), ncol = 2)

idx=idx+1
plot_func(ecg_filter(uih_samples[idx, ]), predictions1[idx, ])
Sys.sleep(0.25)
plot_func(ecg_filter(uih_samples[idx, ]), predictions2[idx, ])
ann <- ann_continuous2wfdb(predictions2[idx, ])
ann[ann$type %in% c('p','N','t'),]
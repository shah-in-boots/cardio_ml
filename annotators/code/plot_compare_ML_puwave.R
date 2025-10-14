# Compare ML annotators to ecgpuwave

# predict -----------------------------------------------------------------
library(keras)
# setwd('C:/Users/darre/OneDrive/Documents/UICOM Research/annotators/code')
source('annotator_prep_functions.R')
load('../models/model_log.RData')
load('../uih_ecgs_supp.RData')


lead <- 12
rows <- c(881) # lead iii: 358, 722
number_of_derivs <- c(2)
location <- 'local'

normalize <- TRUE
leads <- c('I','II','III','AVR','AVL','AVF','V1','V2','V3','V4','V5','V6')

uih_samples <- do.call(rbind, lapply(1:length(uih_ecgs_supp), function(idx)
  uih_ecgs_supp[[idx]]$signal[[leads[lead]]]))

filtered <- array(0, dim(uih_samples))
for (i in 1:nrow(uih_samples)) {
  filtered[i, ] <- ecg_filter(uih_samples[i, ])
}

if (normalize) {
  for (i in 1:nrow(uih_samples)) {
    filtered[i, ] <- (filtered[i, ] - min(filtered[i, ])) / (max(filtered[i, ]) - min(filtered[i, ])) * 100
  }
}


# Predict ML:
input <- array(NA, c(dim(filtered), number_of_derivs + 1))
for (i in 1:nrow(filtered)) {
  input[i, , ] <- add_derivs(signal = filtered[i, ], number_of_derivs = number_of_derivs)
}
  
  
model <- load_model_tf(paste0('../models/', model_log$name[rows], '.h5'))
uih_predictions <- model %>% predict(input)
uih_predictions_integer <- array(0, c(nrow(uih_predictions), ncol(uih_predictions)))
for (i in 1:nrow(uih_predictions)) {
  uih_predictions_integer[i, ] <- max.col(uih_predictions[i, , ])
}
#convert from dimension value 1,2,3,4 to 0,1,2,3
uih_predictions_integer <- uih_predictions_integer - 1
predictions_ML <- uih_predictions_integer


# Prep ecgpuwave annotations
predictions_ecgpuwave <- array(NA,dim(predictions_ML))
for (i in 1:nrow(predictions_ecgpuwave)) {
  predictions_ecgpuwave[i,] <- ann_wfdb2continuous2(uih_ecgs_supp[[i]]$annotation)
}

idx=0


# samples <- model_log$training_samples[rows[2]][[1]]
# sort(samples$samples[samples$dataset == 'uih'])
# 1  4  6  8 10 11 12 14 16 17 18 19 20 21 22 23 24 25
# plot --------------------------------------------------------------------

location <- 'local'

# Automatically opens plot in browser if running on laptop locally
if (location == 'local') {
  idx=idx+1
  plot1 <- plot_func2(ecg_filter(uih_samples[idx, ]), predictions_ML[idx, ],linewidth = 0.25,pointsize = 0.5)
  temp_file1 <- tempfile(fileext = ".html")
  saveWidget(plot1, temp_file1, selfcontained = TRUE)
  browseURL(temp_file1)
  
  if (length(rows) > 1) {
    plot2 <- plot_func(ecg_filter(uih_samples[idx, ]), predictions_ecgpuwave[idx, ],linewidth = 0.25,pointsize = 0.5)
    temp_file2 <- tempfile(fileext = ".html")
    saveWidget(plot2, temp_file2, selfcontained = TRUE)
    Sys.sleep(0.25)
    browseURL(temp_file2)
  }
  
} else {
  idx=idx+1
  plot_func2(c(ecg_filter(uih_samples[idx, ])), c(predictions_ML[idx, ]),linewidth = 0.25,pointsize = 0.5)
  # Sys.sleep(2)
  plot_func(ecg_filter(uih_samples[idx, ]), predictions_ecgpuwave[idx, ],linewidth = 0.25,pointsize = 0.5)
  isoelec_find(signal = uih_samples[idx, ], predictions_ML[idx, ])
}

print(idx)
ann <- ann_continuous2wfdb(predictions2[idx, ])
# ann[ann$type %in% c('p','N','t'),]


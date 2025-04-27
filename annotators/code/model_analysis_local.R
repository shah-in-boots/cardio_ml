# View Performance of model trained on cluster
# Description: largely outdated, as stats found in this script are added to the model_log.RData file automatically now.
# Prep --------------------------------------------------------------------
library(keras)
load('../models/model_log.RData')
load('../ludb_set.RData')
source('annotator_prep_functions.R')

# Specify model of choice from model_log
lead <- 1
row <- 1
mask_value <- -1

zero_prior_sig <- TRUE # set signal prior to first LUDB annotation, after last annotation to 0. Set as 'TRUE' if comparing to LUDB annotations

# Load model
model <- load_model_tf(paste0('../models/',model_log$name[row]))

testing_samples <- 1:200
testing_samples <- testing_samples[-unlist(model_log$training_samples[row])]

# Extract only the testing samples, and the lead (first 'lead' is just index values 1:5000, so add 1)

# Signal:
# Sort signal, put into matrix form, and filter
testing_signal_list <- lapply(testing_samples, function(idx) ludb_set[[idx]]$signal[[lead+1]])
n_rows <- length(testing_signal_list)  
n_cols <- length(unlist(testing_signal_list[1]))
testing_signal <- matrix(NA, nrow = n_rows, ncol = n_cols)
for (i in seq_along(testing_signal_list)) {
  testing_signal[i, ] <- ecg_filter(testing_signal_list[[i]])
}

# Annotations
# Sort annotations
testing_annotations <- t(sapply(testing_samples, function(idx) ann_wfdb2continuous2(ludb_set[[idx]]$annotation[[lead]])))

# If requested, zero out signal and annotations in terminal 1000 indicies (rouhgly before first/last annotations are)
if (zero_prior_sig) {
  testing_signal[,c(1:1000, 4001:5000)] <- mask_value
  testing_annotations[,c(1:1000, 4001:5000)] <- 0
  
  # testing_annotations <- lapply(1:length(testing_annotations), function(idx) {
  #   testing_annotations[[idx]][testing_annotations[[idx]]$sample > 1000 & testing_annotations[[idx]]$sample < 4001]
  # })
}

# Predict
predictions <- model %>% predict(testing_signal)
predictions_integer <- array(0,c(nrow(predictions),ncol(predictions))) 
for (i in 1: nrow(predictions)) {
  predictions_integer[i,] <- max.col(predictions[i,,])
}
#convert from dimension value 1,2,3,4 to 0,1,2,3
predictions_integer <- predictions_integer - 1

# Analysis
confusion <- confusion_analysis(predictions = predictions_integer, actual = testing_annotations)
wave_counter <- check_ann_prog_ludb_helper(predictions_integer = predictions_integer, testing_annotations = testing_annotations)

model_log$wave_counter[row] <- I(list(wave_counter))

# Correct --------------------------------------------------------------------
# NOTE: THIS IS UNECESSARY NOW

# Correct annotations in testing_annotations which had wave labels cut off at around index 1000 and 4000
# for (i in 1:length(testing_annotations)) {
#   ann <- testing_annotations[[i]]
#   end_ind <- nrow(ann)
#   
#   if (ann$type[end_ind] %in% c('p','N','t')) {
#     new_row <- data.frame(
#       time = "00:00:08.000",
#       sample = 4000,
#       type = ")",
#       subtype = 0,
#       channel = 0,
#       number = 0
#     )
#     
#     ann <- rbind(ann,new_row)
#     testing_annotations[[i]] <- ann
#     class(testing_annotations[[i]]) <- c("annotation_table","data.table","data.frame")
#     print(paste(i))
#   }
#   
#   if (ann$type[1] %in% c('p','N','t')) {
#     new_row <- data.frame(
#       time = "00:00:02.000",
#       sample = 1000,
#       type = "(",
#       subtype = 0,
#       channel = 0,
#       number = 0
#     )
#     
#     ann <- rbind(ann,new_row)
#     ann <- ann[order(ann$sample),]
#     testing_annotations[[i]] <- ann
#     class(testing_annotations[[i]]) <- c("annotation_table","data.table","data.frame")
#     print(paste(i))
#   }
# }


# Plot --------------------------------------------------------------------

sample <- 43
truth_plot <- plot_func(testing_signal[sample,],testing_annotations[sample,])
predicted_plot <- plot_func(testing_signal[sample,],predictions_integer[sample,])

subplot(truth_plot,predicted_plot,nrows = 2)

# Add duplicated column: if prior or following wave is the same, add it

# test --------------------------------------------------------------------
output <- c()
for (sample in 1:nrow(testing_annotations)) {
  output <- rbind(output,check_ann_prog_ludb(ludb = testing_annotations[sample,], 
                                             predicted = predictions_integer[sample,])) 
}
(sum(output$qrs_unpair) + sum(output$qrs_dupl)) / nrow(output)
(sum(output$p_unpair) + sum(output$p_dupl)) / nrow(output)
(sum(output$t_unpair) + sum(output$t_dupl)) / nrow(output)

length(which(output$qrs_unpair != 0 | output$qrs_dupl != 0)) / nrow(output)
length(which(output$p_unpair != 0 | output$p_dupl != 0)) / nrow(output)
length(which(output$t_unpair != 0 | output$t_dupl != 0)) / nrow(output)

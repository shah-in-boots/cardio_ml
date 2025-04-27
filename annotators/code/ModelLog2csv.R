# Info: reads the model_log.RData file. Simplifies it to be used in an excel document for easier viewing

library(dplyr)
library(purrr)

load('../models/model_log.RData')

# Remove unneeded columns
model_log <- model_log %>%
  select(-time, -lead, -name, -`ann_style`, -normalize, -training_samples)

# model_log$dilate_range <- model_log$dilate_range[2]

# Simplify dilation range to a single value (the max value of the 2 value range)
model_log$dilate_range <- unlist(lapply(1:nrow(model_log),function(idx) model_log$dilate_range[[idx]][2]))

# Extract sensitivity and PPV values from testing set:
# lapply(1:nrow(model_log),function(idx) model_log$confusion[[idx]]$matrix[,c(1,3)])

model_log <- model_log %>%
  mutate(
    sens_0 = map_dbl(confusion, ~ .x$matrix["Class: 0", "Sensitivity"]),
    sens_p = map_dbl(confusion, ~ .x$matrix["Class: 1", "Sensitivity"]),
    sens_qrs = map_dbl(confusion, ~ .x$matrix["Class: 2", "Sensitivity"]),
    sens_t = map_dbl(confusion, ~ .x$matrix["Class: 3", "Sensitivity"]),
    ppv_0 = map_dbl(confusion, ~ .x$matrix["Class: 0", "Pos Pred Value"]),
    ppv_p = map_dbl(confusion, ~ .x$matrix["Class: 1", "Pos Pred Value"]),
    ppv_qrs = map_dbl(confusion, ~ .x$matrix["Class: 2", "Pos Pred Value"]),
    ppv_t = map_dbl(confusion, ~ .x$matrix["Class: 3", "Pos Pred Value"])
  ) %>%
  select(-confusion)  # Remove wave_counter column

# Extract missed_waves column, specifically the 'missed_waves_per_sample' subcolumn
model_log <- model_log %>%
  mutate(
    missed_p = map_dbl(wave_counter, ~ .x$missed_waves_per_sample[["p"]][1]),
    missed_qrs = map_dbl(wave_counter, ~ .x$missed_waves_per_sample[["qrs"]][1]),
    missed_t = map_dbl(wave_counter, ~ .x$missed_waves_per_sample[["t"]][1])
  ) %>%
  select(-wave_counter)  # Remove wave_counter column

# Add column for the total waves missed
model_log <- model_log %>%
  mutate(missed_total = missed_p + missed_qrs + missed_t)

write.csv(model_log,'../model_log.csv')

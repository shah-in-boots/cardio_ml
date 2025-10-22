# data cleaning ---------------------------------------------------------------
# Info: reads the model_log.RData file. Simplifies it to be used in an excel document for easier viewing

library(dplyr)
library(purrr)

setwd("C:/Users/darre/OneDrive/Documents/UICOM Research/annotators/code")
load('../models/model_log.RData')

# Remove unneeded columns
model_log <- model_log %>%
  select(-time, 
         -name, 
         -ann_style, 
         # -normalize, 
         -training_samples)

model_log <- model_log %>%
  select(type, lead, everything())

leads <- c("i","ii","iii","avr","avl","avf","v1","v2","v3","v4","v5","v6")
model_log$lead <- leads[model_log$lead]

# model_log$dilate_range <- model_log$dilate_range[2]

# Simplify dilation range to a single value (the max value of the 2 value range)
model_log$dilate_range <- unlist(lapply(1:nrow(model_log),function(idx) model_log$dilate_range[[idx]][2]))

# Extract sensitivity and PPV values from testing set:
# lapply(1:nrow(model_log),function(idx) model_log$confusion[[idx]]$matrix[,c(1,3)])

# Define helper function:
library(purrr)
safe_extract <- possibly(
  function(x, cls, metric) x$matrix[cls, metric],
  otherwise = NA_real_
)

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
  
  mutate(
    sens_0_uih   = map_dbl(confusion_uih, ~ safe_extract(.x, "Class: 0", "Sensitivity")),
    sens_p_uih   = map_dbl(confusion_uih, ~ safe_extract(.x, "Class: 1", "Sensitivity")),
    sens_qrs_uih = map_dbl(confusion_uih, ~ safe_extract(.x, "Class: 2", "Sensitivity")),
    sens_t_uih   = map_dbl(confusion_uih, ~ safe_extract(.x, "Class: 3", "Sensitivity")),
    ppv_0_uih    = map_dbl(confusion_uih, ~ safe_extract(.x, "Class: 0", "Pos Pred Value")),
    ppv_p_uih    = map_dbl(confusion_uih, ~ safe_extract(.x, "Class: 1", "Pos Pred Value")),
    ppv_qrs_uih  = map_dbl(confusion_uih, ~ safe_extract(.x, "Class: 2", "Pos Pred Value")),
    ppv_t_uih    = map_dbl(confusion_uih, ~ safe_extract(.x, "Class: 3", "Pos Pred Value"))
  ) %>%
  select(-confusion,-confusion_uih)  # Remove wave_counter column

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
  mutate(missed_total = missed_p + missed_qrs + missed_t) %>%
  
  mutate(f1_0 = round(2 * (sens_0 * ppv_0) / (sens_0 + ppv_0),3)) %>%
  mutate(f1_p = round(2 * (sens_p * ppv_p) / (sens_p + ppv_p),3)) %>%
  mutate(f1_qrs = round(2 * (sens_qrs * ppv_qrs) / (sens_qrs + ppv_qrs),3)) %>%
  mutate(f1_t = round(2 * (sens_t * ppv_t) / (sens_t + ppv_t),3)) %>%
  
  mutate(f1_0_uih = round(2 * (sens_0_uih * ppv_0_uih) / (sens_0_uih + ppv_0_uih),3)) %>%
  mutate(f1_p_uih = round(2 * (sens_p_uih * ppv_p_uih) / (sens_p_uih + ppv_p_uih),3)) %>%
  mutate(f1_qrs_uih = round(2 * (sens_qrs_uih * ppv_qrs_uih) / (sens_qrs_uih + ppv_qrs_uih),3)) %>%
  mutate(f1_t_uih = round(2 * (sens_t_uih * ppv_t_uih) / (sens_t_uih + ppv_t_uih),3))

# Rename columns:
model_log <- model_log %>% rename(p_length = p_wave_length, t_length = t_wave_length)
model_log$p_length <- model_log$p_length / 500 * 1000 # length in ms
model_log$t_length <- model_log$t_length / 500 * 1000 # length in ms

# write.csv(model_log,'../model_log.csv')


# full tables -------------------------------------------------------------
idx=1#+idx
library(formattable)
idx=7
leads <- c("i","ii","iii","avr","avl","avf","v1","v2","v3","v4","v5","v6")
lead <- leads[idx]

full_table <- model_log
full_table$row <- 1:nrow(model_log)

full_table <- full_table %>% filter(lead == !!lead) %>% 
  filter(Rpeaks_prog_total < 0.45) %>%
  filter(uih_prog < 0.45) %>%
  arrange(max_noise, filter, rounds, bilstm_layers, normalize, type, epochs) %>%
  filter(!is.na(uih_prog)) %>%
  select(-lead, -dropout, -filters, -Rpeaks_prog_P, -Rpeaks_prog_QRS, -Rpeaks_prog_T,
         -sens_0,-sens_p,-sens_qrs,-sens_t,-ppv_0,-ppv_p,-ppv_qrs,-ppv_t,
         -sens_0_uih,-sens_p_uih,-sens_qrs_uih,-sens_t_uih,-ppv_0_uih,-ppv_p_uih,-ppv_qrs_uih,-ppv_t_uih,
         -missed_p,-missed_qrs,-missed_t, -dilate_range) %>%
  select(rounds, bilstm_layers, type, max_noise, epochs, 
         f1_0_uih, f1_p_uih, f1_qrs_uih, f1_t_uih,
         f1_0, f1_p, f1_qrs, f1_t,
         row,
         uih_prog, everything())

# # Table (in R)
formattable(full_table, list(
  uih_prog = color_tile("green", "red"),
  Rpeaks_prog_total = color_tile("green", "red"),
  f1_0 = color_tile("red", "green"),
  f1_p = color_tile("red", "green"),
  f1_qrs = color_tile("red", "green"),
  f1_t = color_tile("red", "green"),
  f1_0_uih = color_tile("red", "green"),
  f1_p_uih = color_tile("red", "green"),
  f1_qrs_uih = color_tile("red", "green"),
  f1_t_uih = color_tile("red", "green")
  ))
    # confidence, sens 4x, ppv 4x

idx
leads[idx]


#       25%   50%   75% 
# 0   0.955 0.960 0.965 
# p   0.760 0.827 0.860 
# qrs 0.888 0.910 0.920
# t   0.844 0.874 0.890

  

# clean tables ------------------------------------------------------------
library(dplyr)
library(tidyverse)
library(ggplot2)
library(formattable)
library(scales)
# Define the value sets for each column
max_noise <- c(0.01, 0.03, 0.05)
filter <- c(TRUE, FALSE)
rounds <- c(3, 4)
bilstm_layers <- c(200, 300)
epochs <- c(10, 15, 20, 30, 40)
normalize <- c(TRUE,FALSE)


# run ---------------------------------------------------------------------
# Create a dataframe with all possible combinations, sorted in the specified order
idx <- idx+1
# idx=2
leads <- c("i","ii","iii","avr","avl","avf","v1","v2","v3","v4","v5","v6")
lead <- leads[idx]
lead_table <- expand.grid(
  max_noise = max_noise,
  filter = filter,
  rounds = rounds,
  bilstm_layers = bilstm_layers,
  epochs = epochs,
  normalize = normalize

)



# Perform a left join to match rows between lead_table and model_log
lead_table <- lead_table %>% 
  left_join(model_log %>% filter(lead == !!lead) %>% select(normalize,
                                                            max_noise,
                                                            filter, 
                                                            rounds, 
                                                            bilstm_layers, 
                                                            epochs, 
                                                            uih_prog, 
                                                            Rpeaks_prog_total),
            by = c("max_noise", 
                   "filter", 
                   "rounds", 
                   "bilstm_layers", 
                   "epochs",
                   "normalize"))

lead_table <- lead_table %>%
  arrange(max_noise, filter, rounds, bilstm_layers, epochs,normalize) %>%
  filter(!is.na(uih_prog)) # remove NA rows

# # Table (in R)
# formattable(lead_table, list(
#   uih_prog = color_tile("green", "red"),
#   Rpeaks_prog_total = color_tile("green", "red")
# ))


formattable(lead_table, list(
  uih_prog = formatter("span", 
                       style = ~ style(display = "block",
                                       background.color = colorRampPalette(c("green", "red"))(100)[rescale(uih_prog, from = c(0.1, 0.6), to = c(1, 100))])),
  Rpeaks_prog_total = formatter("span", 
                                style = ~ style(display = "block",
                                                background.color = colorRampPalette(c("green", "red"))(100)[rescale(Rpeaks_prog_total, from = c(0.1, 0.6), to = c(1, 100))])))
)

idx
leads[idx]


# Print best models -------------------------------------------------------


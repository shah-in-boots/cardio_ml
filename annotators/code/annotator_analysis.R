# Inputs -----------------------------------------------------------------
# PPlotting 12-lead analysis, for UICOM Researhc forum 2025
library(keras)
source('annotator_prep_functions.R')
load('../models/model_log.RData')
load('../uih_ecgs_supp.RData')


lead <- 'V6' # lead list: I II III AVR AVL AVF V1 V2 V3 V4 V5 V6

leads <- c('I', 'II', 'III', 'AVR', 'AVL', 'AVF', 'V1', 'V2', 'V3', 'V4', 'V5', 'V6')
models <- c(861,856,851,846,836,841,826,821,866,871,876,881)

number_of_derivs <- 2

output <- setNames(vector("list", (length(leads)+1)), c(leads,'puwave'))
uih_supp_Rpeak_std <- setNames(vector("list", (length(leads)+1)), c(leads,'puwave'))
# rows <- c(881) # lead iii: 358, 722
# Prepare and predict -----------------------------------------------------
for (lead in leads) {
  lead_number <- which(leads == lead)
  rows <- models[lead_number]
  
  
  
  uih_samples <- do.call(rbind, lapply(1:length(uih_ecgs_supp), function(idx)
    uih_ecgs_supp[[idx]]$signal[[lead]]))
  
  # Filter and normalize:
  filtered <- array(0, dim(uih_samples))
  for (i in 1:nrow(uih_samples)) {
    filtered[i, ] <- ecg_filter(uih_samples[i, ])
  }
  for (i in 1:nrow(uih_samples)) {
    filtered[i, ] <- (filtered[i, ] - min(filtered[i, ])) / (max(filtered[i, ]) - min(filtered[i, ])) * 100
  }
  
  
  # Predict ML:
  input <- array(NA, c(dim(filtered), number_of_derivs + 1))
  for (i in 1:nrow(filtered)) {
    input[i, , ] <- add_derivs(signal = filtered[i, ], number_of_derivs = number_of_derivs)
  }
  
  
  model <- load_model_tf(paste0('../models/', model_log$name[rows], '.h5'))
  uih_predictions <- model %>% predict(input)
  predictions_ML <- array(0, c(nrow(uih_predictions), ncol(uih_predictions)))
  for (i in 1:nrow(uih_predictions)) {
    predictions_ML[i, ] <- max.col(uih_predictions[i, , ])
  }
  #convert from dimension value 1,2,3,4 to 0,1,2,3
  predictions_ML <- predictions_ML - 1
  
  # Fill gaps
  predictions_ML_filled <- predictions_ML
  for (i in 1:nrow(predictions_ML)) {
    predictions_ML_filled[i, ] <- fill_wave_gaps(predictions_ML[i, ], 20)
  }
  
  # Prep ecgpuwave annotations
  predictions_ecgpuwave <- array(NA,dim(predictions_ML))
  for (i in 1:nrow(predictions_ecgpuwave)) {
    predictions_ecgpuwave[i,] <- ann_wfdb2continuous2(uih_ecgs_supp[[i]]$annotation)
  }
# }
  
  
# Analyze: relative to Rpeak -----------------------------------------------------------------
# for (lead in leads) {
  ML_Rpeak_analysis <- do.call(rbind, lapply(1:nrow(input[, , 1]), function(idx) {
    round(check_ann_prog_RPeaks(input[idx, , 1], predictions_ML[idx, ]), 2)
  }))
  
  ML_Rpeak_filled_analysis <- do.call(rbind, lapply(1:nrow(input[, , 1]), function(idx) {
    round(check_ann_prog_RPeaks(input[idx, , 1], predictions_ML_filled[idx, ]), 2)
  }))
  
  puwave_Rpeak_analysis <- do.call(rbind, lapply(1:nrow(input[, , 1]), function(idx) {
    round(check_ann_prog_RPeaks(input[idx, , 1], predictions_ecgpuwave[idx, ]), 2)
  }))
  
  ML_Rpeak_analysis_summary <- as.data.frame(t(colSums(ML_Rpeak_analysis))) / nrow(ML_Rpeak_analysis)
  ML_Rpeak_filled_analysis_summary <- as.data.frame(t(colSums(ML_Rpeak_filled_analysis))) / nrow(ML_Rpeak_filled_analysis)
  puwave_Rpeak_analysis_summary <- as.data.frame(t(colSums(puwave_Rpeak_analysis))) / nrow(puwave_Rpeak_analysis)
  
  output[[lead]] <- list(
    raw    = ML_Rpeak_analysis_summary,
    filled = ML_Rpeak_filled_analysis_summary
  )
  
  std_dev <- sapply(ML_Rpeak_filled_analysis, sd)
  uih_supp_Rpeak_std[[lead]] <- list(std_dev)
  
  if (lead == 'I') {
    output[['puwave']] <- list(
      raw = puwave_Rpeak_analysis_summary
    )
    uih_supp_Rpeak_std[['puwave']] <- sapply(puwave_Rpeak_analysis, sd)
  }
  print(paste('done with lead',lead))
  
}

uih_supp_Rpeak_analysis <- output


# Analyze: isoelectric ----------------------------------------------------
for (lead in leads) {
  # First, create ECG signal that is not normalized (so values can be in units of mV)
  filtered <- array(0, dim(uih_samples))
  for (i in 1:nrow(uih_samples)) {
    filtered[i, ] <- ecg_filter(uih_samples[i, ])
  }
  
  distance_from_isoelec <- data.frame(p_ML = array(NA,nrow(filtered)),
                                      qrs_ML = array(NA,nrow(filtered)),
                                      t_ML = array(NA,nrow(filtered)),
                                      p_pu = array(NA,nrow(filtered)),
                                      qrs_pu = array(NA,nrow(filtered)),
                                      t_pu = array(NA,nrow(filtered))
                                      )
  for (idx in 1:nrow(filtered)) {
    # Find Isoelectric point
    isoelec <- isoelec_find(signal = filtered[idx,],annotations = predictions_ML_filled[idx, ])
    
    # Calculate ML P-waves
    pwaves_ML <- make_wave_table(predictions_ML_filled[idx, ], wave_value = 1)
    pwaves_ML$onset_dist = NA
    pwaves_ML$offset_dist = NA
    for (i in 1:nrow(pwaves_ML)) {
      onset <- pwaves_ML$wave_on[i]
      offset <- pwaves_ML$wave_off[i]
      pwaves_ML$onset_dist[i] <- abs(filtered[idx,onset] - isoelec) 
      pwaves_ML$offset_dist[i] <- abs(filtered[idx,offset] - isoelec)
    }
    # Calculate ML QRS waves_ML
    qrswaves_ML <- make_wave_table(predictions_ML_filled[idx, ], wave_value = 2)
    qrswaves_ML$onset_dist = NA
    qrswaves_ML$offset_dist = NA
    for (i in 1:nrow(qrswaves_ML)) {
      onset <- qrswaves_ML$wave_on[i]
      offset <- qrswaves_ML$wave_off[i]
      qrswaves_ML$onset_dist[i] <- abs(filtered[idx,onset] - isoelec) 
      qrswaves_ML$offset_dist[i] <- abs(filtered[idx,offset] - isoelec)
    }
    # Calculate ML T-waves_ML
    twaves_ML <- make_wave_table(predictions_ML_filled[idx, ], wave_value = 3)
    twaves_ML$onset_dist = NA
    twaves_ML$offset_dist = NA
    for (i in 1:nrow(twaves_ML)) {
      onset <- twaves_ML$wave_on[i]
      offset <- twaves_ML$wave_off[i]
      twaves_ML$onset_dist[i] <- abs(filtered[idx,onset] - isoelec) 
      twaves_ML$offset_dist[i] <- abs(filtered[idx,offset] - isoelec)
    }
    
    
    
    # Calculate puwave P-waves
    pwaves_pu <- make_wave_table(predictions_ecgpuwave[idx, ], wave_value = 1)
    pwaves_pu$onset_dist = NA
    pwaves_pu$offset_dist = NA
    for (i in 1:nrow(pwaves_pu)) {
      onset <- pwaves_pu$wave_on[i]
      offset <- pwaves_pu$wave_off[i]
      pwaves_pu$onset_dist[i] <- abs(filtered[idx,onset] - isoelec) 
      pwaves_pu$offset_dist[i] <- abs(filtered[idx,offset] - isoelec)
    }
    # Calculate puwave QRS waves_pu
    qrswaves_pu <- make_wave_table(predictions_ecgpuwave[idx, ], wave_value = 2)
    qrswaves_pu$onset_dist = NA
    qrswaves_pu$offset_dist = NA
    for (i in 1:nrow(qrswaves_pu)) {
      onset <- qrswaves_pu$wave_on[i]
      offset <- qrswaves_pu$wave_off[i]
      qrswaves_pu$onset_dist[i] <- abs(filtered[idx,onset] - isoelec) 
      qrswaves_pu$offset_dist[i] <- abs(filtered[idx,offset] - isoelec)
    }
    # Calculate puwave T-waves_pu
    twaves_pu <- make_wave_table(predictions_ecgpuwave[idx, ], wave_value = 3)
    twaves_pu$onset_dist = NA
    twaves_pu$offset_dist = NA
    for (i in 1:nrow(twaves_pu)) {
      onset <- twaves_pu$wave_on[i]
      offset <- twaves_pu$wave_off[i]
      twaves_pu$onset_dist[i] <- abs(filtered[idx,onset] - isoelec) 
      twaves_pu$offset_dist[i] <- abs(filtered[idx,offset] - isoelec)
    }
    
    # ****pwaves_ML: dataframe of the pwaves of an individual ECG. Distance from isoelec is meant for all ECGs combined 
    distance_from_isoelec$p_ML[idx] <- mean(pwaves_ML$onset_dist) + mean(pwaves_ML$offset_dist)
    distance_from_isoelec$qrs_ML[idx] <- mean(qrswaves_ML$onset_dist) + mean(qrswaves_ML$offset_dist)
    distance_from_isoelec$t_ML[idx] <- mean(twaves_ML$onset_dist) + mean(twaves_ML$offset_dist)
    distance_from_isoelec$p_pu[idx] <- mean(pwaves_pu$onset_dist) + mean(pwaves_pu$offset_dist)
    distance_from_isoelec$qrs_pu[idx] <- mean(qrswaves_pu$onset_dist) + mean(qrswaves_pu$offset_dist)
    distance_from_isoelec$t_pu[idx] <- mean(twaves_pu$onset_dist) + mean(twaves_pu$offset_dist)
    
 
  }
  colMeans(distance_from_isoelec, na.rm=TRUE)
}
  
# print(ML_Rpeak_analysis_summary)
# print(ML_Rpeak_filled_analysis_summary)
# print(puwave_Rpeak_analysis_summary)




# Graph relative Rpeak (grouped by lead) -------------------------------------------
library(plotly)
library(dplyr)
library(tidyr)

load('../uih_supp_Rpeak_analysis.RData')
# Average across all leads
uih_supp_Rpeak_analysis[['mean']]$filled <- round(colMeans(do.call(rbind,lapply(1:12, function(idx) {
  mat <- uih_supp_Rpeak_analysis[[idx]]$filled}))),4)

# Desired order
lead_order <- c("puwave","mean","I","II","III","AVR","AVL","AVF",
                "V1","V2","V3","V4","V5","V6")

# 1. Build a tidy data frame directly
df <- bind_rows(
  lapply(lead_order, function(lead) {
    tmp <- uih_supp_Rpeak_analysis[[lead]]$filled
    tmp$lead <- lead
    tmp
  }),
  .id = NULL
)

# 2. Reshape to long format
df_long <- df %>%
  pivot_longer(
    cols = c(
             missed_QRS, missed_P, missed_T,
             duplicate_QRS, duplicate_P, duplicate_T
    ),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    # Force the x-axis order
    lead = factor(lead, levels = lead_order)
  )

# 3. Plot grouped bar chart
fig <- plot_ly(
  df_long,
  x = ~lead,
  y = ~value,
  color = ~metric,
  type = "bar"
) %>%
  layout(
    barmode = "group",
    xaxis = list(title = "ECG Lead"),
    yaxis = list(title = "Value"),
    legend = list(title = list(text = "Metric"))
  )

fig






# graph relative Rpeak (grouped by wave) ----------------------------------
library(dplyr)
lead <- 'T'

frame_title <- paste0('missed_',lead)
graph_title <- paste0('Missed ',lead,'-wave rate')
if (lead == 'QRS') {
  graph_title <- paste0('Missed ',lead,' Complex rate')
}
load('../uih_supp_Rpeak_analysis.RData')
# Average across all leads
uih_supp_Rpeak_analysis[['mean']]$filled <- round(colMeans(do.call(rbind,lapply(1:12, function(idx) {
  mat <- uih_supp_Rpeak_analysis[[idx]]$filled}))),4)

uih_supp_Rpeak_analysis$puwave$filled <- uih_supp_Rpeak_analysis$puwave$raw

lead_order <- c("puwave","mean","I","II","III","AVR","AVL","AVF",
                "V1","V2","V3","V4","V5","V6")

# Desired order of leads
# 1. Build a data frame with just missed_P
df <- bind_rows(
  lapply(lead_order, function(lead) {
    tmp <- uih_supp_Rpeak_analysis[[lead]]$filled
    data.frame(
      lead = lead,
      graph_col = tmp[[frame_title]]
    )
  }),
  .id = NULL
)
error_bars <- bind_rows(
  lapply(lead_order[3:14], function(lead) {
    tmp <- uih_supp_Rpeak_std[[lead]][1]
    data.frame(
      lead = lead,
      graph_col = tmp[[1]][[frame_title]]
    )
  }),
  .id = NULL
)
error_bars <- rbind(error_bars,c('mean',mean(error_bars$graph_col)))
error_bars <- rbind(error_bars,c('puwave',uih_supp_Rpeak_std$puwave[frame_title]))

# 2. Ensure correct order
df$lead <- factor(df$lead, levels = lead_order)
error_bars$lead <- factor(error_bars$lead, levels = lead_order)
error_bars <- error_bars %>%
  arrange(factor(lead, levels = lead_order))

library(RColorBrewer)
fire_engine <- "#D50032"
navy_pier <- "#001E62"
chicago_blue <- "#41B6E6"

# Build a color vector based on lead
colors <- ifelse(df$lead == "puwave", fire_engine,
                 ifelse(df$lead == "mean", navy_pier, chicago_blue))

# Plot
fig <- plot_ly(
  df,
  x = ~lead,
  y = ~graph_col,
  type = "bar",
  name = frame_title,
  marker = list(color = colors)   # assign per-bar colors
  
) %>%
  layout(
    title = list(
      text = graph_title,          # use your variable as the plot title
      x = 0.5,                     # center the title
      xanchor = "center",
      yanchor = "top",
      font = list(size = 56)       # adjust title font size
    ),
    xaxis = list(
      title = "ECG Lead",
      titlefont = list(size = 42), # increase x-axis label size
      tickfont  = list(size = 24)  # increase tick label size
    ),
    yaxis = list(
      title = graph_title,
      titlefont = list(size = 42), # increase y-axis label size
      tickfont  = list(size = 24)  # increase tick label size
    ),
    barmode = "group"
  )

fig


# plot 2 ------------------------------------------------------------------

library(dplyr)
library(plotly)

# Ensure both dfs are aligned and ordered
df <- df %>%
  mutate(lead = factor(lead, levels = lead_order)) %>%
  arrange(lead)

error_bars <- error_bars %>%
  mutate(lead = factor(lead, levels = lead_order)) %>%
  arrange(lead)

# Merge into one dataframe
plot_df <- df %>%
  left_join(error_bars, by = "lead", suffix = c("", "_err"))

# Assign colors per lead
plot_df <- plot_df %>%
  mutate(color = ifelse(lead == "puwave", fire_engine,
                        ifelse(lead == "mean", navy_pier, chicago_blue)))

# Plot with per-bar error bars
fig <- plot_ly(
  plot_df,
  x = ~lead,
  y = ~graph_col,
  type = "bar",
  name = frame_title,
  marker = list(color = ~color),
  error_y = ~list(
    type = "data",
    array = graph_col_err,   # error values from merged df
    color = color,           # match bar color
    thickness = 2,
    width = 5
  )
) %>%
  layout(
    title = list(
      text = graph_title,
      x = 0.5,
      xanchor = "center",
      yanchor = "top",
      font = list(size = 56)
    ),
    xaxis = list(
      title = "ECG Lead",
      titlefont = list(size = 42),
      tickfont  = list(size = 24)
    ),
    yaxis = list(
      title = graph_title,
      titlefont = list(size = 42),
      tickfont  = list(size = 24),
      range = c(0, NA)
    ),
    barmode = "group"
  )

fig


# Analyze: distance to isoelectric line -----------------------------------
# Distance from P/QRS/T wave offset and onset from the P-T interval (or median val of entire ECG)

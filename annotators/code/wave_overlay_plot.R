# prep --------------------------------------------------------------------

load('../uih_ecgs_supp.RData')
ecgs <- uih_ecgs_supp[1:10]
ann <- predict_ecgs(ecgs)
# "I"   "II"  "III" "AVR" "AVL" "AVF" "V1"  "V2"  "V3"  "V4"  "V5"  "V6" 
# ann_wfdb2continuous2(ann[1]$annotation[['I']])


# Combine -----------------------------------------------------------------
idx <- 4
lead <- 5
wave <- 1

leads <- c("I", "II", "III", "AVR", "AVL", "AVF", 
           "V1", "V2", "V3", "V4", "V5", "V6")

# Apply your function to each lead, producing a list of logical vectors
lead_vectors <- lapply(leads, function(ld) {
  ann_wfdb2continuous2(ann[[idx]]$annotation[[ld]]) == wave
})

# Convert to a matrix (12 x 5000) and take column means
avg_vector <- colMeans(do.call(rbind, lead_vectors))

# avg_vector is now a single vector of length 5000 with values between 0 and 1

# Signal for plots
sig <- ecg_filter(ann[[idx]]$signal[[leads[lead]]])

# Plot avg --------------------------------------------------------------------
library(plotly)

x_vals <- 1:5000
y_vals <- sig
color <- round(avg_vector, 3)

# Build hover text vector (length 5000)
hover_text <- paste0(
  "Time: ", round(x_vals, 3),
  "<br>Signal: ", round(y_vals, 3),
  "<br>Avg: ", color
)

fig <- plot_ly(
  x = x_vals,
  y = y_vals,
  type = 'scatter',
  mode = 'markers',
  marker = list(
    color = color,
    colorscale = 'Viridis',
    showscale = FALSE,   # hide colorbar
    size = 6
  ),
  hoverinfo = "skip"    # only show custom text
)

# Add light gray line trace
fig <- fig %>%
  add_trace(
    x = x_vals,
    y = y_vals,
    type = 'scatter',
    mode = 'lines',
    line = list(color = 'lightgray', width = 1),
    name = "Signal shape",
    text = hover_text,     # attach hover text
    hoverinfo = "text",   # skip hover for line
    name = "Points"
  ) %>%
  layout(
    title = paste("Lead", leads[lead], "with avg_vector coloring"),
    xaxis = list(title = "Time (s)"),
    yaxis = list(title = "Signal amplitude")
  )

fig


# entropy -----------------------------------------------------------------

# Suppose preds is an N x 12 matrix of integers (0–3)
# rows = time indices, cols = leads/models

wfdb <- ann[[idx]]
# array <- array(NA,c(5000,length(wfdb$annotation)))


calc_ecg_entropy <- function(wfdb) {
  # Input: wfdb object, with 12 lead annotations
  if (any(class(wfdb) == 'egm')) {
    ann <- wfdb$annotation
  }
  array <- array(NA,c(5000,length(ann)))
  
  # Convert from wfdb to array format
  for (i in 1:length(wfdb$annotation)) {
    array[,i] <- ann_wfdb2continuous2(wfdb$annotation[[i]])
  }
  entropy_vec <- apply(array, 1, calc_entropy)
  
}

calc_entropy <- function(row) {
  # Input: single annotation time point, all 12 leads
    # tabulate counts for classes 0:3
    counts <- tabulate(factor(row, levels = 0:3))
    probs <- counts / sum(counts)
    # Shannon entropy
    H <- -sum(probs[probs > 0] * log(probs[probs > 0]))
    # normalize to [0,1]
    H / log(length(probs))
}


# Edit --------------------------------------------------------------------
# if less than 5 indices, remove

# Can you detect_QRS to find P, QRS and T
#   for more rigorous P/T wave detection, if multiple waves found within window (ie if pt is tachy), take closest to QRSh

# if missing, 
    # (especially if it's a whole wave, ie no consecutive overlap)
    # Note as low confidence
    # add P-wave using other leads

# if a lead is missing P-waves, note that as low confidence (or QRS, T)


# confusion -----------------------------------------------------------------
library(dplyr)
load('../models/model_log.RData')

lead <- 6
which(model_log$lead == lead)

uih_conf_indicies <- which(model_log$lead == lead & !is.na(model_log$confusion_uih))
lapply(uih_conf_indicies, function(idx) {model_log$confusion_uih[[idx]]$table})

lapply(uih_conf_indicies, function(idx) {model_log$confusion_uih[[idx]]$matrix[,7]})

conf_indicies <- which(model_log$lead == lead & !is.na(filt$confusion))
lapply(conf_indicies, function(idx) {filt$confusion[[idx]]$table})


# temp --------------------------------------------------------------------
source('annotator_prep_functions.R')
load('../models/model_log.RData')
ecgs <- generate_ecgs(20,'sinus')
# lead <- 5
ecgs_top <- predict_ecgs(ecgs,lead=lead)
ecgs_new <- predict_ecgs(ecgs,lead=lead,model_number = 916)
leads <- c('I','II','III','AVR','AVL','AVF','V1','V2','V3','V4','V5','V6')
idx=0

idx <- idx+1
sig <- c(ecg_filter(ecgs_top[[idx]]$signal[[leads[lead]]]))
plot_func2(sig,ann_wfdb2continuous2(ecgs_top[[idx]]$annotation[[leads[lead]]]),linewidth = 0.5,pointsize = 0.5)
Sys.sleep(1)
plot_func2(sig,ann_wfdb2continuous2(ecgs_new[[idx]]$annotation[[leads[lead]]]),linewidth = 0.5,pointsize = 0.5)
# Sys.sleep(1)
# plot_func2(sig,ann_wfdb2continuous2(ecgs[[idx]]$annotation),linewidth = 0.5,pointsize = 0.5)
# plot_func2(sig,ann_wfdb2continuous2(ecgs_new2[[idx]]$annotation[[leads[lead]]]),linewidth = 0.5,pointsize = 0.5)


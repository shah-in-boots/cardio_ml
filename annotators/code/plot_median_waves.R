# functions ---------------------------------------------------------------
ann_wfdb2compact <- function(ann_wfdb) {
  library(dplyr)
  class(ann_wfdb) <- "data.frame"
  ann_wfdb <- ann_wfdb %>% mutate(idx = row_number())
  
  # Extract waves
  ann_compact <- ann_wfdb %>%
    mutate(
      type_lead = lead(type),
      type_lag = lag(type),
      sample_lead = lead(sample),
      sample_lag = lag(sample)
    ) %>%
    filter(type %in% c("p", "N", "t"), type_lag == "(", type_lead == ")") %>%
    transmute(
      type,
      onset = sample_lag,
      peak = sample,
      offset = sample_lead
    )
  return(ann_compact)
}

ann_compact2continuous <- function(ann_compact) {
  # Assuming your dataframe is called wave_table
  ann_continuous <- rep(0, 5000)
  
  for (i in seq_len(nrow(ann_compact))) {
    start <- ann_compact$onset[i]
    end <- ann_compact$offset[i]
    
    value <- switch(ann_compact$type[i],
                    "p" = 1,
                    "N" = 2,
                    "t" = 3)
    
    ann_continuous[start:end] <- value
  }
  return(ann_continuous)
}

# code --------------------------------------------------------------------
# Parameters you can tune
window_left  <- 100   # extra points on the left of the peak
window_right <- 100   # extra points on the right of the peak
min_wave_len    <- 10    # minimum contiguous length to accept a T-wave segment
prop_thresh  <- 0.5   # threshold to mark average onset/offset from proportions
wave <- 'T'
color_limits <- 'rel' # abs or 'rel'


signal <- filtered
ML_ann <- predictions_ML
puwave_ann <- predictions_ecgpuwave


waves <- c('P','QRS','T')

wave_val <- which(waves == wave)

# For both pu wave and ML: 
  # Define a t-wave if there is a puwave t-wave OR ML t-wave at a given index

# Create list/dataframe of all T-waves - signal, ML annotation and puwave annotation
#   include signal/ann ~50-100 indices on each side (ie window = 50 or 100)

# For each T-wave, define alignment point 
  # Best method is likely to find the global max
    # Note: exclude the extra "window" indicies, as these could include QRS or P-waves, which have their own maxima

# Will need to transform each T-wave snippet to be same length 
  # pad with same value at each end 


# The signal variable is mV values, vector
# The ML_ann and puwave_ann are vectors with integer values. 0 for no wave, 1 for p-wave, 2 QRS, 3 T-wave





# chatgpt -----------------------------------------------------------------
library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)
library(scales)

# 1) Find contiguous T-wave segments (start, end) in a 0/1 indicator vector
find_segments <- function(vec_bin) {
  if (all(vec_bin == 0) || length(vec_bin) == 0) return(matrix(numeric(0), ncol = 2))
  r <- rle(vec_bin)
  ends <- cumsum(r$lengths)
  starts <- c(1, head(ends, -1) + 1)
  segs <- cbind(starts[r$values == 1], ends[r$values == 1])
  # filter by minimum wave length
  segs <- segs[(segs[,2] - segs[,1] + 1) >= min_wave_len, , drop = FALSE]
  segs
}

# 2) Extract the alignment point (peak) within the true T segment (no padding)
find_wave_peak <- function(sig_row, seg_start, seg_end) {
  seg <- sig_row[seg_start:seg_end]
  which.max(seg) + seg_start - 1
}

# 3) Build a fixed-length window around alignment point; pad with NA to keep equal length
# Returns indices used for the window (same across signal and annotations)
build_window_indices <- function(peak_idx, n_cols, left = window_left, right = window_right) {
  idx <- (peak_idx - left):(peak_idx + right)
  idx[idx < 1] <- NA
  idx[idx > n_cols] <- NA
  idx
}

# 4) Extract a windowed vector from a row with NA padding where indices are NA
extract_window <- function(row_vec, idx) {
  out <- rep(NA_real_, length(idx))
  valid <- !is.na(idx)
  out[valid] <- row_vec[idx[valid]]
  out
}

# Inputs: signal, ML_ann, puwave_ann (matrices: n_samples x 5000)
stopifnot(is.matrix(signal), is.matrix(ML_ann), is.matrix(puwave_ann))
stopifnot(all(dim(signal) == dim(ML_ann)), all(dim(signal) == dim(puwave_ann)))
n_samples <- nrow(signal)
n_cols    <- ncol(signal)

# 1) Define combined wave indicator: T if ML OR rules-based says wave_val
combined_wave <- (ML_ann == wave_val) | (puwave_ann == wave_val)

# 2) For each sample, find segments and choose one wave (e.g., first, or largest)
# Here: choose the longest contiguous wave segment per sample, if any
get_longest_segment <- function(bin_row) {
  segs <- find_segments(vec_bin = bin_row)
  if (nrow(segs) == 0) return(NULL)
  lens <- segs[,2] - segs[,1] + 1
  segs[which.max(lens), ]
}

# 3) Collect aligned windows for signal, ML wave indicator, and rules-based wave indicator
#    Alignment = peak of signal within the wave segment only
sig_windows_list <- list()
ml_windows_list  <- list()
pu_windows_list  <- list()


# chatgpt run -------------------------------------------------------------

for (i in seq_len(n_samples)) {
  seg <- get_longest_segment(bin_row = combined_wave[i, ])
  if (is.null(seg)) next
  
  peak_idx <- find_wave_peak(signal[i, ], seg_start = seg[1], seg_end = seg[2])
  idx_win  <- build_window_indices(peak_idx, n_cols)
  
  # Extract windowed snippets (signal) and binary wave indicators (1 for T, else 0)
  sig_windows_list[[length(sig_windows_list) + 1]] <- extract_window(signal[i, ], idx_win)
  
  ml_wave_bin <- as.numeric(ML_ann[i, ] == wave_val)
  pu_wave_bin <- as.numeric(puwave_ann[i, ] == wave_val)
  
  ml_windows_list[[length(ml_windows_list) + 1]] <- extract_window(ml_wave_bin, idx_win)
  pu_windows_list[[length(pu_windows_list) + 1]] <- extract_window(pu_wave_bin, idx_win)
}

# 4) Stack windows into matrices (rows = waves, cols = timepoints relative to peak)
to_matrix <- function(lst) {
  if (length(lst) == 0) return(matrix(NA_real_, nrow = 0, ncol = window_left + window_right + 1))
  do.call(rbind, lst)
}
sig_mat <- to_matrix(sig_windows_list)
ml_mat  <- to_matrix(ml_windows_list)
pu_mat  <- to_matrix(pu_windows_list)

# 5) Compute averages: signal (mean, ignoring NA), annotations (proportion of 1’s ignoring NA)
avg_signal <- colMeans(sig_mat, na.rm = TRUE)

prop_ml <- colMeans(ml_mat, na.rm = TRUE)  # proportion of samples calling T at each position
prop_pu <- colMeans(pu_mat, na.rm = TRUE)

# 6) Derive average onset/offset from proportions using a threshold
find_onset_offset <- function(prop_vec, thresh = prop_thresh) {
  above <- which(prop_vec >= thresh)
  if (length(above) == 0) return(c(NA_integer_, NA_integer_))
  c(min(above), max(above))
}
ml_on_off <- find_onset_offset(prop_ml)
pu_on_off <- find_onset_offset(prop_pu)

# 7) Build plotting frame
time_rel <- (-window_left):window_right  # time index relative to peak

df_plot <- tibble(
  Time = time_rel/5000,
  Signal = avg_signal,
  ML = prop_ml,
  puwave = prop_pu
)


# custom plot -------------------------------------------------------------

library(plotly)
library(ggplot2)
library(htmlwidgets)
library(dplyr)

df_plot <- tibble(
  Time = time_rel/5000,
  Signal = avg_signal,
  ML = prop_ml,
  puwave = prop_pu
)

if (color_limits == 'abs') {
  color_limits_ml <- c(0,1)
  color_limits_pu <- c(0,1)
  } else if (color_limits == 'rel') {
    
    df_plot <- df_plot %>%
      mutate(
        ML     = round(ML / max(ML, na.rm = TRUE), 2),
        puwave = round(puwave / max(puwave, na.rm = TRUE), 2)
      )
    
    
    color_limits_ml <- c(0,max(df_plot$ML))
    color_limits_pu <- c(0,max(df_plot$puwave))
  }

# p_ml <- ggplotly(ggplot(df_plot, aes(Time, Signal, color = ML)) + 
#   geom_point(size = 1.5) +
#   scale_x_continuous(breaks = seq(0, 10, 1)) + 
#   scale_color_gradient(limits = color_limits_ml, low = "yellow", high = "blue") +
#   theme(legend.position = "none",
#         legend.title = element_blank()) +
#   theme_bw() + geom_line(aes(y = Signal), color = "gray", linewidth = 0.4))
# 
# 
# p_pu <- ggplotly(ggplot(df_plot, aes(Time, Signal, color = puwave)) + 
#   geom_point(size = 1.5) +
#   scale_x_continuous(breaks = seq(0, 10, 1)) + 
#   scale_color_gradient(limits = color_limits_pu, low = "yellow", high = "blue") +
#   theme(legend.position = "none",
#         legend.title = element_blank()) +
#   theme_bw() + geom_line(aes(y = Signal), color = "gray", linewidth = 0.4))
# 
# 
# temp_file1 <- tempfile(fileext = ".html")
# saveWidget(p_ml, temp_file1, selfcontained = TRUE)
# browseURL(temp_file1)
# 
# temp_file2 <- tempfile(fileext = ".html")
# saveWidget(p_pu, temp_file2, selfcontained = TRUE)
# browseURL(temp_file2)


# plot combined, same plot ggplot ------------------------------------------------
# library(ggnewscale)
# 
# 
# p_comb <- ggplot(df_plot, aes(Time)) +
#   # ML layer (offset +5)
#   geom_point(aes(y = Signal + 10, color = ML), size = 1.5) +
#   geom_line(aes(y = Signal + 10), color = "gray", linewidth = 0.4) +
#   scale_color_gradient(limits = color_limits_ml, low = "yellow", high = "blue", name = "ML") +
#   new_scale_color() +
#   # puwave layer
#   geom_point(aes(y = Signal, color = puwave), size = 1.5) +
#   geom_line(aes(y = Signal), color = "gray", linewidth = 0.4) +
#   scale_color_gradient(limits = color_limits_pu, low = "yellow", high = "blue", name = "puwave") +
#   scale_x_continuous(breaks = seq(0, 10, 1)) +
#   theme_bw()
# 
# p_comb


# plot combined, same plot plotly -----------------------------------------
library(plotly)

# Example assumes df_plot has columns: Time, Signal, ML, puwave
if (color_limits == 'rel') {
  title_combined_plot <- paste0("ML and Rules-Based Annotation on Median ",wave,'-wave')
} else if (color_limits == 'abs') {
  title_combined_plot <- paste0("ML and Rules-Based Annotation on Median ",wave,'-wave')
}

point_size <- 12 # 6 normal, 12 for poster

p_comb <- plot_ly() %>%
  # --- ML layer (offset +10) ---
  add_trace(
    data = df_plot,
    x = ~Time,
    y = ~Signal + 10,
    type = "scatter",
    mode = "markers",
    marker = list(
      size = point_size,
      color = ~ML,
      colorscale = list(c(0, "yellow"), c(color_limits_ml[2], "blue")),
      cmin = 0,
      cmax = 1,
      colorbar = list(
        title = "ML",
        y = 0.8,
        len = 0.4,        # shrink vertical length (40% of plot height)
        thickness = 12,   # shrink width (pixels)
        titlefont = list(size = 10),
        tickfont  = list(size = 9)
      )
    ),
    name = "ML"
  ) %>%
  add_trace(
    data = df_plot,
    x = ~Time,
    y = ~Signal + 10,
    type = "scatter",
    mode = "lines",
    line = list(color = "gray", width = 1),
    showlegend = FALSE,
    inherit = FALSE
  ) %>%
  # --- puwave layer (baseline) ---
  add_trace(
    data = df_plot,
    x = ~Time,
    y = ~Signal,
    type = "scatter",
    mode = "markers",
    marker = list(
      size = point_size,
      color = ~puwave,
      colorscale = list(c(0, "yellow"), c(1, "blue")),
      cmin = 0,
      cmax = 1,
      colorbar = list(
        title = "puwave",
        y = 0.2,
        len = 0.4,        # shrink vertical length
        thickness = 12,   # shrink width
        titlefont = list(size = 10),
        tickfont  = list(size = 9)
      )
    ),
    name = "puwave"
  ) %>%
  add_trace(
    data = df_plot,
    x = ~Time,
    y = ~Signal,
    type = "scatter",
    mode = "lines",
    line = list(color = "gray", width = 1),
    showlegend = FALSE,
    inherit = FALSE
  ) %>%
  layout(
    xaxis = list(title = "Time", tick0 = 0, dtick = 1),
    yaxis = list(title = "Signal"),
    title = list(
      text = title_combined_plot,
      titlefont = list(size = 60),
      tickfont  = list(size = 24)
    ),
    margin = list(t = 100),
    font = list(size = 24)
  )

p_comb



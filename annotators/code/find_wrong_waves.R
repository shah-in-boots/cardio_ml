
# df output ---------------------------------------------------------------

validate_ecg_waves <- function(signal, ann_wfdb, 
                                            fs = 500,
                                            terminal_boundary = 250,
                                            pr_window = c(80, 350),   # ms before R
                                            rt_window = c(120, 500)) { # ms after R
  library(dplyr)
  library(EGM)
  
  # Convert ms to samples
  pr_window <- round(pr_window * fs / 1000)
  rt_window <- round(rt_window * fs / 1000)
  
  # Collapse annotations
  ann_compact <- ann_wfdb2compact(ann_wfdb)
  
  # Detect R peaks using detect_QRS, remove those at termini
  rpeaks <- EGM::detect_QRS(signal, frequency = fs)
  rpeaks <- rpeaks[rpeaks < (length(signal) - terminal_boundary)]
  rpeaks <- rpeaks[rpeaks > terminal_boundary]
  
  # --- QRS CHECK SECTION ---
  # Identify QRS from annotations, remove those at termini
  qrs_ann <- ann_compact %>% filter(type == "N")
  qrs_ann <- qrs_ann %>% filter(peak < (length(signal) - terminal_boundary) & peak > terminal_boundary)
  
  # For each rpeak, check if it falls inside any QRS annotation
  rpeak_matches <- sapply(rpeaks, function(r) {
    any(r >= qrs_ann$onset & r <= qrs_ann$offset)
  })
  
  # For each QRS annotation, check if it has a matching rpeak
  qrs_matches <- sapply(1:nrow(qrs_ann), function(i) {
    any(rpeaks >= qrs_ann$onset[i] & rpeaks <= qrs_ann$offset[i])
  })
  
  # Missed = rpeak without QRS
  missed_rpeaks <- rpeaks[!rpeak_matches]
  
  # Duplicate = QRS without rpeak
  duplicate_qrs <- qrs_ann[!qrs_matches, ]
  
  # Collapse duplicate QRS indices into one string
  duplicate_idx <- if (nrow(duplicate_qrs) > 0) {
    paste(round(rowMeans(duplicate_qrs[, c("onset", "offset")])), collapse = ";")
  } else NA
  
  qrs_flags <- data.frame(
    rpeak = rpeaks,
    qrs_flag = ifelse(rpeaks %in% missed_rpeaks, "missed", NA),
    qrs_flag_idx = ifelse(rpeaks %in% missed_rpeaks, as.character(rpeaks), NA),
    stringsAsFactors = FALSE
  )
  
  # Add one row for all duplicate QRS (collapsed)
  if (!is.na(duplicate_idx)) {
    dup_row <- data.frame(
      rpeak = NA,
      qrs_flag = "duplicate",
      qrs_flag_idx = duplicate_idx,
      stringsAsFactors = FALSE
    )
    qrs_flags <- bind_rows(qrs_flags, dup_row)
  }
  # --- END QRS CHECK SECTION ---
  
  # Initialize results as a list of rows
  results <- lapply(rpeaks, function(r) {
    # Define search windows
    pr_range <- (r - pr_window[2]):(r - pr_window[1])
    rt_range <- (r + rt_window[1]):(r + rt_window[2])
    rt_range <- rt_range[rt_range <= 5000]
    
    # Restrict to annotations within windows
    p_candidates <- ann_compact %>% filter(type == "p", onset %in% pr_range)
    t_candidates <- ann_compact %>% filter(type == "t", offset %in% rt_range)
    
    # Flags
    p_flag <- NA
    p_flag_idx <- NA
    if (nrow(p_candidates) == 0) {
      p_flag <- "missing"
      p_flag_idx <- round(mean(pr_range))
    } else if (nrow(p_candidates) > 1) {
      p_flag <- "duplicate"
      p_flag_idx <- paste(p_candidates$peak[-1], collapse = ";")
    }
    
    t_flag <- NA
    t_flag_idx <- NA
    if (nrow(t_candidates) == 0) {
      t_flag <- "missing"
      t_flag_idx <- round(mean(rt_range))
    } else if (nrow(t_candidates) > 1) {
      t_flag <- "duplicate"
      t_flag_idx <- paste(t_candidates$peak[-1], collapse = ";")
    }
    
    return(data.frame(rpeak = r,
                      p_flag = p_flag,
                      p_flag_idx = p_flag_idx,
                      t_flag = t_flag,
                      t_flag_idx = t_flag_idx,
                      stringsAsFactors = FALSE))
  })
  
  results_df <- bind_rows(results)
  
  # Merge QRS flags into results
  # final <- results_df %>%
  #   left_join(qrs_flags, by = "rpeak")
  
  final <- results_df %>%
    left_join(qrs_flags %>% filter(!is.na(rpeak)), by = "rpeak")
  
  # Add the duplicate summary row back in
  dup_row <- qrs_flags %>% filter(is.na(rpeak))
  final <- bind_rows(final, dup_row)
  
  
  return(final)
}




# temp (remove) -----------------------------------------------------------
idx=idx+1
plot_func2(ecg_filter(ann[[idx]]$signal$AVL),ann[[idx]]$annotation$AVL)
flag_missing_or_duplicate_waves(signal = ecg_filter(ann[[idx]]$signal$AVL),ann_wfdb=ann[[idx]]$annotation$AVL)
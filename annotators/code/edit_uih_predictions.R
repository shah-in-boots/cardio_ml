# setup -------------------------------------------------------------------
# To download the files, run these in termal:
# scp -o MACs=hmac-sha2-256 -r dseaney2@131.193.182.245:/home/dseaney2/cardio_darbar_chi_link/common/cohorts/wes_ml/annotator_models/uih_ecgs.RData 'C:\Users\darre\OneDrive\Documents\UICOM Research\annotators'
# scp -o MACs=hmac-sha2-256 -r dseaney2@131.193.182.245:/home/dseaney2/cardio_darbar_chi_link/common/cohorts/wes_ml/annotator_models/uih_predictions 'C:\Users\darre\OneDrive\Documents\UICOM Research\annotators'
# scp -o MACs=hmac-sha2-256 -r dseaney2@131.193.182.245:/home/dseaney2/cardio_darbar_chi_link/common/cohorts/wes_ml/annotator_models/code/annotator_prep_functions.R 'C:\Users\darre\OneDrive\Documents\UICOM Research\annotators'

# setwd("C:/Users/darre/OneDrive/Documents/UICOM Research/annotators") # adjust working directory as needed

# Run this code in R once:
# load('uih_ecgs.RData')
# for (i in 1:length(uih_ecgs)) {
#   attributes(uih_ecgs[[i]]$header)$record_line$record_name <- 0 # single value
#   uih_ecgs[[i]]$header$file_name <- array(0,12) # vector of 12
# }
# deid_uih_ecgs <- uih_ecgs
# save(deid_uih_ecgs,file='deid_uih_ecgs.RData')
# rm(uih_ecgs)

# Then **DELETE** uih_ecgs.RData after deid_uih_ecgs.RData is saved to the folder you're using

# Create a folder for the adjusted annotations:
# if (!dir.exists("uih_predictions_adjusted")) {
#   dir.create("uih_predictions_adjusted")
# }



# functions ---------------------------------------------------------------
# 2 functions:
# 1. convert wfdb table to a more compact version, easier to edit
# 2. convert wfdb compact version to continuous (then use ann_continuous2wfdb function)
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


# Prepare file ---------------------------------------------------------------------
library(fs)
library(htmlwidgets)
setwd("C:/Users/darre/OneDrive/Documents/UICOM Research/annotators") # adjust as needed
source('code/annotator_prep_functions.R') 
load('deid_uih_ecgs.RData') # might need to change this path

sample <- 1+sample # next: 12
lead <- 'AVR'
# I	II	III	AVF	AVL	AVR	V1	V2	V3	V4	V5	V6

load(paste0('uih_predictions/matched_predictions_',lead,'.RData')) # might need to change path. Can load in other leads as needed
annotation_list <- get(paste0('matched_predictions_',lead)) # this allows you to dynamically set a variable. So when 
# lead = 'I', it's the equivalent of: annotation_list <- matched_predictions_I

signal <- deid_uih_ecgs[[sample]]$signal[[lead]]
annotation_wfdb <- annotation_list[[sample]]

# Plot the function:
original_plot <- plot_func(ecg_filter(signal),ann_wfdb2continuous2(annotation_wfdb),pointsize = 0.5)
# NOTE: try playing around with the point size varialbe (default is 1). Linewidth could help too
# NOTE: always use the ecg_filter() function
# Open plot in browser:
temp_file <- tempfile(fileext = ".html")
saveWidget(original_plot, temp_file, selfcontained = TRUE)
browseURL(temp_file)

# Print the annotation matrix to command line:
annotation_compact <- ann_wfdb2compact(annotation_wfdb)
annotation_compact

# Create a csv file: 
directory_name <- 'uih_predictions_adjusted'
file_name <- paste0(sprintf("%03d", c(sample)),'_',lead)
full_name <- fs::path(directory_name,file_name,ext='csv')
if (file.exists(full_name)) {
  stop(paste('Warning, file exists. Run script after this line to continue'))
  }
write.csv(annotation_compact, file = full_name, row.names = TRUE)


# Edit file ---------------------------------------------------------------

# Now open that csv file in excel

# As you make adjustments, manually save the csv (Control + S). Then reload via:
annotation_compact_revised <- read.csv(full_name)


# Clean up:
# remove NA rows
# annotation_compact_revised <- annotation_compact_revised[!is.na(annotation_compact_revised$onset),]
# renumber as needed
# annotation_compact_revised$X <- 1:nrow(annotation_compact_revised)
# write.csv(annotation_compact_revised,file = full_name)


# If there are other symbols (ie V for PVC), remove them for plotting purposes.
# annotation_compact_revised <- annotation_compact_revised[annotation_compact_revised$type %in% c('p', 'N', 't'), ]
annotation_compact_revised$X <- seq_len(nrow(annotation_compact_revised))


annotation_continuous_revised <- ann_compact2continuous(annotation_compact_revised)
annotation_wfdb_revised <- ann_continuous2wfdb(annotation_continuous_revised)

# And replot:
revised_plot <- plot_func(ecg_filter(signal),annotation_continuous_revised,pointsize= 0.5)
temp_file <- tempfile(fileext = ".html")
saveWidget(revised_plot, temp_file, selfcontained = TRUE)
browseURL(temp_file)


# annotation_compact_revised$offset[annotation_compact_revised$type == 'p'] = 
#   annotation_compact_revised$offset[annotation_compact_revised$type == 'p']+23
# checking ----------------------------------------------------------------
# sample <- 6
# file_name <- paste0(sprintf("%03d", c(sample)),'_',lead)
# 
# signal <- deid_uih_ecgs[[sample]]$signal[[lead]]
# annotation_wfdb <- annotation_list[[sample]]
# 
# annotation_compact_revised <- read.csv(paste0('../uih_predictions_adjusted/shreyank/',file_name,'.csv'))
# annotation_continuous_revised <- ann_compact2continuous(annotation_compact_revised)
# 
# 
# plot_func(ecg_filter(signal),ann_wfdb2continuous2(annotation_wfdb),pointsize = 0.5)
# plot_func(ecg_filter(signal),annotation_continuous_revised,pointsize= 0.5)
# 

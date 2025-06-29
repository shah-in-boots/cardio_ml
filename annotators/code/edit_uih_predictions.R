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


# run ---------------------------------------------------------------------
library(fs)
setwd("C:/Users/darre/OneDrive/Documents/UICOM Research/annotators/code") # adjust as needed
source('annotator_prep_functions.R') 
load('../deid_uih_ecgs.RData') # might need to change this path
load('../uih_predictions/matched_predictions_I.RData') # might need to change path. Can load in other leads as needed

sample <- 1
lead <- 'I'

annotation_list <- get(paste0('matched_predictions_',lead)) # this allows you to dynamically set a variable. So when 
# lead = 'I', it's the equivalent of: annotation_list <- matched_predictions_I

signal <- deid_uih_ecgs[[sample]]$signal[[lead]]
annotation_wfdb <- annotation_list[[sample]]

# Plot the function:
plot_func(ecg_filter(signal),ann_wfdb2continuous2(annotation_wfdb),pointsize = 0.5)
# NOTE: try playing around with the point size varialbe (default is 1). Linewidth could help too
# NOTE: always use the ecg_filter() function

# Print the annotation matrix to command line:
annotation_compact <- ann_wfdb2compact(annotation_wfdb)
annotation_compact

# Create a csv file: 
directory_name <- 'uih_predictions_adjusted'
file_name <- paste0(sprintf("%03d", c(sample)),'_',lead)
full_name <- fs::path(directory_name,file_name,ext='csv')
if (file.exists(full_name)) {
  print(paste('Warning, file exists. Run script after this line to continue'))
  break
  }
write.csv(annotation_compact, file = full_name, row.names = TRUE)

# Now open that csv file in excel

# As you make adjustments, manually save the csv (Control + S). Then reload via:
annotation_compact_revised <- read.csv(full_name)

annotation_continuous_revised <- ann_compact2continuous(annotation_compact_revised)
annotation_wfdb_revised <- ann_continuous2wfdb(annotation_continuous_revised)

# And replot:
plot_func(ecg_filter(signal),annotation_continuous_revised,pointsize= 0.5)

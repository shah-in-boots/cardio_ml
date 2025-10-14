# Add UIH Annotations -----------------------------------------------------
# Add UIH annotations
setwd('C:/Users/darre/OneDrive/Documents/UICOM Research/annotators/code')
source('annotator_prep_functions.R')
load('../deid_uih_ecgs.RData')

leads <- c('I','II','III','AVR','AVL','AVF','V1','V2','V3','V4','V5','V6')

# Run this once:
# for (i in seq_along(deid_uih_ecgs)) {
#   deid_uih_ecgs[[i]]$annotation <- setNames(vector("list", length(leads)), leads)
# }

files <- list.files(path = "../uih_predictions_shreyank/", full.names = TRUE)
files <- files[file.info(files)$isdir == FALSE]

for (lead in c(1:2,9:12)) {
  lead_files <- files[grepl(files,pattern=leads[lead])]
  for (sample in 1:length(lead_files)) {
    file_name <- lead_files[sample]
    number <- as.integer(sub("_.*", "", basename(file_name)))
    ann_compact <- read.csv(file_name)
    
    ann_continuous <- ann_compact2continuous(ann_compact)
    ann_wfdb <- ann_continuous2wfdb(ann_continuous)
    
    deid_uih_ecgs[[number]]$annotation[[lead]] <- ann_wfdb
  }
}

save(deid_uih_ecgs,file='../deid_uih_ecgs.RData')

# For each lead 
#   For each sample
#     transform annotation to wfdb format
#     add annotation to matrix


# update cluster ----------------------------------------------------------

# note: will need to add to uih_ecgs on cluster 
load('uih_ecgs.RData')
load('deid_uih_ecgs.RData')
for (i in seq_along(uih_ecgs)) {
  uih_ecgs[[i]]$annotation <- deid_uih_ecgs[[i]]$annotation
}
save(uih_ecgs,file='uih_ecgs.RData')
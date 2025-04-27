# Pick sinus samples
library(fs)
library(EGM)
library(stringr)
library(XML)
library(xml2)
library(dplyr)
options(wfdb_path = '/mmfs1/home/dseaney2/wfdb/bin')
source('../../ECG Segmentation/code/cluster_functions.R') # unused?

muse_path = '/mmfs1/projects/cardio_darbar_chi/common/data/muse/muse.log'
muse_log <- read.csv(muse_path)

wfdb_path = '/mmfs1/projects/cardio_darbar_chi/common/data/wfdb/wfdb.log'
wfdb_log <- read.csv(wfdb_path) 

# Generate random samples -----------------------------------------------------
size <- 100

uih_ecgs <- list(NA,c(size))
counter <- 1
row <- 1
while (row <= size) {
  idx <- sample(1:nrow(wfdb_log), 1, replace = FALSE)
  
  
  dir <- fs::path(base_path,dirname(wfdb_log$PATH[idx]))
  file <- wfdb_log$FILE_NAME[idx]
  
  muse_idx <- which((muse_log$FILE_NAME == wfdb_log$FILE_NAME[idx]))
  xml <- read_xml(fs::path(base_path, dirname(muse_log$PATH[muse_idx]), muse_log$FILE_NAME[muse_idx],ext = 'xml'))
  diagnosis <- xml_text(xml_find_all(xml, ".//DiagnosisStatement"))
  
  if (any(grepl(diagnosis,pattern='Sinus'))) {
    wfdb <- read_wfdb(record = file,record_dir = dir)
    sig <- wfdb$signal$I
    if (length(sig) == 5000) {
      uih_ecgs[row] <- list(wfdb)
      row <- row+1
    }
  }
  
  counter <- counter+1
  
  if (!row %% 10) {
    print(paste0('Added: 'row,', Iterated thru: 'counter))
  }
}

save_path <- '/mmfs1/projects/cardio_darbar_chi/common/cohorts/wes_ml/annotator_models/uih_ecgs.RData'
save(uih_ecgs,file = save_path)

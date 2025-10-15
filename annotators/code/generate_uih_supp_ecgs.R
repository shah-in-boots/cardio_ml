# Based on 'generate_uih_test_ecgs.R'. 
# Now included in annotator_prep_functions

# Pick sinus samples
library(fs)
library(EGM)
library(stringr)
library(XML)
library(xml2)
library(dplyr)
options(wfdb_path = '/mmfs1/home/dseaney2/wfdb/bin')
# source('../../ECG Segmentation/code/cluster_functions.R') # unused?

muse_path = '/mmfs1/projects/cardio_darbar_chi/common/data/muse/muse.log'
muse_log <- read.csv(muse_path)

wfdb_path = '/mmfs1/projects/cardio_darbar_chi/common/data/wfdb/wfdb.log'
wfdb_log <- read.csv(wfdb_path) 

base_path <-  "/mmfs1/projects/cardio_darbar_chi/common/"

# Generate random samples -----------------------------------------------------
size <- 10
dx_pattern <- 'Atrial fib|Afib'#'Sinus'

dir.create('temp')
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
  
  # Check if sinus rhythm
  if (any(grepl(diagnosis,pattern=dx_pattern))) {
    test <- read_wfdb(record = file, record_dir = dir)
    sig <- test$signal$I
    # Check if ECG is 500 Hz sampling
    if (length(sig) == 5000) {
      # Need to move ecgpuwave file to same folder as header and data file
      file.copy(from = fs::path(dir,file,ext='hea'), to = fs::path('temp'))
      file.copy(from = fs::path(dir,file,ext='dat'), to = fs::path('temp'))
      
      # Copy ecgpuwave ECG
      ecgpuwave_dir <- sub("wfdb", "ecgpuwave", dir)
      file.copy(from = fs::path(ecgpuwave_dir,file,ext='ecgpuwave'), to = fs::path('temp'))
      
      # Read ECG 
      wfdb <- read_wfdb(record = file, record_dir = 'temp', annotator = 'ecgpuwave')
      uih_ecgs[row] <- list(wfdb)
      
      row <- row+1
      file.remove(list.files('temp', full.names=TRUE))
    }
  }
  
  counter <- counter+1
  
  if (!row %% 10) {
    print(paste0('Added: ',row,', Iterated thru: ', counter))
  }
}
unlink("temp")
uih_ecgs_supp <- uih_ecgs

save_path <- '/mmfs1/projects/cardio_darbar_chi/common/cohorts/wes_ml/annotator_models/uih_ecgs_supp.RData'
save(uih_ecgs_supp,file = save_path)


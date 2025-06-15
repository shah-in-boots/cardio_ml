# Predict -----------------------------------------------------------------
library(keras)
library(tools)
library(EGM)
dir <- '../../AF_DM/wfdb/matched_controls/'
load('../models/model_log.RData')
source('annotator_prep_functions.R')
options(wfdb_path = 'wsl /usr/local/bin')

testfiles <- file_path_sans_ext(list.files(path = dir, pattern = '.hea')[1:100])
# for (i in 1:length(testfiles)) {
#   read_header(record = testfiles[i], record_dir = dir)
# }

row <- 41
lead <- 1

model <- load_model_tf(paste0('../models/',model_log$name[row]))

signal <- do.call(rbind, lapply(testfiles, function(file) {
  read_signal(record = file, record_dir = dir)$I
}))

filtered <- array(0,dim(signal))
for (i in 1:nrow(signal)) {
  filtered[i,] <- ecg_filter(signal[i,])
}

# Predict
predictions <- model %>% predict(filtered)
predictions_integer <- array(0,c(nrow(predictions),ncol(predictions))) 
for (i in 1: nrow(predictions)) {
  predictions_integer[i,] <- max.col(predictions[i,,])
}
#convert from dimension value 1,2,3,4 to 0,1,2,3
predictions_integer <- predictions_integer - 1


# Stats -------------------------------------------------------------------
prog_count <- array(NA,nrow(predictions_integer))
prog_count_revised <- array(NA,nrow(predictions_integer))

predictions_integer_revised <- array(0,dim(predictions_integer))
for (i in 1:nrow(predictions_integer)) {
  predictions_integer_revised[i,] <- fill_wave_gaps(predictions_integer[i,],20)
}

for (i in 1:nrow(predictions_integer)) {
  prog_count[i] <- check_ann_prog(annotation = predictions_integer[i,])
  prog_count_revised[i] <- check_ann_prog(predictions_integer_revised[i,])
}

which(prog_count_revised == 0)
# [1]   4   5   6   8  10  11  12  14  16  18  19  20  23  24  26  30  32  33  34  35
# [21]  40  42  43  44  47  49  50  51  52  54  55  58  59  60  63  64  65  67  69  70
# [41]  73  80  83  84  85  86  87  88  90  93  94  99 100

# Filter out terminal 250 indices?
# For plotting, could change problematic waves to a different color
# 

# plot --------------------------------------------------------------------
sample <- 73
raw <- plot_func(signal[sample,],predictions_integer[sample,])
filt <- plot_func(filtered[sample,],predictions_integer[sample,])

ann <- ann_continuous2wfdb(predictions_integer[sample,])
ann[ann$type %in% c('p','N','t'),]

  filt
# subplot(raw,filt,nrows = 2)

# Ideas for determining accuracy: compare to find_Rpeaks, and check p > QRS > t > p 
#   Rpeaks anchors the QRS (global), and p qrs t checks local progression
# Helper functions --------------------------------------------------------
check_ann_prog <- function(annotation) {
  ann <- ann_continuous2wfdb(annotation)
  progression <- ann$type[ann$type %in% c('p','N','t')]
  
  # Create logical checks
  check_p <- which(progression == "p")  # Indices of "p"
  check_N <- which(progression == "N")  # Indices of "N"
  check_t <- which(progression == "t")  # Indices of "t"
  
  # Ensure "p" is followed by "N"
  valid_pN <- all(progression[check_p[check_p < length(progression)] + 1] == "N")
  
  # Ensure "N" is followed by "t"
  valid_Nt <- all(progression[check_N[check_N < length(progression)] + 1] == "t")
  
  # Ensure "t" is followed by "p"
  valid_tp <- all(progression[check_t[check_t < length(progression)] + 1] == "p")
  
  # Final check
  if (valid_pN & valid_Nt & valid_tp) {
    # print("The progression follows the expected pattern.")
    output <- 1
  } else {
    # print("The progression does NOT follow the expected pattern.")
    output <- 0
  }
}

fill_wave_gaps <- function(annotation = predictions_integer, max_gap = 10) {
  # Find nonzero indices
  nonzero_indices <- which(annotation != 0)
  
  # Iterate through the nonzero indices
  for (i in seq_along(nonzero_indices)) {
    if (i == length(nonzero_indices)) break  # Stop if at last index
    
    current_val <- annotation[nonzero_indices[i]]
    next_idx <- nonzero_indices[i + 1]
    
    # Check if the next wave is of the same value and gap is ≤ max_gap
    if (annotation[next_idx] == current_val && (next_idx - nonzero_indices[i] <= max_gap)) {
      annotation[(nonzero_indices[i] + 1):(next_idx - 1)] <- current_val
    }
  }
  
  return(annotation)
}


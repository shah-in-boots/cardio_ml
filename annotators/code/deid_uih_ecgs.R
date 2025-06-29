load('uih_ecgs.RData')
for (i in 1:length(uih_ecgs)) {
  attributes(uih_ecgs[[i]]$header)$record_line$record_name <- 0 # single value
  uih_ecgs[[i]]$header$file_name <- array(0,12) # vector of 12
}
deid_uih_ecgs <- uih_ecgs
save(deid_uih_ecgs,file='deid_uih_ecgs.RData')
rm(uih_ecgs)

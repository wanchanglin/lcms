########### Script to process lcms data using xcms ########################################################################

###### Z.Hall, Aug 2014 ##################################################################################################


  spectra_dir <- "C:/Users/charidemoue/Desktop/mzML"
  
  lib_dir <- "C:/MSanalysis/libraries" # the folder where your library files are stored
  
  
  FWHM <- 3 # set approximate FWHM (in seconds) of chromatographic peaks
  
  
  snthresh = 5 # set the signal to noise threshold
  
  
  minfrac = 0.25 # minimum fraction of samples necessary for it to be a valid peak group 
  
  
  adducts = c(H = T, NH4 = T, Na = T, K = T, dH = F, Cl = F, OAc = F) # choose which adduct you want to include in library search
  
  
  ppm.annotate = 5   # choose ppm for annotations
  
  
  ionisation_mode <- "positive"  # either "positive" or "negative" - determines which classes of lipids to search for
  
  
  source("C:/MSanalysis/lcmsprocessing/lcms.R") # source your R script


####################################################################################################################################

  
  

########### Script to process lcms data using xcms ########################################################################

###### Z.Hall, Aug 2014 ##################################################################################################


spectra_dir <- "C:/mzML" # where your mzML files are stored


lib_dir <- "C:/MS analysis/libraries"                # the folder where your library files are stored


FWHM <- 3 # set approximate FWHM (in seconds) of chromatographic peaks


snthresh = 5 # set the signal to noise threshold


minfrac = 0.25 # minimum fraction of samples necessary for it to be a valid peak group 


adducts = c(H = F, NH4 = F, Na = F, K = F, dH = T, Cl = T, OAc = T) # choose which adduct you want to include in library search


ppm.annotate = 5   # choose ppm for annotations


ionisation_mode <- "negative"  # either "positive" or "negative" - determines which classes of lipids to search for


source("C:/MSanalysis/lcmsprocessing") # source your R script


####################################################################################################################################




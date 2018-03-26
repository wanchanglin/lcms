## wl-24-08-2017, Thu:

## ######### Script to process lcms data using xcms ###############
## Z.Hall, Aug 2014

setwd("C:/R_lwc/lcms")
source("lcms_CDF.R") # source your R script

library(xcms)
library(reshape2)

home_dir    <- paste(getwd(), "/", sep="")
lib_dir     <- paste(home_dir,"libraries",sep="") # library files
spectra_dir <- paste(home_dir,"CDF",sep="") # spectra files (mzML files)

FWHM <- 3 # set approximate FWHM (in seconds) of chromatographic peaks
snthresh = 5 # set the signal to noise threshold
minfrac = 0.25 # minimum fraction of samples necessary for it to be a valid peak group
## choose which adduct you want to include in library search
adducts = c(H = T, NH4 = T, Na = T, K = F, dH = F, Cl = F, OAc = F)
ppm.annotate = 15   # choose ppm for annotations
## either "positive" or "negative" - determines which classes of lipids to search for
ionisation_mode <- "positive"

## read in library files
setwd(lib_dir)
read <- read.csv('lib_FA.csv', sep=",", header=T);lookup_FA <- read[,2:4]; row.names(lookup_FA) <- read[,1]
read <- read.csv('lib_class.csv', sep=",", header=T);lookup_lipid_class <- read[,2:3]; row.names(lookup_lipid_class) <- read[,1]
read <- read.csv('lib_element.csv', sep=",", header=T);lookup_element <- read[,2:3]; row.names(lookup_element) <- read[,1]
read <- read.csv('lib_modification.csv', sep=",", header=T);lookup_mod <- read[,2:ncol(read)]; row.names(lookup_mod ) <- read[,1]

FA_expt <-list('10','12','14','15','16','16.1','17','17.1','18','18.1','18.2','18.3','20.3','20.4','20.5','21','22','22.5','22.6','24.1')

## fatty acids to be included in library

## place converted CDF files in working directory
setwd(spectra_dir)
path <- paste(getwd())
## use this if you want to look at all files in directory
files <- list.files(path, recursive=T, full.names=F, pattern=".CDF")

## use this if you want to exclude files using an exlcusion list

##  list_1 <- list.files(path, recursive=T, full.names=F, pattern=".mzML")
##  csv <- read.csv("exclusion_list.csv", header=F)
##  list_2 <- (as.vector(as.character(csv[1:nrow(csv),1])))
##  list_3 <- list_1[- which(list_1 %in% list_2)]
##  files <- list_3

profmethod = "binlin"
## use either "bin" (better for centroid, default), "binlin" (better for profile)

## ################### XCMS and some post-processing #####################
if (ionisation_mode == "positive"){
  ## lipid classes to include - will need to change for negative ion mode
  sel.class<- c(
      T, #TG
      T, #DG
      T, #PC
      F, #PA
      T, #PE
      T, #PS
      F, #PG
      F, #PI
      F, #PIP
      F, #PIP2
      F, #PIP3
      T, #LysoPC
      T, #DG-H20
      T, #CE
      F, #FFA
      T, #SM
      T  #Cer
  )
}

if (ionisation_mode =="negative"){
  sel.class<- c(
      F, #TG
      F, #DG
      T, #PC
      T, #PA
      T, #PE
      T, #PS
      T, #PG
      T, #PI
      T, #PIP
      T, #PIP2
      T, #PIP3
      F, #LysoPC
      F, #DG-H20
      F, #CE
      T, #FFA
      F, #SM
      F  #Cer
  )
}

## =========================================================================
dbase <- makelibrary(sel.class, fixed=F, fixed_FA, lookup_lipid_class,
                     lookup_FA, lookup_element)
## xcms scripts for peak picking and integration
xset <- xcmsSet(files, method="matchedFilter", step = 0.1, sigma=FWHM/2.3548,
                snthresh=snthresh, profmethod=profmethod)
## feature recognition and convert group of files to xcms raw object

xset <- group(xset) # groups peaks, temporarily
## xset2 <- retcor(xset, family = "symmetric", plottype = "mdevden",
##                 smooth = "loess", span = 0.2, missing = 5, extra = 1)
## corrects retention times

xset2 <- retcor(xset, method="obiwarp", profStep=0.1, plottype="deviation")
## different method for retention time correction - more robust but slower

## check for retention time outliers and exclude samples if necessary
rmsd <- sapply(1:length(xset2@rt$corrected),
               function(x) {
                 sqrt(sum((xset2@rt$raw[[x]]-xset2@rt$corrected[[x]]) ** 2))
               })
plot(rmsd)

rmsd <- cbind(files, rmsd)
write.csv(rmsd, "RT_rmsd.csv")

xset2 <- group(xset2, bw = 5,  minfrac = minfrac, mzwid = 0.025)
## groups data based on corrected retention times

## =========================================================================
## make output based on unique peak groups (mz and RT groups)
names <- groupnames(xset2, rtdec=2, mzdec=4)  # define decimal places for mz and RT values
names <- substr(names, 2, 18)

out <- groupval(xset2, method="maxint", value="maxo", intensity="maxo")
## if multiple peaks choose one with highest intensity based on maxo (peak hieght raw)

rownames(out) <- names
out[is.na(out)] <- 0 # replace NAs with 0
write.csv(out, "xcms_peak_height_raw.csv")

out <- groupval(xset2, method="maxint", value="into", intensity="maxo")
rownames(out) <- names
out[is.na(out)] <- 0
write.csv(out, "xcms_peak_area_raw.csv")

## =========================================================================
## deisotoping - working with raw peak area - can also change "out" to use
## for peak height or the fitered results (value = into/intf/maxo/maxf)
out <- groupval(xset2, method="maxint", value="into", intensity="maxo")
## if multiple peaks choose one with highest intensity based on maxo (peak
## hieght raw)

rownames(out) <- names
out[is.na(out)] <- 0 # replace NAs with 0

names_split <- colsplit(as.vector(as.character(names)), "T", c("mz", "RT"))

spectra <- cbind(names_split$mz, out[,1], "", "")
colnames(spectra) <- c("mz.obs", "intensity", "isotope", "modification" )

deisotoped <- deisotoping(ppm=5, no_isotopes=2, prop.1=0.9, prop.2=0.5,
                          spectra=spectra)

indices <- match(rownames(deisotoped), rownames(out))
deisotoped_out <- out[indices, ]

names_deiso <- rownames(deisotoped_out)
names_deiso <- colsplit(as.vector(as.character(names_deiso)), "T", c("mz", "RT"))
deisotoped_out <- cbind(names_deiso$mz, names_deiso$RT, deisotoped_out[,1:ncol(deisotoped_out)])

write.csv(deisotoped_out, "xcms_peak_area_raw_deisotoped.csv")

## =========================================================================
## annotating

annotated<-annotating(deisotoped, adducts, ppm.annotate, dbase)

final_out <- cbind(annotated, deisotoped_out)

write.csv(final_out, "final_annotated_desiotoped.csv")

## =========================================================================
## normalising based on TIC

data <- read.csv("final_annotated_desiotoped.csv")
rownames(data) <- data[,1]
sums <- as.vector(colSums(data[,5:ncol(data)], na.rm=T))
factor <- sums/mean(sums)
TIC_normalised <- cbind(data[,2:4], t(t(data[,5:ncol(data)])/factor))

write.csv(TIC_normalised, "final_annotated_deisotoped_norm.csv")
##hist(factor, breaks=100)

## =========================================================================
## fill in missing peaks by integrating noise
xset3 <- fillPeaks(xset2) # fills in missing peaks

names <- groupnames(xset3, rtdec=2, mzdec=4)  # define decimal places for mz and RT values
names <- substr(names, 2, 18)

out <- groupval(xset3, method="maxint", value="into", intensity="maxo")
## if multiple peaks choose one with highest intensity based on maxo (peak hieght raw)

rownames(out) <- names
names_split <- colsplit(as.vector(as.character(names)), "T", c("mz", "RT"))

## =========================================================================
## deisotope

spectra <- cbind(names_split$mz, out[,1], "", "")
colnames(spectra) <- c("mz.obs", "intensity", "isotope", "modification" )

deisotoped <- deisotoping(ppm=5, no_isotopes=2, prop.1=0.9, prop.2=0.5,
                          spectra=spectra)

indices <- match(rownames(deisotoped), rownames(out))
deisotoped_out <- out[indices, ]

names_deiso <- rownames(deisotoped_out)
names_deiso <- colsplit(as.vector(as.character(names_deiso)), "T", c("mz", "RT"))
deisotoped_out <- cbind(names_deiso$mz, names_deiso$RT, deisotoped_out[,1:ncol(deisotoped_out)])

write.csv(deisotoped_out, "xcms_fill_peak_area_raw_deisotoped.csv")


## =========================================================================
## annotate

annotated <-annotating(deisotoped, adducts, ppm.annotate, dbase)

final_out <- cbind(annotated, deisotoped_out)

write.csv(final_out, "final_annotated_desiotoped_fill.csv")

## =========================================================================
## normalising based on TIC

data <- read.csv("final_annotated_desiotoped_fill.csv")
rownames(data) <- data[,1]
sums <- as.vector(colSums(data[,5:ncol(data)], na.rm=T))
factor <- sums/mean(sums)
TIC_normalised <- cbind(data[,2:4], t(t(data[,5:ncol(data)])/factor))

write.csv(TIC_normalised, "final_annotated_deisotoped_fill_norm.csv")

## hist(factor, breaks=100)


#########################################################################

## plot extracted ion chromatograms for each peak group (corrected RT) if
## desired, can choose specific samples and/or mz

## xic.corr <- getEIC(xset3, rt = "corrected", groupidx = 1:nrow(xset3@groups))
## plot(xic.corr, xset3, groupidx=1:nrow(xset3@groups), sampleidx=1:length(files))

#########################################################################





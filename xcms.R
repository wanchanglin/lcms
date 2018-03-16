## wl-12-03-2018, Mon: apply xcms for -LC-MS 
## to-do
##  1.) Is is neccesary to use 'fillPeaks'?
##  2.) other usages

library(xcms)
## library(faahKO)
library(CAMERA)

FWHM       <- 3 # set approximate FWHM (in seconds) of chromatographic peaks
snthresh   <- 5 # set the signal to noise threshold
minfrac    <- 0.25 # minimum fraction of samples necessary for it to be a valid peak group
profmethod <- "binlin"  # use either "bin" (better for centroid, default), "binlin" (better for profile)

path <- "C:/R_lwc/20180309_EC_SM_AM12_PartIV/mzML"
files <- list.files(path, pattern="mzML",recursive = F, full.names = TRUE)

## ========================================================================
## Run xcms
if (F){
  mt:::tic()
  xset <- xcmsSet(files, method="matchedFilter", step = 0.1,
                  sigma=FWHM/2.3548, snthresh=snthresh, 
                  profmethod=profmethod)

  xset <- group(xset) ## slotNames(xset)

  ## wl-15-03-2018, Thu: Possible memory problem?
  xset <- retcor(xset, method="obiwarp", profStep=0.1, plottype="deviation") 

  xset <- group(xset, bw = 5,  minfrac = minfrac, mzwid = 0.025) 
  ## lwc-05-11-2013: group has three emthods: group.density (default), 
  ##  group.mzClust and group.nearest.

  ## xset <- fillPeaks(xset) 
  ## Note: need mzML files to fill in missing peaks

  mt:::toc() ## wl-12-03-2018, Mon: 15 minutes

  save(xset,file="./test-data/xset.RData")
} else {
  load("./test-data/xset.RData")
}

## ========================================================================
## lwc-08-10-2013: Get peak lists. 
if (F){ ## No need for peaklist
  peakmat <- xcms::peaks(xset)   
  ## Note: there is also function 'peaks' from package mzR.
  dim(peakmat)
}

grpmat             <- groups(xset)    ## is object@groups ##dim(grpmat)
values             <- groupval(xset, method="medret", value="into")
peaklist           <- cbind(grpmat, values)
rownames(peaklist) <- NULL
## lwc-08-10-2013: 'groupval' arguments: 
## 'method': c("medret", "maxint")
## value': c("into","intb","maxo"). 

write.csv(peaklist,file=paste(cdfpath,"/",'xcms_peak.csv',sep=""))

nam <- groupnames(xset, rtdec=2, mzdec=4)  
nam <- substr(nam, 2, 18)

## ------------------------------------------------------------------ 
## or get peak list after annotation by CAMERA
## Note: need the original mzML files
## xsa         <- annotate(xset, cor_eic_th=0)
## peaklist.2  <- getPeaklist(xsa)

## or get peak list directly from CAMERA's hidden function 
## peaklist.1 <- CAMERA:::getPeaks_selection(xset)

## ==========================================================================
## Create complete feature table from CAMERA
## xs     - xcmsSet object
## method - groupval parameter method
## value  - groupval parameter method
## ------------------------------------------------------------------ 
## wl-16-03-2018, Fri: need to change this function to allow more 
## arguments passing.
## ------------------------------------------------------------------ 
getPeaks_selection <- function(xs, method="medret", value="into"){
  if (!class(xs) == "xcmsSet") {
    stop ("Parameter xs is no xcmsSet object\n")
  }

  ## Testing if xcmsSet is grouped
  if (nrow(xs@groups) > 0 && length(xs@filepaths) > 1) {
     ## get grouping information
     groupmat <- groups(xs)
     ## generate data.frame for peaktable
     ts <- data.frame(cbind(groupmat, groupval(xs, method=method, value=value)), row.names = NULL)
     ##rename column names
     cnames <- colnames(ts)
     if (cnames[1] == 'mzmed') {
       cnames[1] <- 'mz'
     } else {
       stop ('Peak information ?!?')
     }
     if (cnames[4] == 'rtmed') {
       cnames[4] <- 'rt'
     } else {
       stop ('Peak information ?!?')
     }
     colnames(ts) <- cnames
  } else if (length(sampnames(xs)) == 1) { #Contains only one sample?
    ts <- xs@peaks
  } else {
    stop ('First argument must be a xcmsSet with group information or contain only one sample.')
  }

  return(as.matrix(ts))
}

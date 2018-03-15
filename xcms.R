## wl-12-03-2018, Mon: apply xcms for -LC-MS 
library(xcms)
library(CAMERA)

FWHM       <- 3 # set approximate FWHM (in seconds) of chromatographic peaks
snthresh   <- 5 # set the signal to noise threshold
minfrac    <- 0.25 # minimum fraction of samples necessary for it to be a valid peak group
profmethod <- "binlin"  # use either "bin" (better for centroid, default), "binlin" (better for profile)

path <- "C:/R_lwc/20180309_EC_SM_AM12_PartIV/mzML"
files <- list.files(path, pattern="mzML",recursive = F, full.names = TRUE)

if (F){
  mt:::tic()
  ## xcms scripts for peak picking and integration
  xset  <- xcmsSet(files, method="matchedFilter", step = 0.1, 
                   sigma=FWHM/2.3548, snthresh=snthresh, profmethod=profmethod)
  xset  <- group(xset)
  xset2 <- retcor(xset, method="obiwarp", profStep=0.1, plottype="deviation") 
  mt:::toc() ## wl-12-03-2018, Mon: 15 minutes

  xset2 <- group(xset2, bw = 5,  minfrac = minfrac, mzwid = 0.025) # groups data based on corrected retention times
  ## lwc-05-11-2013: group has three emthods: group.density (default), 
  ##  group.mzClust and group.nearest.

  save(xset2,file="xset2.RData")
} else {
  load("xset2.RData")
}

## --------------------------------------------------------------------------
## lwc-08-2013: We are only interested in getting peak lists. So do not need 
## annotation. Can directly get peak lists by:
peaklist <- CAMERA:::getPeaks_selection(xset2)

## xsa         <- annotate(grp, cor_eic_th=0)
## peaklist.1  <- getPeaklist(xsa)
write.csv(peaklist,file=paste(cdfpath,"/",'xcms_kidneycortex_081013_neg.csv',sep=""))

## ==========================================================================
## Create complete feature table from CAMERA
## xs     - xcmsSet object
## method - groupval parameter method
## value  - groupval parameter method
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

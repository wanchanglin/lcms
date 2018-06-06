## wl-12-03-2018, Mon: apply xcms for LC-MS 
## wl-19-03-2018, Mon: Use Zoe's parameters seeting for 'xcms'
## wl-20-03-2018, Tue: use 'peakTable' to get peak list

library(xcms)

## file path of mzML files
path  <- "C:/R_lwc/data/20180309_EC_SM_AM12_PartIV/mzML"
files <- list.files(path, pattern="mzML",recursive = F, full.names = TRUE)

## parameters qor 'xcmsSet' (based on Zoe Hall's setting)
FWHM       <- 3 # set approximate FWHM (in seconds) of chromatographic peaks
snthresh   <- 5 # set the signal to noise threshold
## minimum fraction of samples necessary for it to be a valid peak group
minfrac    <- 0.25 
## use either "bin" (better for centroid, default), "binlin" (better for
## profile)
profmethod <- "binlin"  

## ========================================================================
## Run xcms
if (F) {
  xset <- xcmsSet(files, method="matchedFilter", step = 0.1, 
                  sigma=FWHM/2.3548, snthresh=snthresh, 
                  profmethod=profmethod)

  xset <- group(xset) ## slotNames(xset)

  ## wl-15-03-2018, Thu: Possible memqry problem?
  xset <- retcor(xset, method="obiwarp", profStep=0.1, plottype="deviation") 

  xset <- group(xset, bw = 5,  minfrac = minfrac, mzwid = 0.025) 
  ## lwc-05-11-2013: group has three emthods: group.density (default), 
  ##  group.mzClust and group.nearest.

  xset <- fillPeaks(xset) 
  ## Note: need mzML files to fill in missing peaks

  save(xset,file="./test-data/xset.RData")
} else {
  load("./test-data/xset.RData")
}

## ========================================================================
## lwc-08-10-2013: Get peak lists. 
peakmat   <- xcms::peaks(xset)   
## Note: two arguments in 'groupval', 'value' and 'intensity' will use this
## matrix's intensity columns

grpmat    <- groups(xset)    ## is object@groups ##dim(grpmat)
## values <- groupval(xset, method="medret", value="into")
values    <- groupval(xset, method="maxint", value="into", intensity="maxo")
## values.1  <- groupval(xset, method="maxint", value="maxo", intensity="maxo")
## wl-20-03-2018, Tue: it seems that only 'value' to decide the intensity to
## be returned. For details, see source code of 'groupval'

peaklist  <- data.frame(cbind(grpmat, values),row.names = NULL)
colnames(peaklist) <-  gsub("mzmed","mz",colnames(peaklist))
colnames(peaklist) <-  gsub("rtmed","rt",colnames(peaklist))

## wl-20-03-2018, Tue: or use 'peakTable' directly.
## peaklist  <- peakTable(xset,method="maxint", value="into", intensity="maxo")

## keep only mz and rt in peak list
peaklist  <- subset(peaklist, select=-c(mzmin,mzmax,rtmin,rtmax,npeaks,mzML))


## ========================================================================
## ------------------------------------------------------------------ 
## or get peak list after annotation by CAMERA
## Note: need the original mzML files
## library(CAMERA)
## xsa         <- annotate(xset, cor_eic_th=0)
## peaklist.1  <- getPeaklist(xsa)

## or get peak list directly from CAMERA's hidden function 
## peaklist.2 <- CAMERA:::getPeaks_selection(xset)


###########################################################################
## xcms::peakTable
## wl-20-03-2018, Tue: Note that the dot arguments sgould be 'groupval'
## arguments.
setMethod("peakTable", "xcmsSet", function(object, filebase = character(), ...) {

    if (length(sampnames(object)) == 1) {
        return(object@peaks)
    }

    if (nrow(object@groups) < 1) {
        stop ('First argument must be an xcmsSet with group information or contain only one sample.')
    }

    groupmat <- groups(object)


    if (! "value" %in% names(list(...))) {
        ts <- data.frame(cbind(groupmat,groupval(object, value="into",  ...)), row.names = NULL)
    } else {
        ts <- data.frame(cbind(groupmat,groupval(object, ...)), row.names = NULL)
    }

    cnames <- colnames(ts)

    if (cnames[1] == 'mzmed') {
        cnames[1] <- 'mz'
    } else {
        stop ('mzmed column missing')
    }
    if (cnames[4] == 'rtmed') {
        cnames[4] <- 'rt'
    } else {
        stop ('mzmed column missing')
    }

    colnames(ts) <- cnames

    if (length(filebase))
        write.table(ts, paste(filebase, ".tsv", sep = ""), quote = FALSE, 
                    sep = "\t", col.names = NA)

    ts
})

###########################################################################
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

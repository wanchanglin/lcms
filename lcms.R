#' wl-12-03-2018, Mon: Commence
#' wl-15-03-2018, Thu: Tidy R codes
#' wl-19-03-2018, Mon: Use Zoe's parameters setting for 'xcms'
#' wl-20-03-2018, Tue: Use 'peakTable' to get peak list
#' wl-21-03-2018, Wed: Major changes
#' wl-23-03-2018, Fri: Test and debug deisotoping and annotating
#' wl-26-03-2018, Mon: Minor changes
#' wl-17-05-2019, Fri: modification for Galaxy

## ==== General settings ====
rm(list = ls(all = T))

#' flag for command-line use or not. If false, only for debug interactively.
com_f <- F

#' galaxy will stop even if R has warning message
options(warn = -1) #' disable R warning. Turn back: options(warn=0)

#' ------------------------------------------------------------------------
#' Setup R error handling to go to stderr
#' options( show.error.messages=F, error = function (){
#'   cat( geterrmessage(), file=stderr() )
#'   q( "no", 1, F )
#' })

#' we need that to not crash galaxy with an UTF8 error on German LC settings.
loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")

suppressPackageStartupMessages({
  library(optparse)
  library(WriteXLS)
  library(xcms)
})

#' wl-28-08-2018, Tue: Convert a string seperated by comma into character vector
str_vec <- function(x) {
  x <- unlist(strsplit(x, ","))
  x <- gsub("^[ \t]+|[ \t]+$", "", x) #' trim white spaces
}

## ==== Command line or interactive setting ====
if (com_f) {

  #' -----------------------------------------------------------------------
  #' Setup home directory
  #' wl-24-11-2017, Fri: A dummy function for the base directory. The reason
  #' to write such a function is to keep the returned values by
  #' 'commandArgs' with 'trailingOnly = FALSE' in a local environment
  #' otherwise 'parse_args' will use the results of
  #' 'commandArgs(trailingOnly = FALSE)' even with 'args =
  #' commandArgs(trailingOnly = TRUE)' in its argument area.
  func <- function() {
    argv <- commandArgs(trailingOnly = FALSE)
    path <- sub("--file=", "", argv[grep("--file=", argv)])
  }
  #' prog_name <- basename(func())
  tool_dir <- paste0(dirname(func()), "/")

  option_list <-
    list(
      make_option(c("-v", "--verbose"),
        action = "store_true", default = TRUE,
        help = "Print extra output [default]"
      ),
      make_option(c("-q", "--quietly"),
        action = "store_false",
        dest = "verbose", help = "Print little output"
      ),

      #' -------------------------------------------------------------------
      #' input
      make_option("--mzxml_file",
        type = "character",
        help = "mzXML/ mzML file directory or full file list seperated by comma"
      ),
      make_option("--targ_file",
        type = "character",
        help = "Lipid target list with columns of m/z and lipid name"
      ),
      make_option("--samp_name",
        type = "character", default = "",
        help = "Sample names. Default is the names of mz XML file"
      ),
      make_option("--rt_low",
        type = "double", default = 20.0,
        help = "Start time"
      ),
      make_option("--rt_high",
        type = "double", default = 60.0,
        help = "End time"
      ),
      make_option("--mz_low",
        type = "double", default = 200.0,
        help = "Start m/z"
      ),
      make_option("--mz_high",
        type = "double", default = 1200.0,
        help = "End m/z"
      ),
      make_option("--hwidth",
        type = "double", default = 0.01,
        help = "m/z window size/height for peak finder"
      ),

      #' output files (Excel)
      make_option("--sign_file",
        type = "character", default = "signals.tsv",
        help = "Save peak signals (peak table)"
      ),
      make_option("--devi",
        type = "logical", default = TRUE,
        help = "Return m/z deviation results or not"
      ),
      make_option("--devi_file",
        type = "character", default = "deviations.tsv",
        help = "Save m/z deviations"
      ),
      make_option("--indi",
        type = "logical", default = TRUE,
        help = "Return each sample's signal and m/z deviation or not"
      ),
      make_option("--indi_file",
        type = "character", default = "sam_indi.xlsx",
        help = "Save individual sample's signal and m/z deviation in Excel"
      )
    )

  opt <- parse_args(
    object = OptionParser(option_list = option_list),
    args = commandArgs(trailingOnly = TRUE)
  )
  print(opt)
} else {
  tool_dir <- "C:/R_lwc/lcms/"         #' for windows
  #' tool_dir <- "~/my_galaxy/lcms/" #' for linux. must be case-sensitive
  opt <- list(
    #' input

     #' mzxml_file = paste(paste0(tool_dir, "test-data/lcms_neg/ZH_180918_mann_neg_001.mzML"),
     #'                    paste0(tool_dir, "test-data/lcms_neg/ZH_180918_mann_neg_002.mzML"),
     #'                    paste0(tool_dir, "test-data/lcms_neg/ZH_180918_mann_neg_003.mzML"),
     #'                    paste0(tool_dir, "test-data/lcms_neg/ZH_180918_mann_neg_004.mzML"),
     #'                    paste0(tool_dir, "test-data/lcms_neg/ZH_180918_mann_neg_005.mzML"),
     #'                    paste0(tool_dir, "test-data/lcms_neg/ZH_180918_mann_neg_006.mzML"),
     #'                    paste0(tool_dir, "test-data/lcms_neg/ZH_180918_mann_neg_007.mzML"),
     #'                    sep = ","
     #'                    ),
    
    mzxml_file = paste(paste0(tool_dir, "test-data/lcms_pos")),

    #' parameters for 'xcmsSet' (based on Zoe Hall's setting)
    FWHM = 3, # set approximate FWHM (in seconds) of chromatographic peaks
    snthresh = 5, # set the signal to noise threshold
    #' minimum fraction of samples necessary for it to be a valid peak group
    minfrac = 0.25,
    #' use either "bin" (better for centroid, default), "binlin" (better for
    #' profile)
    profmethod = "binlin",

    ppm.annotate = 15, # choose ppm for annotations
    ionisation_mode = "positive", # either "positive" or "negative"
    adducts = c(H = T, NH4 = F, Na = T, K = T, dH = F, Cl = F, OAc = F)

    #' Output
    #' sign_file = paste0(tool_dir, "test-data/res_dimsp/signals.tsv"),
    #' devi = TRUE,
    #' devi_file = paste0(tool_dir, "test-data/res_dimsp/deviations.tsv"),
    #' indi = TRUE,
    #' indi_file = paste0(tool_dir, "test-data/res_dimsp/individuals.xlsx")
  )
}

suppressPackageStartupMessages({
  source(paste0(tool_dir, "lcms_func.R"))
})

lib_dir <- paste0(tool_dir, "libraries/")


#' --------------------------------------------------------------------
#' if (T) {
#'   #' wl-23-03-2018, Fri: from 'annotate' of 'massPix'
#'   if (ionisation_mode == "positive") {
#'     adducts <- c(H = T, NH4 = F, Na = T, K = T, dH = F, Cl = F, OAc = F)
#'   } else {
#'     adducts <- c(H = F, NH4 = F, Na = F, K = F, dH = T, Cl = T, OAc = F)
#'   }
#' } else {
#'   #' wl-23-03-2018, Fri: this the default value for 'annotating'
#'   adducts <- c(H = T, NH4 = T, Na = T, K = F, dH = F, Cl = F, OAc = F)
#' }

## ==== Main process ====

#' process multiple input files seperated by comma
#' wl-04-03-2019, Mon: add file directory option. Note that it is not for
#' galaxy.
if (dir.exists(opt$mzxml_file)) {   ## file directory
  opt$mzxml_file <- list.files(opt$mzxml_file, pattern = "mzml|mzxml", 
                               ignore.case = T, recursive = F, 
                               full.names = TRUE)
} else {  ## multiple files
  opt$mzxml_file <- str_vec(opt$mzxml_file)
} 

#' ========================================================================
#' Run xcms
if (T) {
  xset <- xcmsSet(opt$mzxml_file,
    method = "matchedFilter", step = 0.1,
    sigma = opt$FWHM / 2.3548, snthresh = opt$snthresh,
    profmethod = opt$profmethod
  )

  xset <- group(xset) #' slotNames(xset)

  #' corrects retention times
  xset <- retcor(xset, method = "obiwarp", profStep = 0.1, 
                 plottype = "deviation")
  #' wl-15-03-2018, Thu: Possible memory problem?

  xset <- group(xset, bw = 5, minfrac = opt$minfrac, mzwid = 0.025)
  #' lwc-05-11-2013: group has three methods: group.density (default),
  #'  group.mzClust and group.nearest.

  xset <- fillPeaks(xset)
  #' Note: need mzML files to fill in missing peaks

  save(xset, file = "./test-data/xset_pos.RData")
} else {
  load("./test-data/xset_pos.RData")
}

#' ========================================================================
#' Get peak lists.
peakmat <- xcms::peaks(xset)
#' Note: two arguments in 'groupval', 'value' and 'intensity' will use this
#' matrix's intensity columns

grpmat <- groups(xset) #' is object@groups ##dim(grpmat)
#' values <- groupval(xset, method="medret", value="into")
values <- groupval(xset, method = "maxint", value = "into", intensity = "maxo")
#' values.1  <- groupval(xset, method="maxint", value="maxo", intensity="maxo")
#' wl-20-03-2018, Tue: it seems that only 'value' to decide the intensity to
#' be returned. For details, see source code of 'groupval'

peaklist <- data.frame(cbind(grpmat, values), row.names = NULL)
colnames(peaklist) <- gsub("mzmed", "mz", colnames(peaklist))
colnames(peaklist) <- gsub("rtmed", "rt", colnames(peaklist))

#' peaklist <- peakTable(xset, method = "maxint", value = "into", intensity = "maxo")
#' round mz and rt
#' peaklist <- transform(peaklist, mz = round(mz, 4), rt = round(rt, 2))

#' keep only mz and rt in peak list
#' wl-13-05-2019, Mon: need to change 
peaklist <- subset(peaklist, select = -c(mzmin, mzmax, rtmin, rtmax, npeaks, mzML))


#' =======================================================================
#' Deisotoping

#' spectra <- as.matrix(peaklist[,c("mz","X1")])
spectra <- as.matrix(peaklist[, c(1, 3)]) #' mz and intensity of the 1st sample
spectra <- cbind(spectra, "", "")
colnames(spectra) <- c("mz.obs", "intensity", "isotope", "modification")

deisotoped <- deisotoping(ppm = 5, no_isotopes = 2, prop.1 = 0.9, 
                          prop.2 = 0.5, spectra = spectra)

#' ========================================================================
#' Annotating

#' library files
read <- read.csv(paste(lib_dir, "lib_FA.csv", sep = "/"), sep = ",", header = T)
lookup_FA <- read[, 2:4]
row.names(lookup_FA) <- read[, 1]

read <- read.csv(paste(lib_dir, "lib_class.csv", sep = "/"), sep = ",", header = T)
lookup_lipid_class <- read[, 2:3]
row.names(lookup_lipid_class) <- read[, 1]

read <- read.csv(paste(lib_dir, "lib_element.csv", sep = "/"), sep = ",", header = T)
lookup_element <- read[, 2:3]
row.names(lookup_element) <- read[, 1]

read <- read.csv(paste(lib_dir, "lib_modification.csv", sep = "/"), sep = ",", header = T)
lookup_mod <- read[, 2:ncol(read)]
row.names(lookup_mod) <- read[, 1]

#' make library
dbase <- makelibrary(ionisation_mode, lookup_lipid_class, lookup_FA, lookup_element)

#' annotating
annotated <- annotating(deisotoped, adducts, ppm.annotate, dbase)

#' update peaklist
indices <- match(deisotoped[, "mz.obs"], peaklist[, "mz"])
tmp <- peaklist[indices, ]
rownames(tmp) <- NULL
final_peak <- cbind(annotated, tmp)

#' =======================================================================
#' save results
save(final_peak, deisotoped, annotated, file = "./test-data/peak.RData")

write.csv(final_peak,
  file = "./test-data/lc_ms_data.csv",
  row.names = FALSE, quote = FALSE
)

#' write.table(final_peak, file="./test-data/lc_ms_data.tsv", sep = "\t",
#'             row.names = FALSE, quote = FALSE)


#' =======================================================================
#' Some original codes
#' library(reshape2) ## for colsplit
if (F) {

  #' Deisotoping
  names <- groupnames(xset, rtdec = 2, mzdec = 4)
  names <- substr(names, 2, 18)
  out <- groupval(xset, method = "maxint", value = "into", intensity = "maxo")
  rownames(out) <- names
  names_split <- colsplit(as.vector(as.character(names)), "T", c("mz", "RT"))

  spectra <- cbind(names_split$mz, out[, 1], "", "")
  colnames(spectra) <- c("mz.obs", "intensity", "isotope", "modification")

  deisotoped <- deisotoping(
    ppm = 5, no_isotopes = 2, prop.1 = 0.9, prop.2 = 0.5,
    spectra = spectra
  )

  indices <- match(rownames(deisotoped), rownames(out))
  deisotoped_out <- out[indices, ]

  names_deiso <- rownames(deisotoped_out)
  names_deiso <- colsplit(as.vector(as.character(names_deiso)), "T", c("mz", "RT"))
  deisotoped_out <- cbind(
    names_deiso$mz, names_deiso$RT,
    deisotoped_out[, 1:ncol(deisotoped_out)]
  )

  #' Annotating
  annotated <- annotating(deisotoped, adducts, ppm.annotate, dbase)
  final_out <- cbind(annotated, deisotoped_out)
}

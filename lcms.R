#' wl-12-03-2018, Mon: Commence
#' wl-15-03-2018, Thu: Tidy R codes
#' wl-19-03-2018, Mon: Use Zoe's parameters setting for 'xcms'
#' wl-20-03-2018, Tue: Use 'peakTable' to get peak list
#' wl-21-03-2018, Wed: Major changes
#' wl-23-03-2018, Fri: Test and debug deisotoping and annotating
#' wl-26-03-2018, Mon: Minor changes
#' wl-17-05-2019, Fri: modification for Galaxy
#' wl-21-05-2019, Tue: use 'makelibrary' from 'massPix'
#' wl-22-05-2019, Wed:
#'   - compare deisotoping and annotation with 'massPix'
#'   - move adducts into annotate function
#'   - run in interactive mode
#' wl-23-05-2019, Thu: Fix some bugs.
#' wl-24-05-2019, Fri: debug
#' wl-28-05-2019, Tue: retcor
#'  - 'retcor' consumes a large quantity of memory by 'profStep' and
#'    disable its plot.
#'  - use average of all samples' intensity for deisotyping
#' wl-30-05-2019, Thu: add 'samp_name' especially for galaxy
#' wl-14-10-2019, Mon: fix a bug for planemo test only. 'planemo test' will
#'  put 'test-data' in an random directory. To avoid this, use 'samp_name'
#'  directly in 'peaklist'.
#' wl-26-08-2020, Wed: Review

## ==== General settings ====
rm(list = ls(all = T))
#' set.seed(123)

#' flag for command-line use or not. If false, only for debug interactively.
com_f <- T

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
  library(xcms)
})

#' wl-28-08-2018, Tue: Convert a string separated by comma into character vector
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

      #' processing mzML or mzXML
      make_option("--process",
        type = "logical", default = TRUE,
        help = "Select TRUE to process mzML or mzXML files by xcms
                otherwise use the processed xset R data file"
      ),

      #' input files
      make_option("--mzxml_file",
        type = "character",
        help = "mzXML/mzML file directory or full file list separated by comma"
      ),
      make_option("--xset_file",
        type = "character",
        help = "xcmsSet R data file produced by xcms"
      ),
      make_option("--samp_name",
        type = "character", default = "",
        help = "Sample names. Default is the names of mz XML file"
      ),

      #' process raw data with mzML or mzXML format and get xcmsSet
      make_option("--FWHM",
        type = "double", default = 3,
        help = "Approximate FWHM (in seconds) of chromatographic peaks"
      ),
      make_option("--snthresh",
        type = "double", default = 5,
        help = "The signal to noise threshold"
      ),
      make_option("--profmethod",
        type = "character", default = "binlin",
        help = "Use either 'bin' (better for centroid, default), or
                'binlin' (better for profile)"
      ),
      make_option("--minfrac",
        type = "double", default = 0.25,
        help = "Minimum fraction of samples necessary for it to be
                a valid peak group"
      ),

      #' make library
      make_option("--ionisation_mode", type = "character", default = "positive"),
      make_option("--fixed", type = "logical", default = FALSE),
      make_option("--fixed_FA", type = "double", default = 16),

      #' deisotope
      make_option("--ppm", type = "integer", default = 5),
      make_option("--no_isotope", type = "integer", default = 2),
      make_option("--prop_1", type = "double", default = 0.9),
      make_option("--prop_2", type = "double", default = 0.5),

      #' annotate
      make_option("--ppm_annotate", type = "integer", default = 15),

      #' output files
      make_option("--peak_out",
        type = "character", default = "peak.tsv",
        help = "Annotated peak table"
      ),
      make_option("--rdata", type = "logical", default = TRUE),
      make_option("--rdata_out",
        type = "character", default = "xset.rdata",
        help = "xcmsSet R file after XCMS on mzML or mzXML files."
      )
    )

  opt <- parse_args(
    object = OptionParser(option_list = option_list),
    args = commandArgs(trailingOnly = TRUE)
  )
} else {
  tool_dir <- "~/R_lwc/r_data/cam1/lcms/"
  #' tool_dir <- "~/my_galaxy/lcms/"
  #' tool_dir <- "C:/R_lwc/lcms/" #' for windows
  opt <- list(
    process = F,

    #' input files
    mzxml_file = paste(paste0(tool_dir, "test-data/lcms_neg/ZH_180918_mann_neg_001.mzML"),
                       paste0(tool_dir, "test-data/lcms_neg/ZH_180918_mann_neg_002.mzML"),
                       sep = ","),

    #' mzxml_file = paste0(tool_dir, "test-data/lcms_neg"),
    xset_file = paste0(tool_dir, "test-data/res/xset_neg.rdata"),
    samp_name = "",

    #' process raw data with mzML or mzXML format and get xcmsSet
    FWHM = 3,
    snthresh = 5,
    profmethod = "binlin",
    minfrac = 0.25,

    #' make library
    ionisation_mode = "negative",
    fixed = FALSE,
    fixed_FA = 16,

    #' deisotope
    ppm = 5,
    no_isotope = 2,
    prop_1 = 0.9,
    prop_2 = 0.5,

    #' annotate
    ppm_annotate = 15,

    #' Output
    peak_out = paste0(tool_dir, "test-data/res/peak_neg.tsv"),
    rdata = TRUE,
    rdata_out = paste0(tool_dir, "test-data/res/xset_neg.rdata")
  )
}
print(opt)

suppressPackageStartupMessages({
  source(paste0(tool_dir, "lcms_func.R"))
})

## ==== Pre-processing ====

#' library files
lib_dir <- paste0(tool_dir, "libraries/")

read <- read.csv(paste(lib_dir, "lib_FA.csv", sep = "/"), sep = ",",
                 header = T)
lookup_FA <- read[, 2:4]
row.names(lookup_FA) <- read[, 1]

read <- read.csv(paste(lib_dir, "lib_class.csv", sep = "/"), sep = ",",
                 header = T)
lookup_lipid_class <- read[, 2:3]
row.names(lookup_lipid_class) <- read[, 1]

read <- read.csv(paste(lib_dir, "lib_element.csv", sep = "/"), sep = ",",
                 header = T)
lookup_element <- read[, 2:3]
row.names(lookup_element) <- read[, 1]

read <- read.csv(paste(lib_dir, "lib_modification.csv", sep = "/"),
                 sep = ",", header = T)
lookup_mod <- read[, 2:ncol(read)]
row.names(lookup_mod) <- read[, 1]

## ==== Main process ====

#' XCMS
if (opt$process) {

  #' process multiple input files separated by comma
  #' wl-04-03-2019, Mon: add file directory option. Note that it is not for
  #' galaxy.
  if (dir.exists(opt$mzxml_file)) { ## file directory
    opt$mzxml_file <- list.files(opt$mzxml_file,
      pattern = "mzml|mzxml",
      ignore.case = T, recursive = F,
      full.names = TRUE
    )
  } else { ## multiple files
    opt$mzxml_file <- str_vec(opt$mzxml_file)
  }

  #' handle sample names (especially for galaxy)
  if (opt$samp_name == "") {
    opt$samp_name <- opt$mzxml_file
  } else {
    opt$samp_name <- str_vec(opt$samp_name)
  }
  #' extract only sample names (use greedy match)
  opt$samp_name <- gsub(".*/|\\..*$", "", opt$samp_name, perl = T)

  #' Construct xcmsSet objects for peaks
  xset <- xcmsSet(opt$mzxml_file, snames = opt$samp_name,
    method = "matchedFilter", step = 0.1,
    sigma = opt$FWHM / 2.3548, snthresh = opt$snthresh,
    profmethod = opt$profmethod
  )
  #' wl-30-05-2019, Thu: findPeaks.matchedFilter is one of methods.

  xset <- group(xset) #' slotNames(xset)

  #' Corrects retention times
  xset <- retcor(xset,
    method = "obiwarp", profStep = 0.1,
    plottype = "none"
  ) # "deviation"
  #' wl-15-03-2018, Thu: Possible memory problem?
  #' wl-28-05-2019, Tue: 'profStep': step size (in m/z) to use for profile
  #' generation. Smaller steps yield more precision at the cost of greater
  #' memory usage. (from profStep-methods)
  #' wl-30-05-2019, Thu: two methods: 'peakgroups' and 'obiwarp'

  #' Group peaks from different samples together
  xset <- group(xset, bw = 5, minfrac = opt$minfrac, mzwid = 0.025)
  #' lwc-05-11-2013: three methods: group.density (default),
  #'  group.mzClust and group.nearest.

  #' Integrate areas of missing peaks
  xset <- fillPeaks(xset)
  #' Note: need mzML files to fill in missing peaks

  #' Peak lists. (peak area: value = "into"; peak height: value = "maxo")
  peaklist <- peakTable(xset,
                        method = "maxint", value = "into",
                        intensity = "maxo"
                        ) #' peak area

  #' round mz and rt, and kepp them
  peaklist <- transform(peaklist, mz = round(mz, 4), rt = round(rt, 2))
  peaklist <- peaklist[,c("mz","rt",opt$samp_name)]
  #' wl-14-10-2019, Mon: directory will be random in planemo test.
  #' peaklist <- subset(peaklist, select = -c(mzmin, mzmax, rtmin, rtmax, npeaks))
  #' peaklist <- peaklist[, -3] #' remove mzml directory name

  if (opt$rdata) {
    save(xset, peaklist, file = opt$rdata_out)
  }

} else {
    load(opt$xset_file)
}

#' deisotoping
#' wl-28-05-2019, Tue: Use average of all samples.
#' mz and intensity of the 1st sample
#' spectra <- as.matrix(peaklist[, c(1, 3)])
tmp <- round(apply(peaklist[, -c(1, 2)], 1, mean), digits = 3)
spectra <- as.matrix(cbind(mz.obs = peaklist[, 1], intensity = tmp,
                           isotope = "", modification = ""))

deisotoped <- deisotoping(ppm = opt$ppm, no_isotope = opt$no_isotope,
                          prop.1 = opt$prop_1, prop.2 = opt$prop_2,
                          spectra = spectra)

#' make library
dbase <- makelibrary(ionisation_mode = opt$ionisation_mode,
                     fixed = opt$fixed,
                     fixed_FA = opt$fixed_FA,
                     lookup_lipid_class = lookup_lipid_class,
                     lookup_FA = lookup_FA,
                     lookup_element = lookup_element)

#' annotating
annotated <- annotating(ionisation_mode = opt$ionisation_mode,
                        deisotoped = deisotoped,
                        ppm.annotate = opt$ppm_annotate,
                        dbase = dbase)

#' update and normalise peaklist
indices <- match(deisotoped[, "mz.obs"], peaklist[, "mz"])
tmp <- peaklist[indices, ]
rownames(tmp) <- NULL
final_peak <- cbind(annotated, tmp)

#' normalising based on TIC
tmp <- norm_tic(final_peak[, -c(1:3)], dim = 2)
final_peak <- cbind(final_peak[, 1:3], tmp)

#' save results
write.table(final_peak, file = opt$peak_out, sep = "\t", row.names = F)

#' write.table(final_peak, file = "peak_neg.tsv", sep = "\t", row.names = F)

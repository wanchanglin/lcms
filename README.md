# LCMS for Galaxy #

Galaxy tool for LC-MS data conversion and annotation for lipdomics analysis.


The original codes by Zoe Hall is located in 
[lcms_processing](https://github.com/hallz/lcms_processing_). 

## Progress ##

There are three R files:

- `xcms.R`: Convert mzML files using R package `xcms`. The parameters are
  chosen based on Zoe Hall's `lcms.R`. The R data set will be used for
  annotation and peak table. Also in this R script, several methods for peak
  list retrieve are given.
- `anno.R`: Annotation for lipdomics. The peak table produced includes a
  column for the results of annotation. Note that this script is slightly
  modified from `lcms.R`.
- `lcms.R`: Zoe Hall's original R code with a little tidy up.



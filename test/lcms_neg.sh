#!/bin/bash
# wl-25-05-2019, Sat: use file names for mzxml_file

Rscript --vanilla ../lcms.R \
  --process T \
  --mzxml_file "../test-data/lcms_neg/ZH_180918_mann_neg_001.mzML, \
                ../test-data/lcms_neg/ZH_180918_mann_neg_002.mzML" \
  --ionisation_mode "negative" \
  --peak_out "../test-data/res/peak_neg.tsv" \
  --rdata TRUE\
  --rdata_out "../test-data/res/xset_neg.rdata"

#!/bin/bash
# wl-25-05-2019, Sat: use file names for mzxml_file

Rscript --vanilla ../lcms.R \
  --process T \
  --mzxml_file "../test-data/lcms_neg/ZH_180918_mann_neg_001.mzML, \
                ../test-data/lcms_neg/ZH_180918_mann_neg_002.mzML, \
                ../test-data/lcms_neg/ZH_180918_mann_neg_003.mzML, \
                ../test-data/lcms_neg/ZH_180918_mann_neg_004.mzML, \
                ../test-data/lcms_neg/ZH_180918_mann_neg_005.mzML, \
                ../test-data/lcms_neg/ZH_180918_mann_neg_006.mzML, \
                ../test-data/lcms_neg/ZH_180918_mann_neg_007.mzML" \
  --ionisation_mode "negative" \
  --peak_out "../test-data/peak_neg.tsv" \
  --rdata TRUE\
  --rdata_out "../test-data/xset_neg.rdata"

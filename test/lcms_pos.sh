#!/bin/bash
# wl-05-06-2019, Wed: use file names for mzxml_file

Rscript --vanilla ../lcms.R \
  --process T \
  --mzxml_file "../test-data/lcms_pos/ZH_190918_mann_liver_pos_001.mzML, \
                ../test-data/lcms_pos/ZH_190918_mann_liver_pos_002.mzML, \
                ../test-data/lcms_pos/ZH_190918_mann_liver_pos_003.mzML, \
                ../test-data/lcms_pos/ZH_190918_mann_liver_pos_004.mzML" \
  --ionisation_mode "positive" \
  --peak_out "../test-data/res/peak_pos.tsv" \
  --rdata TRUE\
  --rdata_out "../test-data/res/xset_pos.rdata"

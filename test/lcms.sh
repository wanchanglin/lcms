#!/bin/bash
# wl-25-05-2019, Sat: use file names for mzxml_file

Rscript --vanilla ../lcms.R \
  --process T \
  --mzxml_file "../test-data/lcms_pos/ZH_190918_mann_liver_pos_001.mzML, \
                ../test-data/lcms_pos/ZH_190918_mann_liver_pos_002.mzML, \
                ../test-data/lcms_pos/ZH_190918_mann_liver_pos_003.mzML, \
                ../test-data/lcms_pos/ZH_190918_mann_liver_pos_004.mzML, \
                ../test-data/lcms_pos/ZH_190918_mann_liver_pos_005.mzML, \
                ../test-data/lcms_pos/ZH_190918_mann_liver_pos_006.mzML, \
                ../test-data/lcms_pos/ZH_190918_mann_liver_pos_007.mzML" \
  --ionisation_mode "positive" \
  --peak_out "../test-data/peak_pos_fils.tsv" \
  --rdata TRUE\
  --rdata_out "../test-data/xset_pos_file.rdata"

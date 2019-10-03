#!/bin/bash
# wl-25-05-2019, Sat: use directory for mzxml_file

Rscript --vanilla ../lcms.R \
  --process T \
  --mzxml_file "../test-data/lcms_pos" \
  --xset_file "../test-data/res/xset_pos.rdata" \
  --ionisation_mode "positive" \
  --peak_out "../test-data/res/peak_pos.tsv" \
  --rdata TRUE\
  --rdata_out "../test-data/res/xset_pos.rdata"

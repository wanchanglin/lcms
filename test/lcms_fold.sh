#!/bin/bash
# wl-25-05-2019, Sat: use directory for mzxml_file

Rscript --vanilla ../lcms.R \
  --process F \
  --mzxml_file "../test-data/lcms_pos" \
  --xset_file "../test-data/xset_pos.rdata" \
  --ionisation_mode "negative" \
  --peak_out "../test-data/peak_pos_tmp.tsv" \
  --rdata TRUE\
  --rdata_out "../test-data/xset_pos.rdata"

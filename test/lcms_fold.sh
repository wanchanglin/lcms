#!/bin/bash
# wl-25-05-2019, Sat: use directory for mzxml_file

Rscript --vanilla ../lcms.R \
  --process T \
  --mzxml_file "../test-data/lcms_neg" \
  --ionisation_mode "negative" \
  --peak_out "../test-data/peak_neg.tsv" \
  --rdata TRUE\
  --rdata_out "../test-data/xset_neg.rdata"

#!/bin/bash
# wl-01-11-2017, Wed: Rscript test code for Linux
# wl-25-03-2019, Mon: add output directory
# --xset_file "../test-data/xset_pos.rdata" \

Rscript --vanilla ../lcms.R \
  --process T \
  --mzxml_file "../test-data/lcms_neg" \
  --ionisation_mode "negative" \
  --peak_out "../test-data/peak_tmp.tsv" \
  --rdata TRUE\
  --rdata_out "../test-data/xset_tmp.rdata"

#!/bin/bash
# wl-01-11-2017, Wed: Rscript test code for Linux
# wl-25-03-2019, Mon: add output directory

Rscript --vanilla ../lcms.R \
  --process F \
  --mzxml_file "../test-data/lcms_pos" \
  --xset_file "../test-data/xset_pos.rdata" \
  --peak_out "../test-data/peak_tmp.tsv" \
  --rdata TRUE\
  --rdata_out "../test-data/xset_tmp.rdata"

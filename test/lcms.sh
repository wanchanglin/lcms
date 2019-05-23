#!/bin/bash
# wl-01-11-2017, Wed: Rscript test code for Linux
# wl-25-03-2019, Mon: add output directory

Rscript --vanilla ../lcms.R \
  --process TRUE \
  --mzxml_file "../test-data/lcms_pos" \
  --peak_out "../test-data/peak_deb.tsv" \
  --rdata TRUE\
  --rdata_out "../test-data/xset_pos_deb.rdata"

#!/bin/bash
# wl-14-10-2019, Mon: test R data

Rscript --vanilla ../lcms.R \
  --process F \
  --xset_file "../test-data/res/xset_pos.rdata" \
  --ionisation_mode "positive" \
  --peak_out "../test-data/res/peak_pos.tsv" 

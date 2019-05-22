rem wl-22-05-2019, Wed:

Rscript --vanilla ../lcms.R ^
  --process TRUE ^
  --mzxml_file "../test-data/lcms_pos" ^
  --peak_out "../test-data/peak.tsv" ^
  --rdata TRUE ^
  --rdata_out "../test-data/xset_pos_deb.rdata" ^

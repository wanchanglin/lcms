rem wl-22-05-2019, Wed:

Rscript --vanilla ../lcms.R ^
  --process F ^
  --mzxml_file "../test-data/lcms_pos" ^
  --xset_file "../test-data/xset_pos.rdata" ^
  --peak_out "../test-data/peak_tmp.tsv" ^
  --rdata TRUE ^
  --rdata_out "../test-data/xset_tmp.rdata" ^

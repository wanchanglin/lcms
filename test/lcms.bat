rem wl-22-05-2019, Wed:

Rscript --vanilla ../lcms.R ^
  --process F ^
  --mzxml_file "../test-data/lcms_pos" ^
  --xset_file "../test-data/res/xset_pos.rdata" ^
  --peak_out "../test-data/res/peak_pos.tsv" ^
  --rdata TRUE ^
  --rdata_out "../test-data/res/xset_pos.rdata" ^

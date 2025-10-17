# run example (make sure you have 1000g and hgdp in your workfolder)
python Get_data.py -outfile=MMP20 -chrom=chr11 -start=102_576_832 -end=102_625_332
Rscript plot.R MMP20 chr11 102576832 102625332
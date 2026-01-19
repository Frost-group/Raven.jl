set datafile separator "\t"
set key outside
set xlabel "S"
set ylabel "D"
set term pngcairo size 1000,650 enhanced
set output "SvsD.png"

plot for [i=0:7] "data/attempt1.tsv" index i using 2:3 with lines title sprintf("β block %d", i)
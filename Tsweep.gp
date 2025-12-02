set term pngcairo size 1200,800 font ",12"
set datafile separator whitespace
set datafile commentschars "#"
set grid
set key right top

file = "data/defect_Tsweep.tsv"

# Choose which temperatures to plot (strings that match the T column numerically)
# Edit this line to match what you ran:
Temps = "1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0 10.0"

# 1) D_tr vs defect fraction, multiple T
set output "plots/defects_Dtr_vs_phi_Tsweep.png"
set xlabel "Defect fraction (ϕ)"
set ylabel "D_{tr}"
plot for [t in Temps] file u ( (abs($1 - t) < 1e-9) ? $2 : 1/0 ):3 w lp lw 2 pt 7 title sprintf("T=%s", t)

# 2) Haven vs defect fraction
set output "plots/defects_Haven_vs_phi_Tsweep.png"
set ylabel "Haven (D_{tr}/D_{bulk})"
set yrange [0:1.3]
plot for [t in Temps] file u ( (abs($1 - t) < 1e-9) ? $2 : 1/0 ):6 w lp lw 2 pt 7 title sprintf("T=%s", t)

# 3) Normalized Dtr(ϕ)/Dtr(0) vs defect fraction
set output "plots/defects_Dtr_norm_vs_phi_Tsweep.png"
set ylabel "D_{tr}(ϕ) / D_{tr}(0)"
set yrange [0:*]
plot for [t in Temps] file u ( (abs($1 - t) < 1e-9) ? $2 : 1/0 ):10 w lp lw 2 pt 7 title sprintf("T=%s", t)

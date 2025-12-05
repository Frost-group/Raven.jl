### plot_disorder.gnuplot
# Inputs
file      = "data/defect_Ts_combined.tsv"   # your combined file path
Ts        = "0.5 1.0 2.0"                    # list of T values you ran (space-separated)
outprefix = "plots/disorder"                 # where to save images

# General styling
set datafile separator '\t'
set key outside right
set grid
set term pngcairo size 1400,900
set samples 500
set tics nomirror
set border lw 1.2
set rmargin 10

# Axes
set xlabel "Disorder strength  S = β^{2} σ_V^{2} / 3"
# Uncomment if you want log y for diffusion/cond
# set logscale y

# --- Helper macro: plot many Ts on same axes, filtering rows by column(1)==T ---
# x-col = 4 (S), y-col passed as arg
# Float comparison is OK here since you used exact decimals; if needed, use tolerance with abs($1-t)<1e-9
plotTs(ycol, label) = sprintf(\
    "plot for [t in '%s'] '%s' using (($1==t)? $4:1/0):%d with linespoints lw 2 pt 7 title sprintf('T=%%s', t)",\
    Ts, file, ycol)

# Ensure output dir exists (gnuplot won't create it; just a note)
# mkdir -p plots

# 1) Tracer diffusion vs S
set output outprefix."_Dtr.png"
set ylabel "Tracer diffusion  D_{tr}"
eval plotTs(6, "Dtr")

# 2) Collective/bulk diffusion vs S
set output outprefix."_Dbulk.png"
set ylabel "Collective diffusion  D_{bulk}"
eval plotTs(7, "Dbulk")

# 3) Haven ratio vs S
set output outprefix."_Haven.png"
unset logscale y
set ylabel "Haven ratio  H = D_{tr}/D_{bulk}"
set yrange [0:*]
eval plotTs(8, "Haven")

# 4) Tracer correlation factor f_tr vs S
set output outprefix."_ftr.png"
set ylabel "Tracer correlation factor  f_{tr}"
set yrange [0:*]
eval plotTs(9, "f_tr")

# 5) Collective correlation factor f_col vs S
set output outprefix."_fcol.png"
set ylabel "Collective correlation factor  f_{col}"
set yrange [0:*]
eval plotTs(10, "f_col")

# 6) Reduced conductivity vs S
set output outprefix."_sigma_red.png"
set ylabel "Reduced conductivity  C q^{2} D_{tr} / (k_B T)"
# set logscale y
eval plotTs(11, "sigma_red")

unset output
print sprintf("Wrote images to %s_*.png", outprefix)

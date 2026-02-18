# Set the output to a PNG file
set terminal pngcairo enhanced font 'Arial,10'
set output 'sweeps_vs_trMSD_over_sigma.png'
file = 'combined_500K_runs.tsv'

# Set the title and labels
set title '500K, N=400 with varying potential amplitudes'
set xlabel 'Sweeps'
set ylabel 'Tracer MSD'

set xrange [0:10000000]

# Set datafile separator to tab
set datafile separator '\t'

# Define the list of sigma values
sigma_values = "1.0 3.0 5.0 7.0 10.0"

n = words(sigma_values)

# Plot the data with labels for each sigma value
#plot for [i=1:5] 'combined_500K_runs.tsv' every ::(i-1)*lines_per_chunk::i*lines_per_chunk using 2:3 with lines title sprintf("σ = %s", sigma_array[i-1])

plot for [i=0:n-1] file index i using 2:3 with lines \
    title sprintf("σ=%s", word(sigma_values, i+1))
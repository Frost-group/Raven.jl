# List all your TSV files
tsv_files = ["beta_vs_sweep_msd_N400-sigma1.0.tsv", "beta_vs_sweep_msd_N400-sigma3.0.tsv", "beta_vs_sweep_msd_N400-sigma5.0.tsv", "beta_vs_sweep_msd_N400-sigma7.0.tsv", "beta_vs_sweep_msd_N400-sigma10.0.tsv"]  # Add your files here
output_file = "combined_500K_runs.tsv"

# Open the output file for writing
open(output_file, "w") do out_f
    for tsv in tsv_files
        # Read the file and split by two blank lines
        file_content = read(tsv, String)
        chunks = split(file_content, "\n\n")  # Split by two blank lines
        
        first_chunk = chunks[1]  # Get the first chunk
        # Write the chunk to the output file with two blank lines separation
        println(out_f, first_chunk)
        println(out_f, "")  # Blank line separator
        println(out_f, "")  
    end
end

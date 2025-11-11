#!/usr/bin/env julia
# run_morgan.jl
#
# Usage:
# ./run_morgan.jl a b c sweeps lagtime [a_lat] [kB] [T] [q] [outfile]
#
# Example:
# ./run_morgan.jl 30 30 30 50000 20 1.0 1.0 1.0 1.0 data/haven_sweep.tsv

using Pkg
# Activate project
isfile("Project.toml") && Pkg.activate(".")

# Load Raven
using Raven

# argparse
function argparse(T, i, default)
    (length(ARGS) ≥ i) ? parse(T, ARGS[i]) : default
end

if length(ARGS) < 5
    println("Usage: ./run_morgan.jl a b c sweeps lagtime [a_lat] [kB] [T] [q] [outfile]")
    exit(1)
end

a       = parse(Int, ARGS[1])
b       = parse(Int, ARGS[2])
c       = parse(Int, ARGS[3])
sweeps  = parse(Int, ARGS[4])
lagtime = parse(Int, ARGS[5])

a_lat   = argparse(Float64, 6, 1.0)
kB      = argparse(Float64, 7, 1.0)
T       = argparse(Float64, 8, 1.0)
q       = argparse(Float64, 9, 1.0)
outfile = (length(ARGS) ≥ 10) ? ARGS[10] : "data/Morgan_noninteracting.tsv"

# Call Sweep from src/mc.jl\
Raven.Sweep(a, b, c, sweeps, lagtime;
                  a_lat = a_lat, kB = kB, T = T, q = q, outfile = outfile)

println("Done. Wrote ", outfile)
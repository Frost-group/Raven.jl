#!/usr/bin/env julia
# run_defects.jl
#
# Usage:
#   julia --project=. run_defects.jl a b c Nions sweeps lagtime [phi_min] [phi_step] [phi_max] [outfile]
#
# Example:
#   julia --project=. run_defects.jl 30 10 10 150 50000 5 0.0 0.05 0.60 data/defect_sweep.tsv

using Pkg
# Activate project (assumes this file lives in scripts/)
isfile("Project.toml") && Pkg.activate(".")

using Raven
using Printf
using Dates

argparse(::Type{T}, i, default) where {T} = (length(ARGS) ≥ i) ? parse(T, ARGS[i]) : default
argstr(i, default) = (length(ARGS) ≥ i) ? ARGS[i] : default

if length(ARGS) < 6
    println("Usage: run_defects.jl a b c Nions sweeps lagtime [phi_min] [phi_step] [phi_max] [outfile]")
    exit(1)
end

a       = parse(Int, ARGS[1])
b       = parse(Int, ARGS[2])
c       = parse(Int, ARGS[3])
Nions   = parse(Int, ARGS[4])
sweeps  = parse(Int, ARGS[5])
lagtime = parse(Int, ARGS[6])

phi_min = argparse(Float64, 7, 0.0)
phi_step= argparse(Float64, 8, 0.05)
phi_max = argparse(Float64, 9, 0.95)
outfile = argstr(10, "defect_sweep.tsv")

mkpath(dirname(outfile))

φs = collect(phi_min:phi_step:phi_max)
@printf("Defect sweep: a=%d b=%d c=%d  Nions=%d  sweeps=%d  lag=%d  φ∈[%.2f:%.2f:%.2f]\n",
        a,b,c,Nions,sweeps,lagtime,phi_min,phi_step,phi_max)
@printf("→ writing %s\n", outfile)

# Run your function from src/mc.jl
Raven.DefectSweepFixedN(a, b, c, Nions, sweeps, lagtime;
    defect_fracs = φs,
    outfile = outfile)

println("Done at ", Dates.now())

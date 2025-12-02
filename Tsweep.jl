#!/usr/bin/env julia
# run_defects_Ts.jl
#
# Usage:
#   julia --project=. run_defects_Ts.jl a b c Nions sweeps lagtime "Tlist" [phi_min] [phi_step] [phi_max] [outprefix]
# Example:
#   julia --project=. run_defects_Ts.jl 30 10 10 150 50000 5 "0.5,1.0,2.0" 0.0 0.05 0.60 data/defect_Ts

using Pkg
isfile("Project.toml") && Pkg.activate(".")
using Raven, Printf, Dates

argparse(::Type{T}, i, default) where {T} = (length(ARGS) ≥ i) ? parse(T, ARGS[i]) : default
argstr(i, default) = (length(ARGS) ≥ i) ? ARGS[i] : default

if length(ARGS) < 7
    println("Usage: run_defects_Ts.jl a b c Nions sweeps lagtime \"T1,T2,...\" [phi_min] [phi_step] [phi_max] [outprefix]")
    exit(1)
end

a       = parse(Int, ARGS[1])
b       = parse(Int, ARGS[2])
c       = parse(Int, ARGS[3])
Nions   = parse(Int, ARGS[4])
sweeps  = parse(Int, ARGS[5])
lagtime = parse(Int, ARGS[6])
Tlist   = parse.(Float64, split(ARGS[7], ","))

phi_min = argparse(Float64, 8, 0.0)
phi_step= argparse(Float64, 9, 0.05)
phi_max = argparse(Float64, 10, 0.95)
outprefix = argstr(11, "data/defect_Ts")

mkpath(dirname(outprefix))

φs = collect(phi_min:phi_step:phi_max)
@printf("Defect sweep (fixed N) over temperatures: a=%d b=%d c=%d  Nions=%d  sweeps=%d  lag=%d\n", a,b,c,Nions,sweeps,lagtime)
@printf("φ ∈ [%.2f:%.2f:%.2f], Tlist = %s\n", phi_min, phi_step, phi_max, join(Tlist, ", "))

# Combined table across T, keyed by disorder strength (β^2 χ0 / 3)
combined_path = outprefix * "_combined.tsv"
open(combined_path, "w") do io
    println(io, "T\tdefect_frac\tM\tbeta2chi0_over3\tchi0\tDtr\tDbulk\tHaven\tf_tr\tf_col\treduced_conductivity\tDtr_norm\tDbulk_norm")
end

for T in Tlist
    outfile = @sprintf("%s_T%.3f.tsv", outprefix, T)
    @printf("→ T=%.3f  writing %s\n", T, outfile)

    # NOTE: make sure Raven.DefectSweepFixedN forwards kB and T into mcloop_g! (see patch below)
    φs, Mlist, Dtrs, Dbulks, Havens, f_tr, f_col, σred, Dtr_norm, Dbulk_norm, χ0s, Svals =
        Raven.DefectSweepFixedN(a, b, c, Nions, sweeps, lagtime;
                                defect_fracs = φs, T=T, outfile = outfile,
                                return_disorder=true)  # see patch below

    # append to combined file
    open(combined_path, "a") do io
        for k in eachindex(φs)
            @printf(io, "%.12g\t%.6f\t%d\t%.12g\t%.12g\t%.12g\t%.12g\t%.12g\t%.12g\t%.12g\t%.12g\t%.12g\t%.12g\n",
                    T, φs[k], Mlist[k], Svals[k], χ0s[k], Dtrs[k], Dbulks[k], Havens[k], f_tr[k], f_col[k],
                    σred[k], Dtr_norm[k], Dbulk_norm[k])
        end
    end
end

println("Done at ", Dates.now())
println("Combined file: ", combined_path)

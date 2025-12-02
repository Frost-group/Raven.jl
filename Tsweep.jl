#!/usr/bin/env julia
# Tsweep.jl
#
# Usage:
#   julia --project=. Tsweep.jl a b c Nions sweeps lagtime  Tmin  Tstep  Tmax  [phi_min] [phi_step] [phi_max] [outfile]
#
# Example:
#   julia --project=. Tsweep.jl 30 10 10 150 50000 5  0.5 0.5 3.0  0.0 0.05 0.60  data/defect_Tsweep.tsv

using Pkg
isfile("Project.toml") && Pkg.activate(".")

using Raven
using Printf
using Dates

argparse(::Type{T}, i, default) where {T} = (length(ARGS) ≥ i) ? parse(T, ARGS[i]) : default
argstr(i, default) = (length(ARGS) ≥ i) ? ARGS[i] : default

if length(ARGS) < 9
    println("Usage:")
    println("  julia --project=. run_defects_Tsweep.jl a b c Nions sweeps lagtime Tmin Tstep Tmax [phi_min] [phi_step] [phi_max] [outfile]")
    exit(1)
end

a       = parse(Int, ARGS[1])
b       = parse(Int, ARGS[2])
c       = parse(Int, ARGS[3])
Nions   = parse(Int, ARGS[4])
sweeps  = parse(Int, ARGS[5])
lagtime = parse(Int, ARGS[6])

Tmin  = parse(Float64, ARGS[7])
Tstep = parse(Float64, ARGS[8])
Tmax  = parse(Float64, ARGS[9])

phi_min  = argparse(Float64, 10, 0.0)
phi_step = argparse(Float64, 11, 0.05)
phi_max  = argparse(Float64, 12, 0.95)
outfile  = argstr(13, "data/defect_Tsweep.tsv")

mkpath(dirname(outfile))

Ts = collect(Tmin:Tstep:Tmax)
φs = collect(phi_min:phi_step:phi_max)

@printf("Defect sweep (multi-T): a=%d b=%d c=%d  Nions=%d  sweeps=%d  lag=%d\n", a,b,c,Nions,sweeps,lagtime)
@printf("  T ∈ [%.3g:%.3g:%.3g], φ ∈ [%.3g:%.3g:%.3g]\n", Tmin,Tstep,Tmax, phi_min,phi_step,phi_max)
@printf("→ writing %s\n", outfile)

open(outfile, "w") do io
    println(io, "# Raven defects+static field; Coulomb mcloop_g!; one row per (T, φ)")
    println(io, "T\tdefect_frac\tM\tDtr\tDbulk\tHaven\tf_tr\tf_col\treduced_conductivity\tDtr_norm\tDbulk_norm")

    for T in Ts
        @printf(">> T = %.6g\n", T)
        # call your library sweep at this temperature (now that it passes kB,T through)
        φs_out, Mlist, Dtrs, Dbulks, Havens, f_tr, f_col, σ_red, Dtr_norm, Dbulk_norm =
            Raven.DefectSweepFixedN(a, b, c, Nions, sweeps, lagtime;
                                    defect_fracs = φs,
                                    T = T,
                                    outfile = tempname())  # discard its own file; we’ll write combined

        for k in eachindex(φs_out)
            Printf.@printf(io, "%.9g\t%.9g\t%d\t%.12g\t%.12g\t%.12g\t%.12g\t%.12g\t%.12g\t%.12g\t%.12g\n",
                           T, φs_out[k], Mlist[k],
                           Dtrs[k], Dbulks[k], Havens[k], f_tr[k], f_col[k], σ_red[k],
                           Dtr_norm[k], Dbulk_norm[k])
        end
        flush(io)
    end
end

println("Done at ", Dates.now())

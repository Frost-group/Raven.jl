using Pkg
# Activate project
isfile("Project.toml") && Pkg.activate(".")

using Raven
using Printf
using Dates

a, b, c = 30, 30, 10
N, M = 1, 1
steps = 200000
dr, acc_log, _, _, _ = Raven.mcloop_g!(a, b, c, N, M, steps; A=-3.0, sigma=2.0, a_lat=1.0, kB=1.0, T=1.0)
msd_tr, lag = Raven.msd(dr, 10)
println("D_tr ~ ", Raven.tracerD(msd_tr, 3, lag))

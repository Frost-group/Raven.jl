using Pkg
# Activate project
isfile("Project.toml") && Pkg.activate(".")

using Raven
using Printf
using Dates

# create a tiny system with defects
a, b, c = 30, 30, 10
N, M = 100, 150
-,-,-,-, occ, _ = Raven.initialize(a, b, c, N, M)

V = Raven.build_potential(a, b, c, occ; A = -2.0, sigma = 2.0, a_lat = 1.0)

z0 = cld(c, 2)
Raven.write_slice_tsv(V, z0, "data/pes_slice_z$(z0).tsv")
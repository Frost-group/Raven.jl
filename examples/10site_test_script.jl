# 10site_test_script
using Raven, Plots

edges = Raven.readToFeTedges(open("testdata/10chain.edge","r"))

W = Raven.ratematrix(edges, Raven.rateAdiabatic)

P = zeros(10)
P[1] = 1.0

T, hist = Raven.MEloop(W,P)

plt = plot(hist.steps, hist.diffs, xlabel="step", ylabel="|ΔP|^2, norm between consecutive dP/dts", yscale=:log10, lw=2, legend=false, title="Master equation convergence")
plots_dir = joinpath(@__DIR__, "plots")
name = "stepsvdiffs.png"
outpath = joinpath(plots_dir, name)
savefig(plt, outpath)
println("Saved plot to: $outpath")
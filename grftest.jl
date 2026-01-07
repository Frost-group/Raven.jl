using Raven, Random, LsqFit, CSV, Statistics, DataFrames

σ_values = [0.1, 0.5, 1.0, 2.0, 4.0]
β = 1.0

results = []

for σ in σ_values
    st = Raven.initialization(20, 20, 20, 500; σ=σ, ξ=2.0, rng=Random.default_rng())
    times, msd = Raven.run!(st; β=β, sweeps=2000, sample_every=10, rng=Random.default_rng())
    
    coeffs = LsqFit.fit(LsqFit.Line, times[end-100:end], msd[end-100:end])
    
    D_t = coeffs[1] / (2 * 3)
    
    push!(results, (σ, D_t))
end

data_file = "data/disorder_vs_diffusion.csv"
CSV.write(data_file, DataFrame(σ=map(x -> x[1], results), D_t=map(x -> x[2], results)))
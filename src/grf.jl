using Random, FFTW, Statistics, Printf

# ---------------------------
#  GRF potential
# ---------------------------

@inline function kfreq(n::Int, L::Int)
    n <= (L ÷ 2) ? (2π*n/L) : (2π*(n-L)/L)
end

function grfpotential(a::Int, b::Int, c::Int, σ::Float64, ξ::Float64; rng=Random.default_rng())
    η  = randn(rng, a, b, c)
    ηk = FFTW.rfft(η)  # size: (a÷2+1, b, c)

    Pk_sqrt = Array{Float64}(undef, size(ηk)...)

    @inbounds for i in 1:size(ηk, 1)
        kx = 2π*(i-1)/a  # nonnegative frequencies only (rfft)
        for j in 1:b
            ky = kfreq(j-1, b)
            for k in 1:c
                kz = kfreq(k-1, c)
                k2 = kx*kx + ky*ky + kz*kz
                Pk_sqrt[i,j,k] = exp(-0.25 * (ξ^2) * k2)
            end
        end
    end

    Vk = ηk .* Pk_sqrt
    Vk[1,1,1] = 0.0 + 0.0im                # remove k=0 mode
    V = FFTW.irfft(Vk, a)                  # back to real space

    V .-= mean(V)
    V .*= (σ / std(V))                     # enforce target std (approx exact)

    return V
end

# disorder metrics (match your old definition)
@inline function disorder_strength(V, β; d::Int=3)
    vmean = mean(V)
    χ0 = mean((V .- vmean).^2)             # spatial variance over sites
    S  = (β*β) * χ0 / d
    return χ0, S
end

# ---------------------------
#  MC state + moves
# ---------------------------

const NBR = ((1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1))
@inline mod1p(i, L) = mod(i, L) + 1        # 1..L periodic wrap

Base.@kwdef mutable struct State
    a::Int
    b::Int
    c::Int
    pos::Matrix{Int}        # N×3
    occ::Array{Int,3}       # a×b×c, 0 empty else particle id
    disp::Matrix{Int64}     # 3×N unwrapped displacement (Int64 to avoid overflow)
    V::Array{Float64,3}     # potential
end

function initialization(a::Int, b::Int, c::Int, N::Int; σ=1.0, ξ=2.0, rng=Random.default_rng())
    S = a*b*c
    pos  = Matrix{Int}(undef, N, 3)
    occ  = fill(0, a, b, c)
    disp = zeros(Int64, 3, N)

    picks = randperm(rng, S)[1:N]          # OK for moderate sizes
    CI = CartesianIndices((a,b,c))
    @inbounds for id in 1:N
        x,y,z = Tuple(CI[picks[id]])
        pos[id,:] .= (x,y,z)
        occ[x,y,z] = id
    end

    V = (σ > 0) ? grfpotential(a,b,c,Float64(σ),Float64(ξ); rng=rng) : zeros(a,b,c)
    return State(a=a,b=b,c=c,pos=pos,occ=occ,disp=disp,V=V)
end

function attempt!(st::State, id::Int, β::Float64, rng)
    a,b,c = st.a, st.b, st.c
    pos, occ, disp, V = st.pos, st.occ, st.disp, st.V

    @inbounds begin
        x,y,z = pos[id,1], pos[id,2], pos[id,3]
        dx,dy,dz = NBR[rand(rng, 1:length(NBR))]

        x2 = mod1p(x-1 + dx, a)
        y2 = mod1p(y-1 + dy, b)
        z2 = mod1p(z-1 + dz, c)

        occ[x2,y2,z2] != 0 && return false

        ΔE = V[x2,y2,z2] - V[x,y,z]
        if (ΔE <= 0.0) || (rand(rng) < exp(-β*ΔE))
            occ[x,y,z] = 0
            occ[x2,y2,z2] = id
            pos[id,:] .= (x2,y2,z2)
            disp[1,id] += dx; disp[2,id] += dy; disp[3,id] += dz
            return true
        end
    end
    return false
end

# ---------------------------
#  MSD logging (fixed-lag, time-origin averaged) without big dr_log
# ---------------------------

mutable struct Ring
    buf::Array{Int64,3}   # (L,3,N)
    t::Int
end

Ring(N::Int, L::Int) = Ring(zeros(Int64, L, 3, N), 0)

# Push displacement snapshot; return MSD over a fixed lag if available, else NaN.
function push_and_msd!(ring::Ring, disp::Matrix{Int64}, lag_samples::Int)
    ring.t += 1
    L = size(ring.buf, 1)
    idx = (ring.t-1) % L + 1
    @inbounds ring.buf[idx, :, :] .= disp

    ring.t <= lag_samples && return NaN

    idx0 = (ring.t-1-lag_samples) % L + 1
    N = size(disp,2)

    total = 0.0
    @inbounds for i in 1:N
        dx = Float64(disp[1,i] - ring.buf[idx0,1,i])
        dy = Float64(disp[2,i] - ring.buf[idx0,2,i])
        dz = Float64(disp[3,i] - ring.buf[idx0,3,i])
        total += dx*dx + dy*dy + dz*dz
    end
    return total / N
end

function run!(st::State; β=1.0, sweeps=10000, sample_every=10, lag_sweeps=200, rng=Random.default_rng())
    N = size(st.pos, 1)

    lag_samples = max(1, cld(lag_sweeps, sample_every))
    eff_lag_sweeps = lag_samples * sample_every

    ring = Ring(N, lag_samples + 1)

    times = Int[]
    msd0  = Float64[]        # MSD from origin (cheap)
    msdτ  = Float64[]        # time-origin averaged MSD at fixed lag τ

    total_attempts = 0
    total_accepts  = 0

    @inbounds for s in 1:sweeps
        for _ in 1:N
            id = rand(rng, 1:N)
            total_attempts += 1
            total_accepts  += attempt!(st, id, β, rng) ? 1 : 0
        end

        if s % sample_every == 0
            push!(times, s)

            # origin-based MSD: <|disp|^2>
            acc = 0.0
            for i in 1:N
                dx = Float64(st.disp[1,i]); dy = Float64(st.disp[2,i]); dz = Float64(st.disp[3,i])
                acc += dx*dx + dy*dy + dz*dz
            end
            push!(msd0, acc / N)

            # fixed-lag, time-origin averaged MSD
            push!(msdτ, push_and_msd!(ring, st.disp, lag_samples))
        end
    end

    acc_ratio = total_accepts / max(1, total_attempts)
    return (times=times, msd0=msd0, msdτ=msdτ, lag=eff_lag_sweeps, acc_ratio=acc_ratio)
end

# diffusion estimate from fixed-lag msdτ series (ignore early NaNs)
function D_from_msdlag(msdτ::Vector{Float64}, lag_sweeps::Int; d::Int=3)
    vals = filter(isfinite, msdτ)
    isempty(vals) && return NaN
    return mean(vals) / (2 * d * lag_sweeps)
end

# ---------------------------
#  Scan + TSV logging: disorder vs diffusion
# ---------------------------

function scan_disorder(outfile::AbstractString;
    σ_values = [0.1, 0.5, 1.0, 2.0, 4.0],
    ξ::Float64 = 2.0,
    β::Float64 = 1.0,
    a::Int = 20, b::Int = 20, c::Int = 20,
    N::Int = 500,
    sweeps::Int = 20000,
    sample_every::Int = 10,
    lag_sweeps::Int = 200,
    seed::Int = 1
)
    mkpath(dirname(outfile))

    rng = MersenneTwister(seed)

    open(outfile, "w") do io
        println(io, "sigma\txi\tbeta\tchi0\tS\tD\tacc_ratio\tlag_sweeps\ta\tb\tc\tN\tsweeps\tsample_every\tseed")
        for σ in σ_values
            # use a fresh sub-RNG per σ so each run is reproducible but different
            r = MersenneTwister(rand(rng, UInt))

            st = initialization(a,b,c,N; σ=σ, ξ=ξ, rng=r)

            χ0, S = disorder_strength(st.V, β; d=3)

            out = run!(st; β=β, sweeps=sweeps, sample_every=sample_every, lag_sweeps=lag_sweeps, rng=r)
            D = D_from_msdlag(out.msdτ, out.lag; d=3)

            @printf(io, "%.6g\t%.6g\t%.6g\t%.12g\t%.12g\t%.12g\t%.6g\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
                    σ, ξ, β, χ0, S, D, out.acc_ratio, out.lag, a, b, c, N, sweeps, sample_every, seed)

            @printf("σ=%.3g  χ0=%.3g  S=%.3g  D=%.3g  acc=%.3f  (lag=%d sweeps)\n",
                    σ, χ0, S, D, out.acc_ratio, out.lag)
        end
    end

    @printf("Wrote %s\n", outfile)
    return nothing
end

function production_run(outfile; S_max = 20, ξ = 2.0, βs = [0.0005, 0.001, 0.025, 0.05, 0.1, 0.25, 0.5, 1.0],
    a = 20, b = 20, c = 20,
    N = 500,
    sweeps = 100000,
    sample_every = 10,
    lag_sweeps = 200,
    seed = 1
)
    mkpath(dirname(outfile))

    S_values = collect(0:0.1:S_max)

    for β in βs
        σs = sqrt(3S_values)/βs[β]

        system = initialization(a, b, c, N; σ=σs[β], ξ=ξ, rng=seed)

        out = run!(system; β=β, sweeps=sweeps, sample_every=sample_every, lag_sweeps=lag_sweeps, rng=seed)

        D = D_from_msdlag(out.msdτ, out.lag; d=3)

        #logging here!
end

function scan_disorder2(outfile::AbstractString;
    σ_values = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 1.0, 2.0, 4.0],
    β_values = [0.5, 1.0, 2.0, 4.0],
    ξ::Float64 = 2.0,
    a::Int = 20, b::Int = 20, c::Int = 20,
    N::Int = 500,
    sweeps::Int = 100_000,
    sample_every::Int = 100,
    lag_sweeps::Int = 200,
    seed::Int = 42
)
    mkpath(dirname(outfile))

    σs = sort(unique(Float64.(σ_values)))
    βs = Float64.(β_values)

    # master RNG for reproducibility
    rng_master = MersenneTwister(seed)

    # IMPORTANT: same disorder realisation per σ across β
    σ_seed = Dict{Float64,UInt}()
    for σ in σs
        σ_seed[σ] = rand(rng_master, UInt)
    end

    open(outfile, "w") do io
        println(io, "sigma\txi\tbeta\tT\tchi0\tS\tu\tD\tD_over_D0\tacc_ratio\tlag_sweeps\ta\tb\tc\tN\tsweeps\tsample_every\tseed")

        for β in βs
            T = 1.0/β
            D0 = NaN

            for σ in σs
                # deterministic init per σ (same V + same initial positions for all β)
                rng_init = MersenneTwister(σ_seed[σ])
                st = initialization(a,b,c,N; σ=σ, ξ=ξ, rng=rng_init)

                χ0, S = disorder_strength(st.V, β; d=3)
                u = β * sqrt(χ0)   # common “disorder parameter” used in this literature

                # MC rng (separate stream so init doesn't affect dynamics)
                rng_mc = MersenneTwister(rand(rng_master, UInt))
                out = run!(st; β=β, sweeps=sweeps, sample_every=sample_every, lag_sweeps=lag_sweeps, rng=rng_mc)
                D = D_from_msdlag(out.msdτ, out.lag; d=3)

                if σ == 0.0
                    D0 = D
                end
                DoverD0 = (isfinite(D0) && D0 != 0.0) ? D / D0 : NaN

                @printf(io, "%.6g\t%.6g\t%.6g\t%.6g\t%.12g\t%.12g\t%.12g\t%.12g\t%.12g\t%.6g\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
                        σ, ξ, β, T, χ0, S, u, D, DoverD0, out.acc_ratio, out.lag, a,b,c,N,sweeps,sample_every,seed)

                @printf("β=%.3g (T=%.4g)  σ=%.3g  S=%.3g  u=%.3g  D=%.3g  D/D0=%.3g  acc=%.3f\n",
                        β, T, σ, S, u, D, DoverD0, out.acc_ratio)
            end


            # blank line between β blocks (gnuplot likes this)
            println(io)
        end
    end

    @printf("Wrote %s\n", outfile)
    return nothing
end


# Example:
# scan_disorder("data/grf_disorder_vs_diffusion.tsv"; σ_values=[0.1,0.5,1,2,4], ξ=2.0, β=1.0, a=20,b=20,c=20,N=500,sweeps=20000,sample_every=10,lag_sweeps=500,seed=42)


function scan_disorder_ugrid(outfile::AbstractString;
    u_values = [0.0, 0.05, 0.1, 0.2, 0.4, 0.8, 1.6, 3.2, 6.4, 12.8, 25.6],
    β_values = [0.0005, 0.001, 0.025, 0.05, 0.1, 0.25, 0.5, 1.0],
    ξ::Float64 = 2.0,
    a::Int = 20, b::Int = 20, c::Int = 20,
    N::Int = 100,
    sweeps::Int = 100_000,
    sample_every::Int = 100,
    lag_sweeps::Int = 200,
    seed::Int = 42
)
    mkpath(dirname(outfile))

    us = Float64.(u_values)
    βs = Float64.(β_values)

    rng_master = MersenneTwister(seed)

    # Same disorder realisation per u-index across β (seeded by point index)
    u_seed = [rand(rng_master, UInt) for _ in 1:length(us)]

    open(outfile, "w") do io
        println(io, "sigma\txi\tbeta\tT\tchi0\tS\tu\tD\tD_over_D0\tacc_ratio\tlag_sweeps\ta\tb\tc\tN\tsweeps\tsample_every\tseed")

        for β in βs
            T = 1.0 / β
            D0 = NaN

            for (i, u_target) in pairs(us)
                σ = (u_target == 0.0) ? 0.0 : u_target / β

                # deterministic init per u-point (same GRF realisation across β)
                rng_init = MersenneTwister(u_seed[i])
                st = initialization(a,b,c,N; σ=σ, ξ=ξ, rng=rng_init)

                χ0, S = disorder_strength(st.V, β; d=3)
                u = β * sqrt(χ0)   # measured u (≈ u_target, small finite-size noise)

                rng_mc = MersenneTwister(rand(rng_master, UInt))
                out = run!(st; β=β, sweeps=sweeps, sample_every=sample_every, lag_sweeps=lag_sweeps, rng=rng_mc)
                D = D_from_msdlag(out.msdτ, out.lag; d=3)

                if i == 1 && u_target == 0.0
                    D0 = D
                end
                DoverD0 = (isfinite(D0) && D0 != 0.0) ? D / D0 : NaN

                @printf(io, "%.6g\t%.6g\t%.6g\t%.6g\t%.12g\t%.12g\t%.12g\t%.12g\t%.12g\t%.6g\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
                        σ, ξ, β, T, χ0, S, u, D, DoverD0, out.acc_ratio, out.lag, a,b,c,N,sweeps,sample_every,seed)

                @printf("β=%.3g (T=%.4g)  u_target=%.3g  σ=%.3g  S=%.3g  u=%.3g  D=%.3g  D/D0=%.3g  acc=%.3f\n",
                        β, T, u_target, σ, S, u, D, DoverD0, out.acc_ratio)
            end

            println(io) # blank line between β blocks
        end
    end

    @printf("Wrote %s\n", outfile)
    return nothing
end

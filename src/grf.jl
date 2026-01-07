import LinearAlgebra, IterativeSolvers, Random, Printf, FFTW, Statistics

function kfreq(n::Int, L::Int)
    n <= L ÷ 2 ? (2π*n/L) : (2π*(n-L)/L)
end

function grfpotential(a,b,c,σ,ξ::Float64=2.0,rng=Random.default_rng())
    η = randn(rng, a, b, c)
    ηk = FFTW.rfft(η)

    Pk_sqrt = Array{Float64}(undef, size(ηk)...)

    for i in 1:size(ηk, 1)
        kx = 2π*(i-1)/a
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

    Vk[1, 1, 1] = 0.0 + 0.0im

    V = FFTW.irfft(Vk, a)

    V .-= Statistics.mean(V)
    V .*= (σ / Statistics.std(V))

    return V
end

function initialization(a,b,c,N; σ=1.0, ξ=2.0, rng=Random.default_rng())
    S = a*b*c
    pos = Matrix{Int}(undef, N, 3)
    occ = fill(0, a, b, c)
    disp = fill(0, 3, N)

    picks = Random.randperm(rng, S)[1:N]
    CI = CartesianIndices((a, b, c))
    for id in 1:N
        x, y, z = Tuple(CI[picks[id]])
        pos[id, :] .= (x, y, z)
        occ[x, y, z] = id
    end

    V = (σ > 0) ? grfpotential(a,b,c,σ,ξ,rng) : zeros(a,b,c)
    return (a=a, b=b, c=c, pos=pos, occ=occ, disp=disp, V=V)
end

const NBR = ((1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1))
@inline mod1p(i, L) = mod(i, L) + 1 # 1..L periodic wrap

function attempt!(st, id, β, rng)
    a,b,c = st.a, st.b, st.c
    pos, occ, disp, V = st.pos, st.occ, st.disp, st.V

    x,y,z = pos[id,1], pos[id,2], pos[id,3]
    dx,dy,dz = NBR[rand(rng, 1:length(NBR))]

    x2 = mod1p(x-1 + dx, a)
    y2 = mod1p(y-1 + dy, b)
    z2 = mod1p(z-1 + dz, c)

    occ[x2, y2, z2] != 0 && return false

    ΔE = V[x2,y2,z2] - V[x,y,z]

    if ΔE <= 0 || rand(rng) < exp(-β*ΔE)
        occ[x,y,z] = 0
        occ[x2,y2,z2] = id
        pos[id,:] .= (x2,y2,z2)

        # displacement update
        disp[1,id] += dx; disp[2,id] += dy; disp[3,id] += dz
        return true
    end
    return false
end

function run!(st; β=1.0, sweeps=10000, sample_every=10, rng=Random.default_rng())
    N = size(st.pos, 1)
    times = Int[]
    msd_t = Float64[]

    for s in 1:sweeps
        for _ in 1:N
            id = rand(rng, 1:N)
            attempt!(st, id, β, rng)
        end

        if s % sample_every == 0
            push!(times, s)
            dx2 = @views st.disp[1,:].^2
            dy2 = @views st.disp[2,:].^2
            dz2 = @views st.disp[3,:].^2
            push!(msd_t, Statistics.mean(dx2 .+ dy2 .+ dz2))
        end
    end

    return times, msd_t
end
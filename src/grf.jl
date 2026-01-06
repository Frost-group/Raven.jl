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
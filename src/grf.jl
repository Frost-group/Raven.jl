import LinearAlgebra, IterativeSolvers, Random, Printf, FFTW

η = randn(Nx, Ny)
ηk = rfft(η)

Sk = exp(-(k^2)*ξ^2/2)
Vk = ηk .* sqrt.(Sk)

V  = irfft(Vk, Nx)          # back to real-space field (Nx, Ny)
V .-= mean(V)               # optional: zero-mean
V ./= std(V)                # optional: unit variance; then scale by σ
V .*= σ

function grfpotential(a,b,c,σ,ξ::Float64=2.0,rng=Random.default_rng())
    η = randn(rng, a, b, c)
    ηk = rfft(η)

    Pk_sqrt = Array{Float64}(undef, size(ηk)...)

    for i in 1:size(ηk, 1)
        kx = 2π*(i-1)/a
        for j in 1:b
            ky = kfreq(k-1, c)
            for k in 1:c
                kz = kfreq(k-1, c)
                k2 = kx*kx + ky*ky + kz*kz
                Pk_sqrt[i,j,k] = exp(-0.25 * (ξ^2) * k2)
            end
        end
    end

    Vk = ηk .* Pk_sqrt

    Vk[1, 1, 1] = 0.0 + 0.0im

    V = irfft(Vk, a)

    V .-= mean(V)
    V .*= (σ / std(V))

    return V
end
import LinearAlgebra, IterativeSolvers, Random, Printf, Statistics

# Note: this code uses multiple dispatch quite often so it is advised to read the comments carefully if you want to follow what's going on!

@inline function pbcΔ(d, L)
    x = mod(d, L)
    return (x <= L>>1) ? x : x - L
end

@inline function r2_pbc(x, y, z, xk, yk, zk, Lx, Ly, Lz, a_lat)
    dx = pbcΔ(x - xk, Lx); dy = pbcΔ(y - yk, Ly); dz = pbcΔ(z - zk, Lz)
    return (dx*a_lat)^2 + (dy*a_lat)^2 + (dz*a_lat)^2
end

function defect_centers(occ)
    Lx, Ly, Lz = size(occ)
    centers = Vector{NTuple{3, Int}}()
    @inbounds for x in 1:Lx, y in 1:Ly, z in 1:Lz
        if occ[x, y, z] == -1
            push!(centers, (x, y, z))
        end
    end
    return centers
end

function build_potential(a, b, c, occ; A, sigma, a_lat)
    centers = defect_centers(occ)
    V = zeros(Float64, a, b, c)
    isempty(centers) && return V
    inv2σ2 = 1.0/(2*sigma^2)
    @inbounds for x in 1:a, y in 1:b, z in 1:c
        v = 0.0
        for (xk, yk, zk) in centers
            r2 = r2_pbc(x, y, z, xk, yk, zk, a, b, c, a_lat)
            v += A * exp(-r2 * inv2σ2)
        end
        V[x, y, z] = v
    end
    return V
end

function build_coulomb_potential(a,b,c, occ; q_def = +1.0, k_e = 1.0, a_lat=1.0, eps2 = 0.25)
    centers = defect_centers(occ)
    V = zeros(Float64, a,b,c)
    isempty(centers) && return V
    @inbounds for x in 1:a, y in 1:b, z in 1:c
        v = 0.0
        for (xk, yk, zk) in centers
            r2 = r2_pbc(x,y,z,xk,yk,zk, a,b,c, a_lat)
            v += k_e * q_def / sqrt(r2 + eps2)
        end
        V[x, y, z] = v
    end
    return V
end

function disorder_strength(V; kB=1.0, T=1.0, d=3)
    β = 1.0 / (kB * T)
    vmean = Statistics.mean(V)
    χ0 = Statistics.mean((V .- vmean).^2)
    S = β^2 * χ0 / d        # spatial variance over lattice sites
    return χ0, S            # d=3, Witkoskie-Yang-Cao (2002) β^2 * χ0 / 3
end

function write_slice_tsv(V, z0, path)
    a, b, c = size(V); @assert 1 <= z0 <= c
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "x\ty\tV")
        @inbounds for y in 1:b, x in 1:a
            Printf.@printf(io, "%d\t%d\t%.9g\n", x, y, V[x, y, z0])
        end
    end
end

function initialize(a, b, c, N)
    S = a * b * c # lattice size
    pos = Matrix{Int}(undef, N, 3) # N rows, columns = (x, y, z)
    occ = fill(0, a, b, c)
    disp = fill(0, 3, N) # relative position tracker per sweep, also later used for MSD computation.
    picks = Random.randperm(S)[1:N] # shuffle matrix S to select the first N entries to determine where ions reside in.
    CI = CartesianIndices((a, b, c))
    for id in 1:N
        x, y, z = Tuple(CI[picks[id]])
        pos[id, :] .= (x, y, z)
        occ[x, y, z] = id
    end
    return a, b, c, pos, occ, disp
end

function initialize(a, b, c, N, M) # Extra variable M for the number of defects introduced in the system.
    S = a * b * c
    pos = Matrix{Int}(undef, N, 3)
    occ = fill(0, a, b, c)
    disp = fill(0, 3, N)
    picks = Random.randperm(S)[1:N + M]
    CI = CartesianIndices((a, b, c))
    for defect in 1:M
        x, y, z = Tuple(CI[picks[defect]])
        # pos[defect, :] .= (x, y, z) we don't make a defect list for now.
        occ[x, y, z] = -1
    end
    for id in 1:N
        x, y, z = Tuple(CI[picks[M + id]])
        pos[id, :] .= (x, y, z)
        occ[x, y, z] = id
    end
    return a, b, c, pos, occ, disp
end

# Alright let's make the mcstep! function do explicit arithemtic and integrate it with the neighbor function.
function mcstep!(a, b, c, pos, occ, disp)
    @inbounds begin
        ion = rand(1:size(pos, 1))
        x = pos[ion, 1]; y = pos[ion, 2]; z = pos[ion, 3]

        # pick a direction and compute neighbor without building a 6x3 matrix
        d = rand(1:6)
        if d == 1
            nx = (x == a ? 1 : x+1); ny = y; nz = z; dx = 1; dy = 0; dz = 0
        elseif d == 2
            nx = (x == 1 ? a : x-1); ny = y; nz = z; dx =-1; dy = 0; dz = 0
        elseif d == 3
            ny = (y == b ? 1 : y+1); nx = x; nz = z; dx = 0; dy = 1; dz = 0
        elseif d == 4
            ny = (y == 1 ? b : y-1); nx = x; nz = z; dx = 0; dy =-1; dz = 0
        elseif d == 5
            nz = (z == c ? 1 : z+1); nx = x; ny = y; dx = 0; dy = 0; dz = 1
        elseif d == 6
            nz = (z == 1 ? c : z-1); nx = x; ny = y; dx = 0; dy = 0; dz =-1
        end

        if occ[nx, ny, nz] == 0
            occ[x, y, z] = 0
            occ[nx, ny, nz] = ion
            pos[ion, 1] = nx; pos[ion, 2] = ny; pos[ion, 3] = nz
            disp[1, ion] += dx; disp[2, ion] += dy; disp[3, ion] += dz
            return ion
        else
            return 0
        end
    end
end

function mcstep!(a, b, c, pos, occ, disp, V, β)
    @inbounds begin
        ion = rand(1:size(pos, 1))
        x = pos[ion, 1]; y = pos[ion, 2]; z = pos[ion, 3]
        d = rand(1:6)
        if d == 1
            nx = (x == a ? 1 : x+1); ny = y; nz = z; dx = 1; dy = 0; dz = 0
        elseif d == 2
            nx = (x == 1 ? a : x-1); ny = y; nz = z; dx =-1; dy = 0; dz = 0
        elseif d == 3
            ny = (y == b ? 1 : y+1); nx = x; nz = z; dx = 0; dy = 1; dz = 0
        elseif d == 4
            ny = (y == 1 ? b : y-1); nx = x; nz = z; dx = 0; dy =-1; dz = 0
        elseif d == 5
            nz = (z == c ? 1 : z+1); nx = x; ny = y; dx = 0; dy = 0; dz = 1
        elseif d == 6
            nz = (z == 1 ? c : z-1); nx = x; ny = y; dx = 0; dy = 0; dz =-1
        end

        if occ[nx, ny, nz] != 0
            return 0
        end

        ΔU = V[nx, ny, nz] - V[x, y, z] # energy change for the hop
        if (ΔU <= 0.0) || (rand() < exp(-β*ΔU))
            occ[x, y, z] = 0
            occ[nx, ny, nz] = ion
            pos[ion, 1] = nx; pos[ion, 2] = ny; pos[ion, 3] = nz
            disp[1, ion] += dx; disp[2, ion] += dy; disp[3, ion] += dz
            return ion
        else
            return 0
        end
    end
end

function mcloop!(a, b, c, N, steps)
    a, b, c, pos, occ, disp = initialize(a, b, c, N)
    attempts = N
    sweeps = cld(steps, attempts)
    dr_log = zeros(Int32, sweeps + 1, 3, N)
    acc = zeros(Int32, N) # cumulative acceptance ratio
    acc_log = zeros(Int32, sweeps + 1, N) # log for acceptance
    dr_log[1, :, :] .= 0
    acc_log[1, :] .= 0

    total_attempts = 0
    total_accepts = 0

    for sweep in 1:sweeps
        for attempt in 1:attempts
            total_attempts += 1
            moved_ion = mcstep!(a, b, c, pos, occ, disp)
            if moved_ion != 0
                total_accepts += 1
                acc[moved_ion] += 1
            #println("attempt $(attempt + attempts * (sweep - 1))")
            end
        end
        dr_log[sweep + 1, :, :] .= disp
        acc_log[sweep + 1, :] .= acc
        #println("sweep $sweep / $sweeps")
    end
    return dr_log, acc_log, pos, occ, (total_accepts, total_attempts)
end

function mcloop!(a, b, c, N, M, steps)
    a, b, c, pos, occ, disp = initialize(a, b, c, N, M)
    attempts = N
    sweeps = cld(steps, attempts)
    dr_log = zeros(Int32, sweeps + 1, 3, N)
    acc = zeros(Int32, N) # cumulative acceptance ratio
    acc_log = zeros(Int32, sweeps + 1, N) # log for acceptance
    dr_log[1, :, :] .= 0
    acc_log[1, :] .= 0

    total_attempts = 0
    total_accepts = 0

    for sweep in 1:sweeps
        for attempt in 1:attempts
            total_attempts += 1
            moved_ion = mcstep!(a, b, c, pos, occ, disp)
            if moved_ion != 0
                total_accepts += 1
                acc[moved_ion] += 1
            #println("attempt $(attempt + attempts * (sweep - 1))")
            end
        end
        dr_log[sweep + 1, :, :] .= disp
        acc_log[sweep + 1, :] .= acc
        #println("sweep $sweep / $sweeps")
    end
    return dr_log, acc_log, pos, occ, (total_accepts, total_attempts)
end

function mcloop_g!(a, b, c, N, M, steps; q_def=+1.0, k_e=1.0, a_lat=1.0, eps2=0.25, kB=1.0, T=1.0)
    β = 1.0/(kB*T)
    a, b, c, pos, occ, disp = initialize(a, b, c, N, M)
    V = build_coulomb_potential(a,b,c,occ; q_def=q_def, k_e=k_e, a_lat=a_lat, eps2=eps2)

    attempts = N
    sweeps   = cld(steps, attempts)
    dr_log   = zeros(Int32, sweeps+1, 3, N)
    acc      = zeros(Int32, N)
    acc_log  =
     zeros(Int32, sweeps+1, N)
    dr_log[1, :, :] .= 0; acc_log[1, :] .= 0

    total_attempts=0; total_accepts=0
    for sweep in 1:sweeps
        for attempt in 1:attempts
            total_attempts += 1
            moved = mcstep!(a, b, c, pos, occ, disp, V, β)
            if moved != 0
                total_accepts += 1
                acc[moved] += 1
            end
        end
        dr_log[sweep+1, :, :] .= disp
        acc_log[sweep+1, :] .= acc
    end
    return dr_log, acc_log, pos, occ, (total_accepts, total_attempts)
end


function msd(dr, lag) # faster implementation of MSD calculation not relying on the norm function: no copying, faster computation.
    N, _, ioncount = size(dr)
    n = lag
    total = 0.0
    count = 0
    for i in 1:N-n, j in 1:ioncount
        dx = dr[i+n, 1, j] - dr[i, 1, j]
        dy = dr[i+n, 2, j] - dr[i, 2, j]
        dz = dr[i+n, 3, j] - dr[i, 3, j]
        total += dx*dx + dy*dy + dz*dz
        count += 1
    end
    return total / count, n
end

function tracerD(msd, dim, lag)
    Dtr = msd / (2 * dim * lag)
    return Dtr
end

function bulkmsd(dr, lag)
    N, _, ioncount = size(dr)
    n = lag
    total = 0.0
    count = 0
    for i in 1:(N-n)
        sx = 0.0; sy = 0.0; sz = 0.0
        for j in 1:ioncount
            sx += dr[i+n, 1, j] - dr[i, 1, j]
            sy += dr[i+n, 2, j] - dr[i, 2, j]
            sz += dr[i+n, 3, j] - dr[i, 3, j]
        end
        total += sx*sx + sy*sy + sz*sz
        count += 1
    end
    return total / count, ioncount
end

function bulkD(msd, dim, lag, N)
    Dbulk = msd / (2 * dim * lag * N)
    return Dbulk
end

function DtrSweep(a, b, c, sweeps, outfile="Dtr_sweep.tsv")
    N = a * b * c
    percentiles = 0.0 : 0.05 : 1.0
    ioncount = round.(Int, N .* percentiles)
    ioncount[1] = 1
    Dtrs = Vector{Float64}(undef, length(percentiles))
    for i in 1:length(ioncount)
        steps = sweeps * ioncount[i]
        dr, _, _ = mcloop!(a, b, c, ioncount[i], steps)
        msd, lag = Raven.msd(dr, 20)
        Dtrs[i] = tracerD(msd, 3, lag)
        Printf.@printf("pct=%.2f (N=%d)  Dtr=%g\n", percentiles[i], ioncount[i], Dtrs[i])
    end

    open(outfile, "w") do io
        println(io, "percentile\tDtr")
        for i in 1:length(ioncount)
            Printf.@printf(io, "%.2f\t%.12g\n", percentiles[i], Dtrs[i])
        end
    end
    Printf.@printf("Wrote %s\n", outfile)

    return percentiles, Dtrs
end

function DbulkSweep(a, b, c, sweeps, outfile="Dbulk_sweep.tsv")
    N = a * b * c
    percentiles = 0.0 : 0.05 : 1.0
    ioncount = round.(Int, N .* percentiles)
    ioncount[1] = 1
    Dbulks = Vector{Float64}(undef, length(percentiles))
    for i in 1:length(ioncount)
        steps = sweeps * ioncount[i]
        dr, _, _ = mcloop!(a, b, c, ioncount[i], steps)
        msd, N = Raven.bulkmsd(dr, 20)
        Dbulks[i] = bulkD(msd, 3, 20, N)
        Printf.@printf("pct=%.2f (N=%d)  Dbulk=%g\n", percentiles[i], ioncount[i], Dbulks[i])
    end

    open(outfile, "w") do io
        println(io, "percentile\tDbulk")
        for i in 1:length(ioncount)
            Printf.@printf(io, "%.2f\t%.12g\n", percentiles[i], Dbulks[i])
        end
    end
    Printf.@printf("Wrote %s\n", outfile)

    return percentiles, Dbulks
end

function Sweep(a, b, c, sweeps, lagtime;
                         a_lat::Float64=1.0, kB::Float64=1.0, T::Float64=1.0, q::Float64=1.0,
                         outfile::AbstractString="Morgan_noninteracting.tsv")

    SITES = a * b * c
    percentiles = collect(0.0 : 0.05 : 1.0)
    ioncount = round.(Int, SITES .* percentiles); ioncount[1] = 1

    Dtrs    = Vector{Float64}(undef, length(percentiles))
    Dbulks  = Vector{Float64}(undef, length(percentiles))
    Havens  = Vector{Float64}(undef, length(percentiles))
    f_tr    = Vector{Float64}(undef, length(percentiles))
    f_col   = Vector{Float64}(undef, length(percentiles))
    sigma_reds = Vector{Float64}(undef, length(percentiles))

    for k in eachindex(ioncount)
        Nions = ioncount[k]
        steps = sweeps * Nions                   # same # sweeps at each loading
        dr, acc_log, _, _, _ = mcloop!(a, b, c, Nions, steps)   # dr::(S,3,N), acc_log::(S,N)

        S, D, N = size(dr); @assert D == 3
        τ = lagtime; @assert 1 <= τ <= S-1

        # --- tracer diffusion via your helpers ---
        msd_tr, lag = msd(dr, τ)
        Dtr = tracerD(msd_tr, 3, lag)

        # --- bulk (Einstein–Helfand) via your helpers ---
        msd_col, Nin = bulkmsd(dr, τ)           # Nin should equal Nions
        Dbulk = bulkD(msd_col, 3, τ, Nin)

        H = Dtr / Dbulk

        # --- correlation factors (ratio-of-sums over windows) ---
        # tracer:   f_tr  = (Σ_s,i ||Δr_i||^2) / (Σ_s,i Δh_i * a^2)
        # collective: f_col = (Σ_s ||Σ_i Δr_i||^2) / (Σ_s Σ_i Δh_i * a^2)
        sum_d2_tr  = 0.0; sum_dh_tr  = 0.0
        sum_d2_col = 0.0; sum_dh_col = 0.0

        @inbounds for s in 1:(S-τ)
            # tracer pieces
            for i in 1:N
                dh = acc_log[s+τ,i] - acc_log[s,i]
                if dh > 0
                    dx = dr[s+τ,1,i] - dr[s,1,i]
                    dy = dr[s+τ,2,i] - dr[s,2,i]
                    dz = dr[s+τ,3,i] - dr[s,3,i]
                    sum_d2_tr += dx*dx + dy*dy + dz*dz
                    sum_dh_tr += dh
                end
            end
            # collective pieces
            sx=0.0; sy=0.0; sz=0.0; dh_sum=0.0
            for i in 1:N
                sx += dr[s+τ,1,i] - dr[s,1,i]
                sy += dr[s+τ,2,i] - dr[s,2,i]
                sz += dr[s+τ,3,i] - dr[s,3,i]
                dh_sum += (acc_log[s+τ,i] - acc_log[s,i])
            end
            if dh_sum > 0
                sum_d2_col += sx*sx + sy*sy + sz*sz
                sum_dh_col += dh_sum
            end
        end

        ftr  = (sum_dh_tr  > 0) ? sum_d2_tr  / (sum_dh_tr  * a_lat*a_lat) : NaN
        fcol = (sum_dh_col > 0) ? sum_d2_col / (sum_dh_col * a_lat*a_lat) : NaN

        # --- reduced conductivity: C q^2 Dtr / (kB T) ---
        C = Nions / SITES
        sigma_red = C * (q*q) * Dtr / (kB * T)

        Dtrs[k]   = Dtr
        Dbulks[k] = Dbulk
        Havens[k] = H
        f_tr[k]   = ftr
        f_col[k]  = fcol
        sigma_reds[k] = sigma_red

        Printf.@printf("pct=%.2f (N=%d)  Dtr=%.6g  Dbulk=%.6g  Haven=%.6g  f_tr=%.6g  f_col=%.6g  σ_red=%.6g\n",
                       percentiles[k], Nions, Dtr, Dbulk, H, ftr, fcol, sigma_red)
    end

    open(outfile, "w") do io
        println(io, "percentile\tDtr\tDbulk\tHaven\tf_tr\tf_col\treduced_conductivity")
        for k in eachindex(percentiles)
            Printf.@printf(io, "%.2f\t%.12g\t%.12g\t%.12g\t%.12g\t%.12g\t%.12g\n",
                           percentiles[k], Dtrs[k], Dbulks[k], Havens[k], f_tr[k], f_col[k], sigma_reds[k])
        end
    end
    Printf.@printf("Wrote %s\n", outfile)

    return percentiles, Dtrs, Dbulks, Havens, f_tr, f_col, sigma_reds
end


function Sweep(a, b, c, defects, sweeps, lagtime;
                        a_lat::Float64=1.0, kB::Float64=1.0, T::Float64=1.0, q::Float64=1.0,
                        outfile::AbstractString="defect_sweep.tsv")
        
    SITES = a * b * c
    percentiles = collect(0.0 : 0.05 : 1.0)
    ioncount = round.(Int, SITES .* percentiles); ioncount[1] = 1
    Dtrs    = Vector{Float64}(undef, length(percentiles))
    Dbulks  = Vector{Float64}(undef, length(percentiles))
    Havens  = Vector{Float64}(undef, length(percentiles))
    f_tr    = Vector{Float64}(undef, length(percentiles))    
    f_col   = Vector{Float64}(undef, length(percentiles))
    sigma_reds = Vector{Float64}(undef, length(percentiles))

    for k in eachindex(ioncount)
        Nions = ioncount[k]
        steps = sweeps * Nions
        dr, acc_log, _, _, _ = mcloop!(a, b, c, Nions, defects, steps)

        S, D, N = size(dr); @assert D == 3
        τ = lagtime; @assert 1 <= τ <= S-1
        # tracer diffusion
        msd_tr, lag = msd(dr, τ)
        Dtr = tracerD(msd_tr, 3, lag)

        # bulk diffusion
        msd_col, Nin = bulkmsd(dr, τ)
        Dbulk = bulkD(msd_col, 3, τ, Nin)

        # Haven ratio
        H = Dtr / Dbulk
    
        # --- correlation factors (ratio-of-sums over windows) ---
        # tracer:   f_tr  = (Σ_s,i ||Δr_i||^2) / (Σ_s,i Δh_i * a^2)
        # collective: f_col = (Σ_s ||Σ_i Δr_i||^2) / (Σ_s Σ_i Δh_i * a^2)
        sum_d2_tr  = 0.0; sum_dh_tr  = 0.0
        sum_d2_col = 0.0; sum_dh_col = 0.0

        @inbounds for s in 1:(S-τ)
            # tracer pieces
            for i in 1:N
                dh = acc_log[s+τ,i] - acc_log[s,i]
                if dh > 0
                    dx = dr[s+τ,1,i] - dr[s,1,i]
                    dy = dr[s+τ,2,i] - dr[s,2,i]
                    dz = dr[s+τ,3,i] - dr[s,3,i]
                    sum_d2_tr += dx*dx + dy*dy + dz*dz
                    sum_dh_tr += dh
                end
            end
            # collective pieces
            sx=0.0; sy=0.0; sz=0.0; dh_sum=0.0
            for i in 1:N
                sx += dr[s+τ,1,i] - dr[s,1,i]
                sy += dr[s+τ,2,i] - dr[s,2,i]
                sz += dr[s+τ,3,i] - dr[s,3,i]
                dh_sum += (acc_log[s+τ,i] - acc_log[s,i])
            end
            if dh_sum > 0
                sum_d2_col += sx*sx + sy*sy + sz*sz
                sum_dh_col += dh_sum
            end
        end

        ftr  = (sum_dh_tr  > 0) ? sum_d2_tr  / (sum_dh_tr  * a_lat*a_lat) : NaN
        fcol = (sum_dh_col > 0) ? sum_d2_col / (sum_dh_col * a_lat*a_lat) : NaN

        # --- reduced conductivity: C q^2 Dtr / (kB T) ---
        C = Nions / SITES
        sigma_red = C * (q*q) * Dtr / (kB * T)

        Dtrs[k]   = Dtr
        Dbulks[k] = Dbulk
        Havens[k] = H
        f_tr[k]   = ftr
        f_col[k]  = fcol
        sigma_reds[k] = sigma_red

        Printf.@printf("pct=%.2f (N=%d)  Dtr=%.6g  Dbulk=%.6g  Haven=%.6g  f_tr=%.6g  f_col=%.6g  σ_red=%.6g\n",
                       percentiles[k], Nions, Dtr, Dbulk, H, ftr, fcol, sigma_red)
    end

    open(outfile, "w") do io
        println(io, "percentile\tDtr\tDbulk\tHaven\tf_tr\tf_col\treduced_conductivity")
        for k in eachindex(percentiles)
            Printf.@printf(io, "%.2f\t%.12g\t%.12g\t%.12g\t%.12g\t%.12g\t%.12g\n",
                           percentiles[k], Dtrs[k], Dbulks[k], Havens[k], f_tr[k], f_col[k], sigma_reds[k])
        end
    end
    Printf.@printf("Wrote %s\n", outfile)

    return percentiles, Dtrs, Dbulks, Havens, f_tr, f_col, sigma_reds
end

function DefectSweepFixedN(a, b, c, Nions, sweeps, lagtime;
    defect_fracs = 0.0:0.05:0.95,
    a_lat::Float64 = 1.0,
    kB::Float64 = 1.0,
    T::Float64 = 1.0,
    q::Float64 = 1.0,
    # keep these in sync with mcloop_g! so the disorder field matches the MC energy surface:
    q_def::Float64 = +1.0,
    k_e::Float64   = 1.0,
    eps2::Float64  = 0.25,
    outfile::AbstractString = "data/defect_sweep.tsv"
)
    SITES = a * b * c
    @assert 1 <= Nions <= SITES "Nions must be between 1 and total sites."

    φs    = collect(defect_fracs)
    Mlist = round.(Int, φs .* SITES)

    Dtrs    = similar(φs, Float64)
    Dbulks  = similar(φs, Float64)
    Havens  = similar(φs, Float64)
    f_tr    = similar(φs, Float64)
    f_col   = similar(φs, Float64)
    σ_red   = similar(φs, Float64)
    χ0s     = similar(φs, Float64)   # <— allocate OUTSIDE the loop
    Svals   = similar(φs, Float64)   # <— allocate OUTSIDE the loop

    for k in eachindex(φs)
        M = Mlist[k]

        # Jammed: no vacancies ⇒ zero diffusion
        if Nions + M >= SITES
            Dtrs[k]=0.0; Dbulks[k]=0.0; Havens[k]=NaN
            f_tr[k]=NaN; f_col[k]=NaN; σ_red[k]=0.0
            χ0s[k]=NaN;  Svals[k]=NaN
            Printf.@printf("φ=%.2f (M=%d)  JAMMED ⇒ Dtr=0  Dbulk=0  Haven=NaN  f_tr=NaN  f_col=NaN  σ_red=0  χ0=NaN  S=NaN\n",
                    φs[k], M)
            continue
        end

        steps = sweeps * Nions
        dr, acc_log, pos, occ, _ = mcloop_g!(a, b, c, Nions, M, steps;
                                             q_def=q_def, k_e=k_e, a_lat=a_lat, eps2=eps2,
                                             kB=kB, T=T)

        # build the SAME potential used for MC acceptance (same params!)
        V = build_coulomb_potential(a, b, c, occ; q_def=q_def, k_e=k_e, a_lat=a_lat, eps2=eps2)

        χ0, Sval = disorder_strength(V; kB=kB, T=T, d=3)
        χ0s[k]   = χ0
        Svals[k] = Sval

        S, D, N = size(dr); @assert D == 3
        τ = lagtime; @assert 1 <= τ <= S-1

        msd_tr, lag = msd(dr, τ)
        Dtr = tracerD(msd_tr, 3, lag)

        msd_col, Nin = bulkmsd(dr, τ)
        Dbulk = bulkD(msd_col, 3, τ, Nin)

        H = Dtr / Dbulk

        # correlation factors
        sum_d2_tr = 0.0; sum_dh_tr = 0.0
        sum_d2_col = 0.0; sum_dh_col = 0.0
        @inbounds for s in 1:(S-τ)
            for i in 1:N
                dh = acc_log[s+τ,i] - acc_log[s,i]
                if dh > 0
                    dx = dr[s+τ,1,i] - dr[s,1,i]
                    dy = dr[s+τ,2,i] - dr[s,2,i]
                    dz = dr[s+τ,3,i] - dr[s,3,i]
                    sum_d2_tr += dx*dx + dy*dy + dz*dz
                    sum_dh_tr += dh
                end
            end
            sx=0.0; sy=0.0; sz=0.0; dh_sum=0.0
            for i in 1:N
                sx += dr[s+τ,1,i] - dr[s,1,i]
                sy += dr[s+τ,2,i] - dr[s,2,i]
                sz += dr[s+τ,3,i] - dr[s,3,i]
                dh_sum += (acc_log[s+τ,i] - acc_log[s,i])
            end
            if dh_sum > 0
                sum_d2_col += sx*sx + sy*sy + sz*sz
                sum_dh_col += dh_sum
            end
        end
        ftr  = (sum_dh_tr  > 0) ? sum_d2_tr  / (sum_dh_tr  * a_lat*a_lat) : NaN
        fcol = (sum_dh_col > 0) ? sum_d2_col / (sum_dh_col * a_lat*a_lat) : NaN

        # reduced conductivity
        C = Nions / SITES
        sigma_red = C * (q*q) * Dtr / (kB * T)

        Dtrs[k]   = Dtr
        Dbulks[k] = Dbulk
        Havens[k] = H
        f_tr[k]   = ftr
        f_col[k]  = fcol
        σ_red[k]  = sigma_red

        Printf.@printf("φ=%.2f (M=%d)  Dtr=%.6g  Dbulk=%.6g  Haven=%.6g  f_tr=%.6g  f_col=%.6g  σ_red=%.6g  χ0=%.6g  S=%.6g\n",
                φs[k], M, Dtr, Dbulk, H, ftr, fcol, sigma_red, χ0, Sval)
    end

    Dtr0   = Dtrs[1]
    Dbulk0 = Dbulks[1]
    norm_or_nan(x, x0) = (isfinite(x0) && x0 != 0.0) ? (x / x0) : NaN
    Dtr_norm   = [norm_or_nan(Dtrs[k],   Dtr0)   for k in eachindex(Dtrs)]
    Dbulk_norm = [norm_or_nan(Dbulks[k], Dbulk0) for k in eachindex(Dbulks)]

    mkpath(dirname(outfile))  # write wherever the caller asked
    open(outfile, "w") do io
        println(io, "defect_frac\tM\tchi0\tbeta2chi0_over3\tDtr\tDbulk\tHaven\tf_tr\tf_col\treduced_conductivity\tDtr_norm\tDbulk_norm")
        for k in eachindex(φs)
            Printf.@printf(io, "%.4f\t%d\t%.12g\t%.12g\t%.12g\t%.12g\t%.12g\t%.12g\t%.12g\t%.12g\t%.12g\t%.12g\n",
                    φs[k], Mlist[k],
                    χ0s[k], Svals[k],
                    Dtrs[k], Dbulks[k], Havens[k], f_tr[k], f_col[k], σ_red[k],
                    Dtr_norm[k], Dbulk_norm[k])
        end
    end
    Printf.@printf("Wrote %s\n", outfile)

    return φs, Mlist, Dtrs, Dbulks, Havens, f_tr, f_col, σ_red, Dtr_norm, Dbulk_norm, χ0s, Svals
end



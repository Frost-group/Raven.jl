import LinearAlgebra, IterativeSolvers, Random, Printf

# Note: this code uses multiple dispatch quite often so it is advised to read the comments carefully if you want to follow what's going on!

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
    outfile::AbstractString = "defect_sweep.tsv"
)
    SITES = a * b * c
    @assert 1 <= Nions <= SITES "Nions must be between 1 and total sites."

    φs = collect(defect_fracs)
    Mlist = round.(Int, φs .* SITES)

    Dtrs    = similar(φs, Float64)
    Dbulks  = similar(φs, Float64)
    Havens  = similar(φs, Float64)
    f_tr    = similar(φs, Float64)
    f_col   = similar(φs, Float64)
    σ_red   = similar(φs, Float64)

    for k in eachindex(φs)
        M = Mlist[k]

        # Jammed: no vacancies ⇒ zero diffusion (avoid running MC)
        if Nions + M >= SITES
            Dtrs[k]   = 0.0
            Dbulks[k] = 0.0
            Havens[k] = NaN
            f_tr[k]   = NaN
            f_col[k]  = NaN
            σ_red[k]  = 0.0
            Printf.@printf("φ=%.2f (M=%d)  JAMMED ⇒ Dtr=0  Dbulk=0  Haven=NaN  f_tr=NaN  f_col=NaN  σ_red=0\n",
                           φs[k], M)
            continue
        end

        steps = sweeps * Nions
        dr, acc_log, _, _, _ = mcloop!(a, b, c, Nions, M, steps)
        S, D, N = size(dr); @assert D == 3
        τ = lagtime; @assert 1 <= τ <= S-1

        # Diffusion coefficients
        msd_tr, lag = msd(dr, τ)
        Dtr = tracerD(msd_tr, 3, lag)

        msd_col, Nin = bulkmsd(dr, τ)
        Dbulk = bulkD(msd_col, 3, τ, Nin)

        H = Dtr / Dbulk

        # Correlation factors
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

        # Reduced conductivity: C q^2 Dtr / (kB T)
        C = Nions / SITES
        sigma_red = C * (q*q) * Dtr / (kB * T)

        Dtrs[k]   = Dtr
        Dbulks[k] = Dbulk
        Havens[k] = H
        f_tr[k]   = ftr
        f_col[k]  = fcol
        σ_red[k]  = sigma_red

        Printf.@printf("φ=%.2f (M=%d)  Dtr=%.6g  Dbulk=%.6g  Haven=%.6g  f_tr=%.6g  f_col=%.6g  σ_red=%.6g\n",
                       φs[k], M, Dtr, Dbulk, H, ftr, fcol, sigma_red)
    end
    Dtr0   = Dtrs[1]
    Dbulk0 = Dbulks[1]

    # safe normalization helpers (avoid divide-by-zero or NaN baseline)
    norm_or_nan(x, x0) = (isfinite(x0) && x0 != 0.0) ? (x / x0) : NaN

    Dtr_norm   = [norm_or_nan(Dtrs[k],   Dtr0)   for k in eachindex(Dtrs)]
    Dbulk_norm = [norm_or_nan(Dbulks[k], Dbulk0) for k in eachindex(Dbulks)]

    # --- write TSV with normalized columns added ---
    mkpath("data")
    open(outfile, "w") do io
        println(io, "defect_frac\tM\tDtr\tDbulk\tHaven\tf_tr\tf_col\treduced_conductivity\tDtr_norm\tDbulk_norm")
        for k in eachindex(φs)
            Printf.@printf(io,
                "%.4f\t%d\t%.12g\t%.12g\t%.12g\t%.12g\t%.12g\t%.12g\t%.12g\t%.12g\n",
                φs[k], Mlist[k],
                Dtrs[k], Dbulks[k], Havens[k], f_tr[k], f_col[k], σ_red[k],
                Dtr_norm[k], Dbulk_norm[k]
            )
        end
    end
    Printf.@printf("Wrote %s\n", joinpath("data", outfile))

    return φs, Mlist, Dtrs, Dbulks, Havens, f_tr, f_col, σ_red, Dtr_norm, Dbulk_norm
end


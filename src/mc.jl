import LinearAlgebra, IterativeSolvers, Random, Printf

# Version 2 initialization code based on flattened data.

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
            return true
        else
            return false
        end
    end
end

function mcloop!(a, b, c, N, steps)
    a, b, c, pos, occ, disp = initialize(a, b, c, N)
    attempts = N
    sweeps = cld(steps, attempts)
    dr_log = zeros(Int16, sweeps + 1, 3, N)
    dr_log[1, :, :] .= 0
    for sweep in 1:sweeps
        for attempt in 1:attempts
            mcstep!(a, b, c, pos, occ, disp)
            #println("attempt $(attempt + attempts * (sweep - 1))")
        end
        dr_log[sweep + 1, :, :] .= disp
        #println("sweep $sweep / $sweeps")
    end
    return dr_log, pos, occ
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

function HavenSweep(a, b, c, sweeps, lagtime, outfile="haven_sweep.tsv")
    SITES = a * b * c
    percentiles = collect(0.0 : 0.05 : 1.0)
    ioncount = round.(Int, SITES .* percentiles)
    ioncount[1] = 1

    Dtrs   = Vector{Float64}(undef, length(percentiles))
    Dbulks = Vector{Float64}(undef, length(percentiles))
    Havens = Vector{Float64}(undef, length(percentiles))

    for k in eachindex(ioncount)
        Nions = ioncount[k]
        steps = sweeps * Nions                     # same # of sweeps across loadings
        dr, _, _ = mcloop!(a, b, c, Nions, steps)  # dr :: (S,3,N)

        # tracer TAMSD → D_tr
        msd_tr, lag = msd(dr, lagtime)
        Dtr = tracerD(msd_tr, 3, lag)

        # bulk (Einstein–Helfand) TAMSD → D_bulk (a.k.a. D_sigma)
        msd_col, Nin_msdbulk = bulkmsd(dr, lagtime)      # your bulkmsd returns (value, Nions)
        Dbulk = bulkD(msd_col, 3, lagtime, Nin_msdbulk)

        H = Dbulk / Dtr

        Dtrs[k]   = Dtr
        Dbulks[k] = Dbulk
        Havens[k] = H

        Printf.@printf("pct=%.2f (N=%d)  Dtr=%.6g  Dbulk=%.6g  Haven=%.6g\n", percentiles[k], Nions, Dtr, Dbulk, H)
    end

    open(outfile, "w") do io
        println(io, "percentile\tDtr\tDbulk\tHaven")
        for k in eachindex(percentiles)
            Printf.@printf(io, "%.2f\t%.12g\t%.12g\t%.12g\n", percentiles[k], Dtrs[k], Dbulks[k], Havens[k])
        end
    end
    Printf.@printf "Wrote %s\n" outfile

    return percentiles, Dtrs, Dbulks, Havens
end

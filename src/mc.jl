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

function DtrSweep(a, b, c, sweeps, outfile="Dtr_sweep.csv")
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
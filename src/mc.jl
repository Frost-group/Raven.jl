import LinearAlgebra, IterativeSolvers, Random

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

function neighbors(a, b, c, pos, id) 
# Note that this function doesn't return the neighbors of, for example, site 7 but the ION with the ID #7!
    x = pos[id, 1]; y = pos[id, 2]; z = pos[id, 3]
    nbrs = Matrix{Int}(undef, 6, 3)
    k = 1
    for (dx, dy, dz) in ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1))
        nx = mod1(x + dx, a)
        ny = mod1(y + dy, b)
        nz = mod1(z + dz, c)
        nbrs[k, :] .= (nx, ny, nz)
        k += 1
    end
    return nbrs
end

function mcstep!(a, b, c, pos, occ, disp; verbose=false)
    ion = rand(1:length(pos[:, 1]))
    oldpos = pos[ion, :]
    nbrs = neighbors(a, b, c, pos, ion)
    newpos = nbrs[rand(1:6), :]
    if occ[newpos[1, 1], newpos[2, 1], newpos[3, 1]] == 0 # check if selected neighbor is empty
        occ[oldpos[1, 1], oldpos[2, 1], oldpos[3, 1]] = 0 # empty the old position
        occ[newpos[1, 1], newpos[2, 1], newpos[3, 1]] = ion # put the ion in the new position
        pos[ion, :] = newpos # update the #ID entry in the ion list

        # load old and new lattice positions
        x, y, z = oldpos[1, 1], oldpos[2, 1], oldpos[3, 1]
        nx, ny, nz = newpos[1, 1], newpos[2, 1], newpos[3, 1]

        # raw differences before PBC correction
        dx = nx - x
        dy = ny - y
        dz = nz - z

        # Directly asking whether we stayed on the same coordinate, the next cell, or the previous one?
        dx = (nx == x) ? 0 : (nx == mod1(x + 1, a) ? 1 : -1)
        dy = (ny == y) ? 0 : (ny == mod1(y + 1, b) ? 1 : -1)
        dz = (nz == z) ? 0 : (nz == mod1(z + 1, c) ? 1 : -1)

        # adding each displacement to dr
        disp[1, ion] += dx
        disp[2, ion] += dy
        disp[3, ion] += dz

        #if verbose
            #println("dx:$dx dy:$dy dz:$dz")
        #end

        return true
    else
        return false
    end
end

function mcloop!(a, b, c, N, steps)
    a, b, c, pos, occ, disp = initialize(a, b, c, N)
    attempts = N
    sweeps = cld(steps, attempts)
    dr_log = zeros(Int, sweeps + 1, 3, N)
    dr_log[1, :, :] .= 0
    for sweep in 1:sweeps
        for attempt in 1:attempts
            mcstep!(a, b, c, pos, occ, disp, verbose=true)
            #println("attempt $(attempt + attempts * (sweep - 1))")
        end
        dr_log[sweep + 1, :, :] .= disp
        #println("sweep $sweep / $sweeps")
    end
    return dr_log, pos, occ
end

function msd(dr, lag)
    N, _, ioncount = size(dr) # number of data points
    n = lag # comparing which steps?
    msd = 1 / ((N - n) * ioncount) * sum(LinearAlgebra.norm(dr[i + n, :, j]- dr[i, :, j])^2 for i in 1:N-n for j in 1:ioncount)
    return msd
end

function tracerD(msd, dim, sweeps)
    Dtr = msd / (2 * dim * sweeps)
    return Dtr
end

function DtrSweep(a, b, c, steps)
    N = a * b * c
    for i in 1:N
        dr, _, _ = mcloop!(a, b, c, i, steps)
        msd = Raven.msd(dr, 1)
        sweeps = round(steps/i)
        Dtr = tracerD(msd, 3, sweeps)
        print("$Dtr \n")
    end
end

# we might need to get rid of the flattening afterall. RIP.
"""
function flatten(X, Y, Z, x, y, z) # for flattening 3D coordinates to indexing i.
    return x + yX + zXY
end

function unflatten(X, Y, Z, i) # for unflattening index i back to 3D coordinates.
    x = i % X
    y = (i ÷ X) % Y
    z = i ÷ (X * Y) 
    return (x, y, z)
end

function checkneighbors(A) # selecting a random site from the ion list and checking if its neighbors are occupied or not.
    id = rand(ionlist)
    I = findfirst(cell -> (id in cell, A))
    neighbors = []
    push![(I[1] ± 1, I[2], I[3]), neighbors]
    push![(I[1], I[2] ± 1, I[3]), neighbors]
    push![(I[1], I[2], I[3] ± 1), neighbors]

function listneighbors(id, ionlist)
    target = reshape(ionlist[id, 2:4], 1, 3) # returns [x, y, z] for ion number id. Update: when we extract xyz, it returns as a 3x1 vector instead of a matrix?
    # The reshape function converts a 3x1 vector into a 1x3 matrix.    
    deltas = [ 1  0  0;
              -1  0  0;
               0  1  0;
               0 -1  0; 
               0  0  1;
               0  0 -1 ] # sooo much cleverer using a delta -- gotta thank Jarv!
    neighbors = target .+ deltas
    return neighbors
end

function flatten(a, b, c, x, y, z) # flattening 3D coordinates to a linear index coordinate.
    LinearIndices((a, b, c))[x, y, z]
end

function unflatten(s, a, b, c)
    Tuple(CartesianIndices((a, b, c))[s])
end
"""
# Version 1 demo before implementing data flattening below.
"""
function initializesystem(a,b,c,f) # creating a physical system of occupancy(0s and 1s, expand later to +1s and -1s). x,y,z for dimensions and O for the occupancy fraction.
    N = round(a * b * c * f)
    N = Int64(N)
    lattice = zeros(Int64, a, b, c)
    ionlist = zeros(Int64, N, 4)
    id=1
    while id - 1 < N
        x, y, z = rand(1:a), rand(1:b), rand(1:c)
        lattice[x, y, z] = id 
        ionlist[id, :] = [id, x, y, z]
        id = id + 1
    end
    return lattice, ionlist
"""
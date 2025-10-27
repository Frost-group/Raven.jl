import LinearAlgebra, IterativeSolvers, Random

# Version 2 initialization code based on flattened data.

function initialize(a, b, c, f)
    S = a * b * c # S is the already 'flattened' 1D index that can be linked 1 to 1 to the original 3D lattice via the 'unflatten' function.
    N = round(a * b * c * f) # number of ions to be generated within the flattened lattice.
    N = Int(N)
    d0 = fill((0.0, 0.0, 0.0), N) # initialization of positions of each ion per row. later used for MSD computation.
    dr = fill((0.0, 0.0, 0.0), N) # relative position tracker per sweep, also later used for MSD computation.
    sel = Random.randperm(S)[1:N] # shuffle flattened matrix S to select the first N entries to determine where ions reside in.
    pos = Vector{Int}(undef, N)
    occ = fill(0, S)
    for id in 1:N
        pos[id] = sel[id]
        occ[sel[id]] = id
    end
    return a, b, c, pos, occ, d0
end   

function flatten(a, b, c, x, y, z) # flattening 3D coordinates to a linear index coordinate.
    LinearIndices((a, b, c))[x, y, z]
end

function unflatten(s, a, b, c)
    Tuple(CartesianIndices((a, b, c))[s])
end

function neighbors(a, b, c, pos, id) 
# Note that this function doesn't return the neighbors of, for example, site 7 but the ION with the ID #7!
    target = pos[id]
    nbrs = Vector{Int}(undef, 6)
    (x, y, z) = unflatten(target, a, b, c)
    k = 1
    for (dx, dy, dz) in ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1))
        nx = mod1(x + dx, a)
        ny = mod1(y + dy, b)
        nz = mod1(z + dz, c)
        nbrs[k] = flatten(a, b, c, nx, ny, nz)
        k += 1
    end
    return nbrs
end

function mcstep(a, b, c, pos, occ, dr)
    stepid = rand(1:length(pos))
    oldpos = pos[stepid]
    nbrs = Raven.neighbors(a, b, c, pos, stepid)
    newpos = rand(nbrs)
    if occ[newpos] == 0 # check if selected neighbor is empty
        occ[oldpos] = 0 # empty the old position
        occ[newpos] = stepid # put the ion in the new position
        pos[stepid] = newpos # update the #ID entry in the ion list

        # load old and new lattice positions
        x, y, z = unflatten(oldpos, a, b, c)
        nx, ny, nz = unflatten(newpos, a, b, c)

        # raw differences before PBC correction
        dx = nx - x
        dy = ny - y
        dz = nz - z

        # PBC handling: when crossing the +x boundary, a -> 1 should be +1 and not 1-a
        if dx == (a - 1); dx = -1 elseif dx == 1 - a; dx = 1 end
        if dy == (b - 1); dy = -1 elseif dy == 1 - b; dy = 1 end
        if dz == (c - 1); dz = -1 elseif dz == 1 - c; dz = 1 end
    
        return true
    else
        return false
    end
end

function mcloop(a, b, c, f, steps)
    a, b, c, pos, occ. d0 = Raven.initialize(a, b, c, f)
    dr = 
    for i in 1:steps
        mcstep(a, b, c, pos, occ)
        print("step $i is complete. \n")
    end
    return dr
end

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
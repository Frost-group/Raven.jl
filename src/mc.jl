import LinearAlgebra, IterativeSolvers, Random

# Version 2 initialization code based on flattened data.

function initialize(a, b, c, f)
    S = a * b * c # S is the already 'flattened' 1D index that can be linked 1 to 1 to the original 3D lattice via the 'unflatten' function.
    N = round(a * b * c * f) # number of ions to be generated within the flattened lattice.
    N = Int64(N)
    sel = Random.randperm(S)[1:N] # shuffle flattened matrix S to select the first N entries to determine where ions reside in.
    pos = Vector{Int}(undef, N)
    occ = fill(0, S)
    for id in 1:N
        pos[id] = sel[id]
        occ[sel[id]] = id
    end
    return a, b, c, pos, occ
end   

function neighbors(id, a, b)
    target = pos[id]
    deltas = [1; -1; a; -a; a * b; -a * b] # defining deltas to look at neighbors in linear representation.
    neighbors = target .+ deltas
    return neighbors
end

function mcstep(lattice, ionlist)
    movingionid=rand(ionlist[:, 1])
    oldpos = reshape(ionlist[movingionid, 2:4], 1, 3)
    candidates = Raven.listneighbors(movingionid, ionlist)
    newpos = reshape(candidates[rand(1:size(candidates)[1]), :], 1, 3)
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
end
"""
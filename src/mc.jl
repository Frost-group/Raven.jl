import LinearAlgebra, IterativeSolvers

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
function checkneighbors(A) # selecting a random site from the ion list and checking if its neighbors are occupied or not.
    id = rand(ionlist)
    I = findfirst(cell -> (id in cell, A))
    neighbors = []
    push![(I[1] ± 1, I[2], I[3]), neighbors]
    push![(I[1], I[2] ± 1, I[3]), neighbors]
    push![(I[1], I[2], I[3] ± 1), neighbors]
"""
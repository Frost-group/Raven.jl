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

"""
function checkneighbors(A) # selecting a random site from the ion list and checking if its neighbors are occupied or not.
    id = rand(ionlist)
    I = findfirst(cell -> (id in cell, A))
    neighbors = []
    push![(I[1] ± 1, I[2], I[3]), neighbors]
    push![(I[1], I[2] ± 1, I[3]), neighbors]
    push![(I[1], I[2], I[3] ± 1), neighbors]
"""
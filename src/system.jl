import LinearAlgebra, IterativeSolvers

function initializesystem(a,b,c,O) # creating a physical system of occupancy(0s and 1s, expand later to +1s and -1s). x,y,z for dimensions and O for the occupancy fraction.
    N = round(a * b * c * O)
    A = [Vector{Int}() for _ in 1:a, _ in 1:b, _ in 1:c]
    ionlist=[]
    id=1
    placed=0
    while placed < N
        x, y, z = rand(1:a), rand(1:b), rand(1:c)
        if isempty(A[x, y, z])
            push!(A[x, y, z], id)
            push!(ionlist, id)
            id = id + 1
            placed = placed + 1
        end
    end
    return A, ionlist
end